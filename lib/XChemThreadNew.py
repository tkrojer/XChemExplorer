import os, sys, glob
from datetime import datetime
from PyQt4 import QtGui, QtCore

import time
import pickle
import cPickle
import base64
import math
import subprocess
from datetime import datetime
import time
import getpass
import csv

sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
from XChemUtils import process
from XChemUtils import parse
from XChemUtils import queue
from XChemUtils import mtztools
from XChemUtils import pdbtools
from XChemUtils import helpers
from XChemUtils import reference
from XChemUtils import misc
import XChemDB
import XChemLog
import XChemMain


class choose_autoprocessing_outcome(QtCore.QThread):
    def __init__(self,
                database,
                visit,
                referenceFileList,
                preferences,
                projectDir,
                rescore,
                xce_logfile):
        QtCore.QThread.__init__(self)
        self.visit = visit
        self.projectDir = projectDir
        self.reference_file_list = reference_file_list
        self.selection_mechanism = preferences['dataset_selection_mechanism']

        self.rescore = rescore
        self.acceptable_unitcell_volume_difference = acceptable_unitcell_volume_difference
        self.acceptable_low_resolution_limit_for_data = acceptable_low_resolution_limit_for_data
        self.acceptable_low_resolution_Rmerge = acceptable_low_resolution_Rmerge

        self.xce_logfile = xce_logfile
        self.Logfile = XChemLog.updateLog(xce_logfile)

        self.db = XChemDB.data_source(os.path.join(database))
        self.allSamples = self.db.collected_xtals_during_visit_for_scoring(visitID,rescore)


    def run(self):
        for sample in sorted(self.allSamples):
            self.Logfile.insert('%s: selecting autoprocessing result' %sample)
            dbList = self.db.all_autoprocessing_results_for_xtal_as_dict(sample)
            self.Logfile.insert('%s: found %s different autoprocessing results' %(sample,str(len(dbList))))
            # 1.) if posssible, only carry forward samples with similar UCvolume and same point group
            dbList = self.selectResultsSimilarToReference(dbList)
            # 2.) if possible, only carry forward samples with low resolution Rmerge smaller than
            #     specified in perferences
            dbList = self.selectResultsWithAcceptableLowResoRmerge(dbList)
            # 3.) Make selection based on speified selection mechanism
            if self.selection_mechanism == 'IsigI*Comp*UniqueRefl':
                dbDict = selectHighestScore(dbList)
            elif self.selection_mechanism == 'highest_resolution':
                dbDict = selectHighestResolution(dbList)
            # 4.) Set new symbolic links in project directory
            self.linkBestAutoProcessingResult(dbDict)
            # 5.) update database
            dbDict('DataProcessingAutoAssigned'] = 'True'
            self.update_data_source(sample,dbDict)

    def selectResultsSimilarToReference(self,dbList):
        dbListOut = []
        for resultDict in dbList:
            try:
                if isinstance(float(resultDict['DataProcessingUnitCellVolume']),float):
                    for reference_file in self.reference_file_list:
                        if not reference_file[4]==0:
                            unitcell_difference=round((math.fabs(reference_file[4]-float(resultDict['DataProcessingUnitCellVolume']))/reference_file[4])*100,1)
                            if unitcell_difference < self.acceptable_unitcell_volume_difference and reference_file[3]==resultDict['DataProcessingLattice']:
                                dbListOut.append(resultDict)
            except ValueError:
                pass
        if dbListOut == []:
            dbListOut = dbList
        return dbListOut


    def selectResultsWithAcceptableLowResoRmerge(self,dbList):
        dbListOut = []
        for resultDict in dbList:
            try:
                if float(resultDict['DataProcessingRmergeLow']) < self.acceptable_low_resolution_Rmerge:
                    dbListOut.append(resultDict)
            except ValueError:
                pass
        if dbListOut == []:
            dbListOut = dbList
        return dbListOut


    def selectHighestScore(self,dbList):
        dbListOut = []
        tmp = []
        for resultDict in dbList:
            try:
                tmp.append(float(resultDict['DataProcessingScore'])
            except ValueError:
                tmp.append(0.0)
        highestScoreDict = dbList.index(max(tmp))
        return highestScoreDict

    def selectHighestResolution(self,dbList):
        dbListOut = []
        tmp = []
        for resultDict in dbList:
            try:
                tmp.append(float(resultDict['DataProcessingResolutionHigh'])
            except ValueError:
                tmp.append(100.0)
        highestResoDict = dbList.index(min(tmp))
        return highestResoDict


    def linkBestAutoProcessingResult(self,xtal,dbDict):
        run =      dbDict['DataCollectionRun']
        visit =    dbDict['DataCollectionVisit']
        autoproc = dbDict['DataProcessingProgram']
        mtzFileAbs = dbDict['DataProcessingPathToMTZfile']
        mtzfile = mtzFileAbs[mtzFileAbs.rfind('/')+1:]
        logFileAbs = dbDict['DataProcessingPathToLOGfile']
        logfile = logFileAbs[logFileAbs.rfind('/')+1:]

        os.chdir(os.path.join(self.projectDir,xtal)
        # MTZ file
        if os.path.isfile(xtal+'.mtz'):
            os.system('/bin/rm %s.mtz' %xtal)
        if os.path.isfile(os.path.join('autoprocessing', visit + '-' + run + autoproc, mtzfile)):
            os.symlink(os.path.join('autoprocessing', visit + '-' + run + autoproc, mtzfile), xtal + '.mtz')
        # LOG file
        if os.path.isfile(xtal+'.log'):
            os.system('/bin/rm %s.log' %xtal)
        if os.path.isfile(os.path.join('autoprocessing', visit + '-' + run + autoproc, logfile)):
            os.symlink(os.path.join('autoprocessing', visit + '-' + run + autoproc, logfile), xtal + '.log')
       





class read_write_autoprocessing_results_from_to_disc(QtCore.QThread):

    """
    major changes:
    - copying of files and updating of DB will be combined in one class
    - fast_dp output will be ignored
    - pkl file will be abandoned
    - results for every autoprocessing output will be recorded in new DB table
    - crystal alignment files will be copied to project directory
    - beamline directory in project directory will not be used anymore; user needs to set data collection dir
    - DB mainTable gets flag if user updated autoprocessing selection
    - checking of reprocessed files needs to be explicit
    Note:
    - still miss reading of pins
    - still miss copying of image (snapshots)

    """
    def __init__(self,
                processedDir,
                database,
                projectDir,
                xce_logfile):
        QtCore.QThread.__init__(self)
        self.processedDir =  processedDir
        self.visit,self.beamline = self.getVisit()
        self.projectDir = projectDir
        self.Logfile = XChemLog.updateLog(xce_logfile)

        self.db = XChemDB.data_source(os.path.join(database))
        self.exisitingSamples = getExistingSamples()

        self.toParse = [
                [   os.path.join(run, 'xia2', '*'),
                    os.path.join('LogFiles', '*aimless.log'),
                    os.path.join('DataFiles', '*free.mtz')],
                [   os.path.join(run, 'autoPROC', '*')
                    '*aimless.log',
                    '*truncate-unique.mtz']
                        ]

    def run(self):
        self.parse_file_system()

    def getExistingSamples(self):
        existingSamples={}
        self.Logfile.insert('reading existing samples from collectionTable')
        allEntries = self.db.execute_statement('select CrystalName,VisitRunAutoproc from collectionTable')
        for item in allEntries:
            if item[0] not in existingSamples:
                existingSamples[item[0]]=[]
            self.Logfile.insert('%s: adding %s' %(item[0],item[1]))
            existingSamples[item[0]].append(item[1])
        return existingSamples

    def getVisit(self):
        visit = 'unknown'
        beamline = 'unknown'
        if 'attic' in self.processedDir:
            try:
                visit=visit_directory.split('/')[6]
                beamline=visit_directory.split('/')[3]
            except IndexError:
                pass
        else:
            try:
                visit=visit_directory.split('/')[5]
                beamline=visit_directory.split('/')[2]
            except IndexError:
                pass
        return visit,beamline


    def createAutoprocessingDir(self,xtal,run,autoproc):
        # create all the directories if necessary
        if not os.path.isdir(os.path.join(self.projectDir,xtal)):
            os.mkdir(os.path.join(self.projectDir,xtal))
        if not os.path.isdir(os.path.join(self.projectDir,xtal,'autoprocessing')):
            os.mkdir(os.path.join(self.projectDir,xtal,'autoprocessing'))
        if not os.path.isdir(
                os.path.join(self.projectDir,xtal,'autoprocessing', self.visit + '-' + run + autoproc)):
            os.mkdir(os.path.join(self.projectDir,xtal,'autoprocessing', self.visit + '-' + run + autoproc))

    def copyMTZandLOGfiles(self,xtal,run,autoproc,mtzfile,logfile):
        mtzNew = ''
        logNew = ''
        os.chdir(os.path.join(self.projectDir,xtal,'autoprocessing', self.visit + '-' + run + autoproc))
        self.Logfile.insert('%s: copying %s' %(xtal,mtzfile))
        os.system('/bin/cp ' + mtzfile + ' .')
        if os.path.isfile(mtzfile[mtzfile.rfind('/')+1:]):
            os.symlink(mtzfile[mtzfile.rfind('/')+1:], xtal + '.mtz')
            mtzNew=os.path.join(self.projectDir,xtal,'autoprocessing', self.visit + '-' + run + autoproc, xtal + '.mtz')
        self.Logfile.insert('%s: copying %s' %(xtal,logfile))
        os.system('/bin/cp ' + logfile + ' .')
        if os.path.isfile(logfile[logfile.rfind('/')+1:]):
            os.symlink(logfile[logfile.rfind('/')+1:], xtal + '.log')
            logNew=os.path.join(self.projectDir,xtal,'autoprocessing', self.visit + '-' + run + autoproc, xtal + '.log')
        return mtzNew,logNew

    def makeJPGdir(self,xtal):
        if not os.path.isdir(os.path.join(self.projectDir,xtal,'jpg')):
            self.Logfile.insert('making jpg directory in '+xtal)
            os.mkdir(os.path.join(self.projectDir,xtal,'jpg'))



    def readProcessingResults(self,xtal,folder,log,mtz,timestamp,current_run):
        db_dict = {}
        if 'ap-run' in folder:
            autoproc = 'autoPROC'
        else:
            autoproc = folder.split('/')[len(folder.split('/'))-1]

        for mtzfile in glob.glob(os.path.join(folder,mtz)):

            if self.visit + '-' + current_run + autoproc in self.exisitingSamples[xtal]:
                self.Logfile.insert(
                    '%s: results from %s already parsed; skipping...' % (xtal, visit + '-' + current_run + autoproc))
                continue

            for logfile in glob.glob(os.path.join(folder, log)):
                self.createAutoprocessingDir(self, xtal, current_run, autoproc)
                mtzNew,logNew = self.copyMTZandLOGfiles(xtal,run,autoproc,mtzfile,logfile)
                db_dict = { 'DataCollectionDate': timestamp,
                            'DataProcessingPathToLogfile': logNew,
                            'DataProcessingPathToMTZfile': mtzNew,
                            'DataProcessingDirectoryOriginal': folder    }
                db_dict.update(parse().read_aimless_logfile(logNew))

        return autoproc,db_dict

    def update_data_collection_table(self,xtal,current_run,autoproc,db_dict):
        condition_dict = {  'CrystalName':  xtal,
                            'DataCollectionVisit':  self.visit,
                            'DataCollectionRun':    current_run,
                            'DataProcessingPipeline':   autoproc    }

        self.db.update_insert_any_table(self, 'collectionTable', db_dict, condition_dict)


    def getProgressSteps(self):
        if len(glob.glob(os.path.join(self.processedDir,'*'))) == 0:
            progress_step = 1
        else:
            progress_step = 100 / float(len(glob.glob(os.path.join(self.processedDir,'*'))))
        return progress_step

    def parse_file_system(self):
        self.Logfile.insert('checking for new data processing results in '+self.processedDir)
        progress = 0
        progress_step = self.getProgressSteps()

        for collected_xtals in sorted(glob.glob(os.path.join(self.processedDir,'*'))):
            if 'tmp' in collected_xtals or 'results' in collected_xtals or 'scre' in collected_xtals:
                continue

            xtal = collected_xtals[collected_xtals.rfind('/')+1:]

            # create directory for crystal aligment images in projectDir
            self.makeJPGdir(xtal)

            for run in sorted(glob.glob(os.path.join(collected_xtals,'*'))):
                current_run=run[run.rfind('/')+1:]
                timestamp=datetime.fromtimestamp(os.path.getmtime(run)).strftime('%Y-%m-%d %H:%M:%S')

                for item in self.toParse:
                    procDir = item[0]
                    logfile = item[1]
                    mtzfile = item[2]
                    for folder in glob.glob(procDir):
                        autoproc,db_dict = self.readProcessingResults(xtal,folder,logfile,mtzfile,timestamp,current_run)
                        if db_dict != {}:
                            self.update_data_collection_table(xtal,current_run,autoproc,db_dict)

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

#        self.emit(QtCore.SIGNAL('create_widgets_for_autoprocessing_results_only'), self.data_collection_dict)
        self.emit(QtCore.SIGNAL("finished()"))
