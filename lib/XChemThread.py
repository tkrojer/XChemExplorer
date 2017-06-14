# last edited: 13/06/2017, 15:00

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


class update_datasource_from_file_system(QtCore.QThread):
    def __init__(self,initial_model_directory,datasource,panddas_directory,xce_logfile):
        QtCore.QThread.__init__(self)
        self.initial_model_directory=initial_model_directory
        self.datasource=datasource
        self.db=XChemDB.data_source(self.datasource)
        self.panddas_directory=panddas_directory
        self.Logfile=XChemLog.updateLog(xce_logfile)

    def run(self):
        self.Logfile.insert('new project directory: '+self.initial_model_directory)
        self.Logfile.insert('updating data source from file system')
        progress_step=1
        if len(glob.glob(os.path.join(self.initial_model_directory,'*'))) != 0:
            progress_step=100/float(len(glob.glob(os.path.join(self.initial_model_directory,'*'))))
        else:
            progress_step=1
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        all_samples_in_datasource=self.db.get_all_samples_in_data_source_as_list()

        for directory in sorted(glob.glob(os.path.join(self.initial_model_directory,'*'))):
            try:
                os.chdir(directory)
            except OSError:
                # this could happen if the user accidentaly left a file in the project directory
                continue
            xtal=directory[directory.rfind('/')+1:]
            if xtal not in all_samples_in_datasource:
                self.Logfile.insert('inserting '+xtal+' into data source')
                self.db.execute_statement("insert into mainTable (CrystalName) values ('%s');" %xtal)
                all_samples_in_datasource.append(xtal)
            compoundID=str(self.db.get_value_from_field(xtal,'CompoundCode')[0])
            db_dict={}
            sample_dict=self.db.get_db_dict_for_sample(xtal)

            dimple_path=''  # will be set to correct path if dimple.pdb is present;
            if os.path.isfile('dimple.pdb'):
                print ''
#                if not os.path.isfile('refine.pdb'):
#                    os.system('/bin/rm refine.pdb')         # this removes broken links that could trip the symlink
#                    os.symlink('dimple.pdb', 'refine.pdb')
            else:
                os.system('/bin/rm dimple.pdb 2> /dev/null')    # this makes sure that any broken link which could rerail PANDDA gets removed
            if os.path.isfile('dimple.mtz'):
                db_dict['DimplePathToMTZ']=os.path.realpath(os.path.join(directory,'dimple.mtz'))
                dimple_mtz=db_dict['DimplePathToMTZ']
                dimple_path=dimple_mtz[:dimple_mtz.rfind('/')]
                if not os.path.isfile('refine.mtz'):
                    os.system('/bin/rm refine.mtz')
                    os.symlink('dimple.mtz', 'refine.mtz')
            else:
                os.system('/bin/rm dimple.mtz 2> /dev/null')    # this makes sure that any broken link which could rerail PANDDA gets removed
            # this should not really be the case, but if a user does not provide an aimless logfile then that's all we can do
            if not os.path.isfile(xtal+'.log'):
                if os.path.isfile(xtal+'.mtz'):
                    mtz_info=mtztools(xtal+'.mtz').get_information_for_datasource()
                    db_dict.update(mtz_info)
                    db_dict['DataCollectionOutcome']='success'
            if not os.path.isfile(xtal+'.free.mtz'):
                os.system('/bin/rm '+xtal+'.free.mtz 2> /dev/null')      # remove possible broken link
                if os.path.isfile(os.path.join(dimple_path,'prepared2.mtz')):
                    os.symlink(os.path.join(dimple_path,'prepared2.mtz'),xtal+'.free.mtz')
                    db_dict['RefinementMTZfree']=xtal+'.free.mtz'
                elif os.path.isfile(os.path.join(dimple_path,'prepared.mtz')):
                    os.symlink(os.path.join(dimple_path,'prepared.mtz'),xtal+'.free.mtz')
                    db_dict['RefinementMTZfree']=xtal+'.free.mtz'
            if os.path.isfile('refine.mtz'):
                if sample_dict['RefinementOutcome']=='None' or sample_dict['RefinementOutcome']=='':
                    refinement_in_progress=False
                    for dirs in glob.glob('*'):
                        if dirs.startswith('Refine_') and os.path.isdir(dirs):
                            db_dict['RefinementOutcome']='3 - In Refinement'
                            refinement_in_progress=True
                            break
                    if not refinement_in_progress:
                        db_dict['RefinementOutcome']='1 - Analysis Pending'
            if os.path.isdir('compound'):
                if sample_dict['CompoundCode']=='None' or sample_dict['CompoundCode']=='':
                    for smiles in glob.glob('compound/*'):
                        if smiles.endswith('smiles'):
                            for line in open(smiles):
                                if len(line.split()) >= 1:
                                    db_dict['CompoundSMILES']=line.split()[0]
                                    db_dict['CompoundCode']=smiles[smiles.rfind('/')+1:smiles.rfind('.')]
                                    compoundID=db_dict['CompoundCode']
                                    break

            for file in glob.glob('*'):
                if file==xtal+'.log' and os.path.isfile(file):
                    db_dict['DataProcessingPathToLogfile']=os.path.join(directory,xtal+'.log')
                    db_dict['DataProcessingLOGfileName']=xtal+'.log'
                    if sample_dict['DataCollectionOutcome']=='None' or sample_dict['DataCollectionOutcome']=='':
                        db_dict['DataCollectionOutcome']='success'
                    aimless_results=parse().read_aimless_logfile(file)
                    db_dict.update(aimless_results)
                if file==xtal+'.mtz' and os.path.isfile(file):
                    db_dict['DataProcessingPathToMTZfile']=os.path.join(directory,xtal+'.mtz')
                    db_dict['DataProcessingMTZfileName']=xtal+'.mtz'
                if file==xtal+'.free.mtz' and os.path.isfile(file):
                    db_dict['RefinementMTZfree']=os.path.join(directory,xtal+'.free.mtz')
                if file==compoundID+'.cif' and os.path.isfile(file):
                    db_dict['RefinementCIF']=os.path.join(directory,compoundID+'.cif')
                if file=='dimple.pdb' and os.path.isfile(file):
                    db_dict['DimplePathToPDB']=os.path.realpath(os.path.join(directory,'dimple.pdb'))
                    pdb_info=parse().PDBheader(file)
                    db_dict['DimpleRcryst']=pdb_info['Rcryst']
                    db_dict['DimpleRfree']=pdb_info['Rfree']
                    db_dict['DimpleResolutionHigh']=pdb_info['ResolutionHigh']
                if file=='refine.pdb' and os.path.isfile(file):
#                    if sample_dict['RefinementOutcome']=='None' or sample_dict['RefinementOutcome']=='':
#                        db_dict['RefinementOutcome']='3 - In Refinement'
                    db_dict['RefinementPDB_latest']=os.path.realpath(os.path.join(directory,'refine.pdb'))
                    pdb_info=parse().PDBheader(file)
                    db_dict['RefinementRcryst']=pdb_info['Rcryst']
                    db_dict['RefinementRfree']=pdb_info['Rfree']
                    db_dict['RefinementSpaceGroup']=pdb_info['SpaceGroup']
                    db_dict['RefinementRmsdBonds']=pdb_info['rmsdBonds']
                    db_dict['RefinementRmsdAngles']=pdb_info['rmsdAngles']
                if file=='refine.mtz' and os.path.isfile(file):
                    db_dict['RefinementMTZ_latest']=os.path.realpath(os.path.join(directory,'refine.mtz'))

            if db_dict != {}:
                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'updating datasource for '+xtal)
                self.db.update_data_source(xtal,db_dict)

            # also need to update PANDDA table...
            pandda_models=self.db.execute_statement("select CrystalName,PANDDA_site_index,PANDDA_site_spider_plot,PANDDA_site_event_map from panddaTable where CrystalName='%s'" %xtal)
            if not pandda_models == []:
                for entry in pandda_models:
                    db_pandda_dict={}
                    db_pandda_dict['PANDDA_site_index']=entry[1]
                    db_pandda_dict['PANDDApath']=self.panddas_directory
                    if entry[3] != None:
                        event_map=os.path.join(self.initial_model_directory,xtal,entry[3].split('/')[len(entry[3].split('/'))-1])
                        if os.path.isfile(event_map):
                            db_pandda_dict['PANDDA_site_event_map']=event_map
                    if entry[2] != None:
                        spider_plot=os.path.join(self.initial_model_directory,xtal,entry[2].split('/')[len(entry[2].split('/'))-3],entry[2].split('/')[len(entry[2].split('/'))-2],entry[2].split('/')[len(entry[2].split('/'))-1])
                        if os.path.isfile(spider_plot):
                            db_pandda_dict['PANDDA_site_spider_plot']=spider_plot
#                            db_pandda_dict['RefinementOutcome']='3 - In Refinement'    # just in case; presence of a spider plot definitely signals that refinement happened
                                                                                        # should probably be not updated! Will overwrite CompChem ready
                    self.Logfile.insert('updating panddaTable for xtal: %s, site: %s' %(entry[0],entry[1]))
                    self.db.update_insert_panddaTable(xtal,db_pandda_dict)

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.Logfile.insert('datasource update finished')


class synchronise_db_and_filesystem(QtCore.QThread):

    '''
        - remove broken links
        - insert new samples in DB
        - update data for existing samples
    '''

    def __init__(self,initial_model_directory,datasource,panddas_directory,xce_logfile,mode):
        QtCore.QThread.__init__(self)
        self.initial_model_directory=initial_model_directory
        self.datasource=datasource
        self.db=XChemDB.data_source(self.datasource)
        self.all_samples_in_datasource=self.db.get_all_samples_in_data_source_as_list()
        self.panddas_directory=panddas_directory
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.mode=mode

    def run(self):
        self.Logfile.insert('synchronising database and filesystem')
        self.Logfile.insert('current project directory: '+self.initial_model_directory)

        #
        # get list of xtals
        #

        self.xtal_list = []
        progress_step=1
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
        if self.mode != 'project_directory':
            # if only a single xtal is synched, then self.mode==xtalID
            self.Logfile.insert('synchronising '+self.mode+' only')
            self.xtal_list.append(self.mode)
        else:
            if len(glob.glob(os.path.join(self.initial_model_directory,'*'))) != 0:
                progress_step=100/float(len(glob.glob(os.path.join(self.initial_model_directory,'*'))))
            self.Logfile.insert('found '+str(len(glob.glob(os.path.join(self.initial_model_directory,'*'))))+' samples in project directory')
            for directory in sorted(glob.glob(os.path.join(self.initial_model_directory,'*'))):
                try:
                    os.chdir(directory)
                except OSError:
                    # this could happen if the user accidentaly left a file in the project directory
                    continue
                xtal=directory[directory.rfind('/')+1:]
                self.xtal_list.append(xtal)

        #
        # go through list
        #

        for xtal in self.xtal_list:
            self.Logfile.insert('directory name: '+xtal+' = sampleID in database')
            os.chdir(os.path.join(self.initial_model_directory,xtal))
            if xtal not in self.all_samples_in_datasource:
                self.Logfile.insert('sampleID not found in database: inserting '+xtal)
                self.db.execute_statement("insert into mainTable (CrystalName) values ('%s');" %xtal)
                self.all_samples_in_datasource.append(xtal)

            db_dict=self.db.get_db_dict_for_sample(xtal)

            db_dict['ProjectDirectory'] = self.initial_model_directory

            db_dict=self.sync_data_processing(xtal,db_dict)

            db_dict=self.sync_dimple_results(xtal,db_dict)

            db_dict=self.sync_compound_information(db_dict)

            db_dict=self.sync_refinement_results(xtal,db_dict)

            if db_dict != {}:
                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'updating datasource for '+xtal)
                self.db.update_data_source(xtal,db_dict)

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.Logfile.insert('database mainTable update finished')
        self.Logfile.insert('updating panddaTable')
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'updating panddaTable')
#        self.sync_pandda_table()
        self.sync_pandda_table_NEW()
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'database panddaTable update finished')
        self.Logfile.insert('database panddaTable update finished')

        self.emit(QtCore.SIGNAL('datasource_menu_reload_samples'))


    def change_absolute_to_relative_links(self,target,filename):
        os.unlink(filename)
        os.symlink(os.path.relpath(target),filename)
        self.Logfile.insert('%s -> %s' %(os.readlink(filename),target))

    def find_file(self,filename,xtal):
        found_file=False
        for files in glob.glob('*'):
            if files == filename:
                # at this point, this could ba a file or a (broken) symbolic link
                try:
                    link=os.readlink(filename)
                    if link.startswith('/'):
                        if os.path.isfile(link.replace(link[:link.find(xtal)],self.initial_model_directory+'/')):
                            target=link.replace(link[:link.find(xtal)],self.initial_model_directory+'/')
                            found_file=True
                            self.change_absolute_to_relative_links(target,filename)
                        elif os.path.isfile(files):
                            # this will leave the absolute path
                            found_file=True
                        else:
                            self.Logfile.insert('removing broken link: '+filename)
                            os.system('/bin/rm %s 2> /dev/null' %filename)
                    else:
                        if not os.path.isfile(link):
                            self.Logfile.insert('removing broken link: '+filename)
                            os.system('/bin/rm %s 2> /dev/null' %filename)
                        else:
                            found_file=True
                except OSError:
                    if os.path.isfile(filename):
                        found_file=True
                break
        return found_file




    def sync_data_processing(self,xtal,db_dict):

        # AIMLESS logfile

        # in case the MTZ file which is used for refinement is different to the one used for refinement
        if os.path.isfile('refine.mtz'):
            if os.path.isfile(xtal+'.free.mtz'):
                freeMTZ=mtztools(xtal+'.free.mtz')
                nREFfree=freeMTZ.get_number_measured_reflections()
                if os.path.isfile(xtal+'.mtz'):
                    procMTZ=mtztools(xtal+'.mtz')
                    nREF=procMTZ.get_number_measured_reflections()
                    CC,errorMessage=freeMTZ.calculate_correlaton_between_intensities_in_mtzfiles(xtal+'.mtz')
                    self.Logfile.insert('%s: calculating CC between %s.free.mtz (%s refl) and %s.mtz (%s refl): %s' %(xtal,xtal,str(nREFfree),xtal,str(nREF),str(CC)))
                    if errorMessage != '':
                        self.Logfile.insert('pointless failed with the following error: %s' %errorMessage)

                    try:
                        if float(CC) < 0.95:
                            self.Logfile.insert('correlation coefficient between the two files is below 0.95; will try to understand from dimple.log which one was used for initial map calculation')
                            if os.path.isfile('dimple/dimple_rerun_on_selected_file/dimple/dimple.log'):
                                foundLine=False
                                mtzin=''
                                for line in open('dimple/dimple_rerun_on_selected_file/dimple/dimple.log'):
                                    if foundLine:
                                        mtzin=line.replace(' ','').replace('\n','').replace('\r','')
                                        self.Logfile.insert('%s was used for inital map calculation' %mtzin)
                                        break
                                    if line.startswith(' --no-cleanup'):
                                        foundLine=True

                                if os.path.isfile(mtzin):
                                    self.Logfile.insert('%s: mtzfile used for refinement is not the same as the one chosen from autoprocessing' %xtal)
                                    self.Logfile.insert('%s: current mtzfile after autoprocessing: %s' %(xtal,os.path.realpath(xtal+'.mtz')))
                                    self.Logfile.insert('%s: removing links for %s.mtz/%s.log' %(xtal,xtal,xtal))
                                    os.system('/bin/rm %s.mtz 2> /dev/null' %xtal)
                                    os.system('/bin/rm %s.log 2> /dev/null' %xtal)
                                    self.Logfile.insert('linking %s to %s.mtz' %(os.path.relpath(mtzin),xtal))
                                    os.symlink(os.path.relpath(mtzin),xtal+'.mtz')
                                    print 'mtzin',mtzin
                                    print "mtzin[mtzin.rfind('/'):]",mtzin[mtzin.rfind('/'):]
                                    print "os.path.join(mtzin[mtzin.rfind('/'):]),'*log')",os.path.join(mtzin[mtzin.rfind('/'):],'*log')
                                    for logfile in glob.glob(os.path.join(mtzin[mtzin.rfind('/'):],'*log')):
                                        self.Logfile.insert('linking %s to %s.log' %(logfile,xtal))
                                        os.symlink(os.path.relpath(logfile),xtal+'.log')
                                        break




#                            for mtzfile in glob.glob('autoprocessing/*/%s.mtz' %xtal):
#                                self.Logfile.insert('checking %s' %mtzfile)
#                                procMTZ=mtztools(mtzfile)
#                                nREF=procMTZ.get_number_measured_reflections()
#                                CC=freeMTZ.calculate_correlaton_between_mtzfiles(mtzfile)
#                                self.Logfile.insert('%s: calculating CC between refine.mtz (%s refl) and %s (%s refl): %s' %(xtal,str(nREFfree),mtzfile.split('/')[1],str(nREF),str(CC)))
#                            self.Logfile.insert('correlation coefficient between the two files is below 0.9; will search autoprocessing directory for file with higher CC')
#                            for mtzfile in glob.glob('autoprocessing/*/%s.mtz' %xtal):
#                                self.Logfile.insert('checking %s' %mtzfile)
#                                procMTZ=mtztools(mtzfile)
#                                nREF=procMTZ.get_number_measured_reflections()
#                                CC=freeMTZ.calculate_correlaton_between_mtzfiles(mtzfile)
#                                self.Logfile.insert('%s: calculating CC between refine.mtz (%s refl) and %s (%s refl): %s' %(xtal,str(nREFfree),mtzfile.split('/')[1],str(nREF),str(CC)))


                    except ValueError:
                        self.Logfile.insert('something went wrong: calculated CC value does not seem to be a floating point number')


#            for mtzfile in glob.glob('autoprocessing/*/%s.mtz' %xtal):
#                procMTZ=mtztools(mtzfile)
#                nREF=procMTZ.get_number_measured_reflections()
#                resoHIGH=procMTZ.get_high_resolution_from_mtz()
#                CC=refineMTZ.calculate_correlaton_between_mtzfiles(mtzfile)
#                self.Logfile.insert('comparing refine.mtz ')
#                if nREF==nREFrefine:
#                    if os.path.isfile(xtal+'.mtz'):
#                        if os.path.realpath(xtal+'.mtz') != os.path.realpath(mtzfile):
#                            self.Logfile.insert('%s: mtzfile used for refinement is not the same as the one chosen from autoprocessing' %xtal)
#                            self.Logfile.insert('%s: current mtzfile after autoprocessing: %s' %(xtal,os.path.realpath(xtal+'.mtz')))
#                            self.Logfile.insert('%s: removing links for %s.mtz/%s.log' %(xtal,xtal,xtal))
#                            os.system('/bin/rm %s.mtz 2> /dev/null' %xtal)
#                            os.system('/bin/rm %s.log 2> /dev/null' %xtal)
#                            self.Logfile.insert('linking %s to %s.mtz' %(mtzfile,xtal))
#                            os.symlink(mtzfile,xtal+'.mtz')
#                            self.Logfile.insert('linking %s to %s.log' %(mtzfile.replace('.mtz','.log'),xtal))
#                            os.symlink(mtzfile.replace('.mtz','.log'),xtal+'.log')
#                            break

        found_logfile=False
        if os.path.isfile(xtal+'.log'):
            found_logfile=True
#            db_dict['DataProcessingPathToLogfile']=os.path.realpath(xtal+'.log').replace(os.getcwd()+'/','')
            db_dict['DataProcessingPathToLogfile']=os.path.realpath(xtal+'.log')
            db_dict['DataProcessingLOGfileName']=xtal+'.log'
            if db_dict['DataCollectionOutcome']=='None' or db_dict['DataCollectionOutcome']=='':
                db_dict['DataCollectionOutcome']='success'
            aimless_results=parse().read_aimless_logfile(xtal+'.log')
            db_dict.update(aimless_results)
        else:
            db_dict['DataProcessingPathToLogfile']=''
            db_dict['DataProcessingLOGfileName']=''

        # MTZ file

#        found_mtzfile=self.find_file(xtal+'.mtz',xtal)
#        if found_mtzfile:
        if os.path.isfile(xtal+'.mtz'):
            db_dict['DataProcessingPathToMTZfile']=os.path.realpath(xtal+'.mtz')
#            db_dict['DataProcessingPathToMTZfile']=os.path.realpath(xtal+'.mtz').replace(os.getcwd()+'/','')
            db_dict['DataProcessingMTZfileName']=xtal+'.mtz'
            if not found_logfile:
                mtz_info=mtztools(xtal+'.mtz').get_information_for_datasource()
                db_dict.update(mtz_info)
                db_dict['DataCollectionOutcome']='success'
        else:
            db_dict['DataProcessingPathToMTZfile']=''
            db_dict['DataProcessingMTZfileName']=''

        return db_dict


    def sync_dimple_results(self,xtal,db_dict):

        # DIMPLE pdb

#        found_dimple_pdb=self.find_file('dimple.pdb',xtal)
#        if found_dimple_pdb:
        if os.path.isfile('dimple.pdb'):
#            db_dict['DimplePathToPDB']=os.path.realpath('dimple.pdb').replace(os.getcwd()+'/','')
            db_dict['DimplePathToPDB']=os.path.realpath('dimple.pdb')
            pdb_info=parse().PDBheader('dimple.pdb')
            db_dict['DimpleRcryst']=pdb_info['Rcryst']
            db_dict['DimpleRfree']=pdb_info['Rfree']
            db_dict['DimpleResolutionHigh']=pdb_info['ResolutionHigh']
            db_dict['DimpleStatus']='finished'
        else:
            db_dict['DimplePathToPDB']=''
            db_dict['DimpleRcryst']=''
            db_dict['DimpleRfree']=''
            db_dict['DimpleResolutionHigh']=''
            db_dict['DimpleStatus']='pending'
        # DIMPLE mtz

        dimple_path=''
#        found_dimple_mtz=self.find_file('dimple.mtz',xtal)
#        if found_dimple_mtz:
        if os.path.isfile('dimple.mtz'):
#            db_dict['DimplePathToMTZ']=os.path.realpath('dimple.mtz').replace(os.getcwd()+'/','')
            db_dict['DimplePathToMTZ']=os.path.realpath('dimple.mtz')
            dimple_mtz=db_dict['DimplePathToMTZ']
            dimple_path=dimple_mtz[:dimple_mtz.rfind('/')]
        else:
            db_dict['DimplePathToMTZ']=''
            db_dict['DimpleStatus']='pending'

        if os.path.isfile(os.path.join(dimple_path,'dimple','dimple_rerun_on_selected_file','dimple_run_in_progress')):
            db_dict['DimpleStatus']='running'

        # MTZ free file

#        found_free_mtz=self.find_file(xtal+'.free.mtz',xtal)
#        if found_free_mtz:
        if os.path.isfile(xtal+'.free.mtz'):
            db_dict['RefinementMTZfree']=os.path.realpath(xtal+'.free.mtz')
#            db_dict['RefinementMTZfree']=os.path.realpath(xtal+'.free.mtz').replace(os.getcwd()+'/','')
        else:
            db_dict['RefinementMTZfree']=''
            if os.path.isfile(os.path.join(dimple_path,'prepared2.mtz')):
                os.symlink(os.path.relpath(os.path.join(dimple_path,'prepared2.mtz')),xtal+'.free.mtz')
#                db_dict['RefinementMTZfree']=os.path.realpath(xtal+'.free.mtz').replace(os.getcwd()+'/','')
                db_dict['RefinementMTZfree']=os.path.realpath(xtal+'.free.mtz')
#                db_dict['RefinementMTZfree']=xtal+'.free.mtz'
            elif os.path.isfile(os.path.join(dimple_path,'prepared.mtz')):
                os.symlink(os.path.relpath(os.path.join(dimple_path,'prepared.mtz')),xtal+'.free.mtz')
                db_dict['RefinementMTZfree']=os.path.realpath(xtal+'.free.mtz')
#                db_dict['RefinementMTZfree']=os.path.realpath(xtal+'.free.mtz').replace(os.getcwd()+'/','')
#               db_dict['RefinementMTZfree']=xtal+'.free.mtz'

        return db_dict


    def sync_compound_information(self,db_dict):
        # only update database if SMILES or compoundID field is blank!

        compoundID=db_dict['CompoundCode']
        if compoundID=='None' or compoundID=='':
            if os.path.isdir('compound'):
                for smiles in glob.glob('compound/*'):
                    if smiles.endswith('smiles'):
                        for line in open(smiles):
                            if len(line.split()) >= 1:
                                db_dict['CompoundCode']=smiles[smiles.rfind('/')+1:smiles.rfind('.')]
                                compoundID=db_dict['CompoundCode']
                                break

        if os.path.isfile(compoundID+'.cif') and os.path.getsize(compoundID+'.cif') > 20:
            db_dict['RefinementCIF']=os.path.realpath(compoundID+'.cif').replace(os.getcwd()+'/','')
            db_dict['RefinementCIFStatus']='restraints generated'
        else:
            os.system('/bin/rm %s.cif 2> /dev/null' %compoundID)
            os.system('/bin/rm compound/%s.cif 2> /dev/null' %compoundID)
            db_dict['RefinementCIF']=''
            db_dict['RefinementCIFStatus']='pending'

        smilesDB=db_dict['CompoundSMILES']
        smiles_found=True
        if smilesDB=='None' or smilesDB=='':
            smiles_found=False
            if os.path.isdir('compound'):
                for smiles in glob.glob('compound/*'):
                    if smiles.endswith('smiles'):
                        for line in open(smiles):
                            if len(line.split()) >= 1:
                                db_dict['CompoundSMILES']=line.split()[0]
                                smilesDB=db_dict['CompoundSMILES']
                                smiles_found=True
                                break

        if not smiles_found:
            db_dict['RefinementCIFStatus']='missing smiles'

        if not os.path.isfile(compoundID+'.pdb') or  os.path.getsize(compoundID+'.pdb') < 20:
            os.system('/bin/rm %s.pdb 2> /dev/null' %compoundID)
            os.system('/bin/rm compound/%s.pdb 2> /dev/null' %compoundID)

        if not os.path.isfile(compoundID+'.png') or  os.path.getsize(compoundID+'.png') < 20:
            os.system('/bin/rm %s.png 2> /dev/null' %compoundID)
            os.system('/bin/rm compound/%s.png 2> /dev/null' %compoundID)

        return db_dict


    def sync_refinement_results(self,xtal,db_dict):

        #
        # REFINE pdb
        #

#        found_refine_pdb=self.find_file('refine.pdb',xtal)
#        if found_refine_pdb and not 'dimple' in os.path.realpath('refine.pdb'):
        if os.path.isfile('refine.pdb') and not 'dimple' in os.path.realpath('refine.pdb'):
            db_dict['RefinementPDB_latest']=os.path.realpath('refine.pdb')
#            db_dict['RefinementPDB_latest']=os.path.realpath('refine.pdb').replace(os.getcwd()+'/','')
            db_dict['RefinementStatus']='finished'
            pdb_info=parse().dict_for_datasource_update('refine.pdb')
            db_dict.update(pdb_info)
            if db_dict['RefinementOutcome']=='None' or db_dict['RefinementOutcome']=='':
                db_dict['RefinementOutcome']='3 - In Refinement'
            elif str(db_dict['RefinementOutcome']).startswith('1'):
                db_dict['RefinementOutcome']='3 - In Refinement'
            elif str(db_dict['RefinementOutcome']).startswith('2'):
                db_dict['RefinementOutcome']='3 - In Refinement'
        else:
            db_dict['RefinementPDB_latest']=''
            db_dict['RefinementStatus']='pending'
            db_dict['RefinementOutcome']='1 - Analysis Pending'
            os.system('/bin/rm refine.pdb 2> /dev/null')

        if os.path.isfile('REFINEMENT_IN_PROGRESS'):
            db_dict['RefinementStatus']='running'

        #
        # REFINE bound pdb
        #

        if os.path.isfile('refine.bound.pdb'):
            db_dict['RefinementBoundConformation']=os.path.realpath('refine.bound.pdb')
#            db_dict['RefinementBoundConformation']='refine.bound.pdb'
        else:
#            if os.path.isfile('refine.pdb'):
#                os.system("giant.strip_conformations pdb=refine.pdb suffix='.bound.pdb'")
#                if os.path.isfile('refine.bound.pdb'):
#                    db_dict['RefinementBoundConformation']=os.path.realpath('refine.bound.pdb')
#                else:
#                    db_dict['RefinementBoundConformation']=''
#            else:
#                db_dict['RefinementBoundConformation']=''
            db_dict['RefinementBoundConformation']=''

        #
        # REFINE mtz
        #

#        found_refine_mtz=self.find_file('refine.mtz',xtal)
#        if found_refine_mtz and not 'dimple' in os.path.realpath('refine.mtz'):
        if os.path.isfile('refine.mtz') and not 'dimple' in os.path.realpath('refine.mtz'):
            db_dict['RefinementMTZ_latest']=os.path.realpath('refine.mtz')
#            db_dict['RefinementMTZ_latest']=os.path.realpath('refine.mtz').replace(os.getcwd()+'/','')
        else:
            db_dict['RefinementMTZ_latest']=''
            os.system('/bin/rm refine.mtz 2> /dev/null')

        return db_dict


    def find_apo_structures_for_PanDDA(self,panddaPATH):

        # first check if structure is already present in DB and if so if all the
        # information concur

        # need to update pandda directory for every exported structure so that
        # we know where to look for the pandda.log file that contains the relevant information

        # update CrystalName_of_pandda_input in DB

        # in DB: update StructureType field accordingly

        # newer pandda versions seem to have severl copies of pandda.log with names like
        # pandda-2016-09-01-2139.log
        panddaLog=glob.glob(os.path.join(panddaPATH,'pandda*log'))
        panddaLog.sort(key=os.path.getmtime)

        panddaVersion='unknown'
        readindApoStructures = False
        apoStructures = []
        apoStructureDict = {}
        apoString=''
        for files in panddaLog:
            for line in open(files):
                if line.startswith('-  Pandda Version'):
                    if len(line.split()) >= 4:
                        panddaVersion=line.split()[3]
                if 'No Statistical Maps Found:' in line:
                    readindApoStructures=True
                if readindApoStructures:
                    if 'Pickling Object: processed_datasets' in line:
                        if line.split() >= 2:
                            # e.g. line.split() -> ['Pickling', 'Object:', 'processed_datasets/NUDT22A-x0055/pickles/dataset.pickle']
                            xtal=line.split()[2].split('/')[1]
                            if os.path.isfile(os.path.join(panddaPATH,'processed_datasets',xtal,xtal+'-pandda-input.pdb')):
                                apoStructures.append(xtal)
                                apoString+=xtal+';'
                if 'Pre-existing statistical maps (from previous runs) have been found and will be reused:' in line:
                    readindApoStructures=False
            apoStructureDict[panddaPATH]=apoStructures

        return apoString[:-1]

#    def sync_pandda_table(self):
#        progress_step=1
#        progress=0
#        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
#
#        # also need to update PANDDA table...
#        pandda_models=self.db.execute_statement("select CrystalName,PANDDA_site_index,PANDDA_site_event_index,PANDDA_site_x,PANDDA_site_y,PANDDA_site_z,PANDDApath,ApoStructures from panddaTable")
#        if len(pandda_models) > 0:
#            progress_step=100/float(len(pandda_models))
#
#        if pandda_models != []:
#            for entry in pandda_models:
#                db_pandda_dict={}
#                xtal=entry[0]
#                site_index=entry[1]
#                event_index=entry[2]
#                panddaPATH=entry[6]
#                apoStructures=entry[7]
#                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'checking %s -> site %s -> event %s ' %(xtal,site_index,event_index))
#                try:
#                    event_x = float(str(entry[3]))
#                    event_y = float(str(entry[4]))
#                    event_z = float(str(entry[5]))
#                except ValueError:
#                    pass
#
#                # do not update pandda path since this one is updated during pandda export!
#                # instead try to get apo semi-colon separated list of apo structures that were used to
#                # calculate event maps; but only if field is blank!
##                db_pandda_dict['PANDDApath']=self.panddas_directory
#                if str(apoStructures)=='None' or apoStructures=='':
#                    if panddaPATH != 'None' or panddaPATH != '':
#                        self.Logfile.insert('trying to find which apo structures were used to calculate the event maps in '+panddaPATH)
#                        db_pandda_dict['ApoStructures']=self.find_apo_structures_for_PanDDA(panddaPATH)
#                    else:
#                        self.Logfile.insert('pandda path for '+xtal+' is empty in database')
#
#                # event map
#
#                found_event_map=False
#                db_pandda_dict['PANDDA_site_event_map']=''
#                for file in glob.glob(os.path.join(self.initial_model_directory,xtal,'*ccp4')):
#                    filename=file[file.rfind('/')+1:]
#                    if filename.startswith(xtal+'-event_'+event_index) and filename.endswith('map.native.ccp4'):
#                        event_map=file
#                        db_pandda_dict['PANDDA_site_event_map']=os.path.realpath(event_map)
##                        db_pandda_dict['PANDDA_site_event_map']=os.path.realpath(event_map).replace(os.getcwd()+'/','')
#                        found_event_map=True
#                        break
#                if not found_event_map:
#                    db_pandda_dict['PANDDA_site_event_map']=''
#
#                db_pandda_dict['PANDDA_site_initial_model']=''
#                for file in glob.glob(os.path.join(self.initial_model_directory,xtal,'*pdb')):
#                    filename=file[file.rfind('/')+1:]
#                    if filename.endswith('-ensemble-model.pdb'):
#                        db_pandda_dict['PANDDA_site_initial_model']=os.path.realpath(file).replace(os.getcwd()+'/','')
#                        break
#
#                db_pandda_dict['PANDDA_site_initial_mtz']=''
#                for file in glob.glob(os.path.join(self.initial_model_directory,xtal,'*mtz')):
#                    filename=file[file.rfind('/')+1:]
#                    if filename.endswith('pandda-input.mtz'):
#                        db_pandda_dict['PANDDA_site_initial_mtz']=os.path.realpath(file).replace(os.getcwd()+'/','')
#                        break
#
#
#                db_pandda_dict['PANDDA_site_ligand_resname'] = ''
#                db_pandda_dict['PANDDA_site_ligand_chain'] = ''
#                db_pandda_dict['PANDDA_site_ligand_sequence_number'] = ''
#                db_pandda_dict['PANDDA_site_ligand_altLoc'] = ''
#                db_pandda_dict['PANDDA_site_ligand_placed'] = 'False'
#                db_pandda_dict['PANDDA_site_spider_plot'] = ''
#                db_pandda_dict['PANDDA_site_ligand_id'] = ''
#
#                db_pandda_dict['PANDDA_site_occupancy'] = ''
#                db_pandda_dict['PANDDA_site_B_average'] = ''
#                db_pandda_dict['PANDDA_site_B_ratio_residue_surroundings'] = ''
#                db_pandda_dict['PANDDA_site_rmsd'] = ''
#                db_pandda_dict['PANDDA_site_RSCC'] = ''
#                db_pandda_dict['PANDDA_site_RSR'] = ''
#                db_pandda_dict['PANDDA_site_RSZD'] = ''
#
#                if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'refine.pdb')):
#                    ligands_in_file=pdbtools(os.path.join(self.initial_model_directory,xtal,'refine.pdb')).find_xce_ligand_details()
#
#                    for ligand in ligands_in_file:
#                        residue_name=   ligand[0]
#                        residue_chain=  ligand[1]
#                        residue_number= ligand[2]
#                        residue_altLoc= ligand[3]
#                        residue_xyz = pdbtools(os.path.join(self.initial_model_directory,xtal,'refine.pdb')).get_center_of_gravity_of_residue_ish(residue_chain, residue_number)
#                        distance = misc().calculate_distance_between_coordinates(residue_xyz[0], residue_xyz[1],residue_xyz[2],event_x, event_y,event_z)
#                        # if coordinate of ligand and event are closer than 7A, then we assume they belong together
#                        if distance < 7:
#                            db_pandda_dict['PANDDA_site_ligand_resname'] = residue_name
#                            db_pandda_dict['PANDDA_site_ligand_chain'] = residue_chain
#                            db_pandda_dict['PANDDA_site_ligand_sequence_number'] = residue_number
#                            db_pandda_dict['PANDDA_site_ligand_altLoc'] = residue_altLoc
#                            db_pandda_dict['PANDDA_site_ligand_placed'] = 'True'
#                            db_pandda_dict['PANDDA_site_ligand_id']=residue_name+'-'+residue_chain+'-'+residue_number
#
#                            if xtal+'/Refine_' in os.path.realpath(os.path.join(self.initial_model_directory,xtal,'refine.pdb')):
#                                tmp=os.path.realpath(os.path.join(self.initial_model_directory,xtal,'refine.pdb'))
#                                spider_plot=os.path.join(tmp[:tmp.rfind('/')],'residue_plots',residue_name+'-'+residue_chain+'-'+residue_number+'.png')
#                                if os.path.isfile(spider_plot):
#                                    db_pandda_dict['PANDDA_site_spider_plot']=os.path.realpath(spider_plot)
##                                    db_pandda_dict['PANDDA_site_spider_plot']=os.path.realpath(spider_plot).replace(os.getcwd()+'/','')
#                                if os.path.isfile(os.path.join(tmp[:tmp.rfind('/')],'residue_scores.csv')):
#
#                                    with open(os.path.join(tmp[:tmp.rfind('/')],'residue_scores.csv'), 'rb') as csv_import:
#                                        csv_dict = csv.DictReader(csv_import)
#                                        for i, line in enumerate(csv_dict):
#                                            residueNameChainNumber = line['']
#                                            if residueNameChainNumber == residue_name+'-'+residue_chain+'-'+residue_number:
#                                                db_pandda_dict['PANDDA_site_occupancy'] = line['Occupancy']
#                                                db_pandda_dict['PANDDA_site_B_average'] = line['Average B-factor (Residue)']
#                                                db_pandda_dict['PANDDA_site_B_ratio_residue_surroundings'] = line['Surroundings B-factor Ratio']
#                                                db_pandda_dict['PANDDA_site_rmsd'] = line['Model RMSD']
#                                                db_pandda_dict['PANDDA_site_RSCC'] = line['RSCC']
#                                                db_pandda_dict['PANDDA_site_RSR'] = line['RSR']
#                                                db_pandda_dict['PANDDA_site_RSZD'] = line['RSZD']
#                            break
#
#                if db_pandda_dict != {}:
##                    self.db.update_panddaTable(xtal, site_index, db_pandda_dict)
#                    self.db.update_site_event_panddaTable(xtal, site_index, event_index, db_pandda_dict)
#                    self.Logfile.insert('updating panddaTable for xtal: %s, site: %s' %(xtal,site_index))
#
#                progress += progress_step
#                self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
#


    def sync_pandda_table_NEW(self):
        progress_step=1
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        # also need to update PANDDA table...
        pandda_models=self.db.execute_statement("select CrystalName,PANDDA_site_index,PANDDA_site_event_index,PANDDA_site_x,PANDDA_site_y,PANDDA_site_z,PANDDApath,ApoStructures from panddaTable")
        if len(pandda_models) > 0:
            progress_step=100/float(len(pandda_models))

        if pandda_models != []:
            for entry in pandda_models:
                db_pandda_dict={}
                xtal=entry[0]
                site_index=entry[1]
                event_index=entry[2]
                panddaPATH=entry[6]
                apoStructures=entry[7]
                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'checking %s -> site %s -> event %s ' %(xtal,site_index,event_index))
                self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
                progress += progress_step

                try:
                    event_x = float(str(entry[3]))
                    event_y = float(str(entry[4]))
                    event_z = float(str(entry[5]))
                except ValueError:
                    pass

                # do not update pandda path since this one is updated during pandda export!
                # instead try to get apo semi-colon separated list of apo structures that were used to
                # calculate event maps; but only if field is blank!
#                db_pandda_dict['PANDDApath']=self.panddas_directory
                if str(apoStructures)=='None' or apoStructures=='':
                    if panddaPATH != 'None' or panddaPATH != '':
                        self.Logfile.insert('trying to find which apo structures were used to calculate the event maps in '+panddaPATH)
                        db_pandda_dict['ApoStructures']=self.find_apo_structures_for_PanDDA(panddaPATH)
                    else:
                        self.Logfile.insert('pandda path for '+xtal+' is empty in database')

                # event map

                found_event_map=False
                db_pandda_dict['PANDDA_site_event_map']=''
                for file in glob.glob(os.path.join(self.initial_model_directory,xtal,'*ccp4')):
                    filename=file[file.rfind('/')+1:]
                    if filename.startswith(xtal+'-event_'+event_index) and filename.endswith('map.native.ccp4'):
                        event_map=file
                        db_pandda_dict['PANDDA_site_event_map']=os.path.realpath(event_map)
#                        db_pandda_dict['PANDDA_site_event_map']=os.path.realpath(event_map).replace(os.getcwd()+'/','')
                        found_event_map=True
                        break
                if not found_event_map:
                    db_pandda_dict['PANDDA_site_event_map']=''

                db_pandda_dict['PANDDA_site_initial_model']=''
                for file in glob.glob(os.path.join(self.initial_model_directory,xtal,'*pdb')):
                    filename=file[file.rfind('/')+1:]
                    if filename.endswith('-ensemble-model.pdb'):
                        db_pandda_dict['PANDDA_site_initial_model']=os.path.realpath(file).replace(os.getcwd()+'/','')
                        break

                db_pandda_dict['PANDDA_site_initial_mtz']=''
                for file in glob.glob(os.path.join(self.initial_model_directory,xtal,'*mtz')):
                    filename=file[file.rfind('/')+1:]
                    if filename.endswith('pandda-input.mtz'):
                        db_pandda_dict['PANDDA_site_initial_mtz']=os.path.realpath(file).replace(os.getcwd()+'/','')
                        break


                db_pandda_dict['PANDDA_site_ligand_resname'] = ''
                db_pandda_dict['PANDDA_site_ligand_chain'] = ''
                db_pandda_dict['PANDDA_site_ligand_sequence_number'] = ''
                db_pandda_dict['PANDDA_site_ligand_altLoc'] = ''
                db_pandda_dict['PANDDA_site_ligand_placed'] = 'False'
                db_pandda_dict['PANDDA_site_spider_plot'] = ''
                db_pandda_dict['PANDDA_site_ligand_id'] = ''

                db_pandda_dict['PANDDA_site_occupancy'] = ''
                db_pandda_dict['PANDDA_site_B_average'] = ''
                db_pandda_dict['PANDDA_site_B_ratio_residue_surroundings'] = ''
                db_pandda_dict['PANDDA_site_rmsd'] = ''
                db_pandda_dict['PANDDA_site_RSCC'] = ''
                db_pandda_dict['PANDDA_site_RSR'] = ''
                db_pandda_dict['PANDDA_site_RSZD'] = ''

                if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'refine.pdb')):
                    ligands_in_file=pdbtools(os.path.join(self.initial_model_directory,xtal,'refine.pdb')).find_xce_ligand_details()
                    if ligands_in_file == []:
                        self.Logfile.warning('%s: could not find any ligands in refine.pdb' %xtal)
                        continue
                    else:
                        self.Logfile.insert('%s: found the following ligands in refine.pdb: %s' %(xtal,str(ligands_in_file)))

                    distanceList=[]
                    for ligand in ligands_in_file:
                        residue_name=   ligand[0]
                        residue_chain=  ligand[1]
                        residue_number= ligand[2]
                        residue_altLoc= ligand[3]
                        residue_xyz = pdbtools(os.path.join(self.initial_model_directory,xtal,'refine.pdb')).get_center_of_gravity_of_residue_ish(residue_chain, residue_number)
                        distance = misc().calculate_distance_between_coordinates(residue_xyz[0], residue_xyz[1],residue_xyz[2],event_x, event_y,event_z)
                        distanceList.append([distance,residue_name,residue_chain,residue_number,residue_altLoc])
                        self.Logfile.insert('%s: calculating distance between event and ligand (%s %s %s): %s' %(xtal,residue_name,residue_chain,residue_number,str(distance)))

                    # now take the ligand that is closest to the event
                    try:
                        smallestDistance = min(distanceList, key=lambda x: x[0])
                    except ValueError:
                        self.Logfile.error('could not determine smallest distance between current ligand pandda events')
                        continue
                    distance =          smallestDistance[0]
                    residue_name =      smallestDistance[1]
                    residue_chain =     smallestDistance[2]
                    residue_number =    smallestDistance[3]
                    residue_altLoc =    smallestDistance[4]
                    self.Logfile.insert('%s: ligand with the shorted distance (%sA) to the current event (id: %s): %s %s %s %s' %(xtal,str(distance),event_index,residue_name,residue_chain,residue_number,residue_altLoc))
                    db_pandda_dict['PANDDA_site_ligand_resname'] = residue_name
                    db_pandda_dict['PANDDA_site_ligand_chain'] = residue_chain
                    db_pandda_dict['PANDDA_site_ligand_sequence_number'] = residue_number
                    db_pandda_dict['PANDDA_site_ligand_altLoc'] = residue_altLoc
                    db_pandda_dict['PANDDA_site_ligand_placed'] = 'True'
                    db_pandda_dict['PANDDA_site_ligand_id']=residue_name+'-'+residue_chain+'-'+residue_number

                    if xtal+'/Refine_' in os.path.realpath(os.path.join(self.initial_model_directory,xtal,'refine.pdb')):
                        tmp=os.path.realpath(os.path.join(self.initial_model_directory,xtal,'refine.pdb'))
                        spider_plot=os.path.join(tmp[:tmp.rfind('/')],'residue_plots',residue_name+'-'+residue_chain+'-'+residue_number+'.png')
                        if os.path.isfile(spider_plot):
                            db_pandda_dict['PANDDA_site_spider_plot']=os.path.realpath(spider_plot)
                        if os.path.isfile(os.path.join(tmp[:tmp.rfind('/')],'residue_scores.csv')):
                            with open(os.path.join(tmp[:tmp.rfind('/')],'residue_scores.csv'), 'rb') as csv_import:
                                csv_dict = csv.DictReader(csv_import)
                                for i, line in enumerate(csv_dict):
                                    residueNameChainNumber = line['']
                                    if residueNameChainNumber == residue_name+'-'+residue_chain+'-'+residue_number:
                                        db_pandda_dict['PANDDA_site_occupancy'] = line['Occupancy']
                                        db_pandda_dict['PANDDA_site_B_average'] = line['Average B-factor (Residue)']
                                        db_pandda_dict['PANDDA_site_B_ratio_residue_surroundings'] = line['Surroundings B-factor Ratio']
                                        db_pandda_dict['PANDDA_site_rmsd'] = line['Model RMSD']
                                        db_pandda_dict['PANDDA_site_RSCC'] = line['RSCC']
                                        db_pandda_dict['PANDDA_site_RSR'] = line['RSR']
                                        db_pandda_dict['PANDDA_site_RSZD'] = line['RSZD']
#                    break

                if db_pandda_dict != {}:
#            self.db.update_panddaTable(xtal, site_index, db_pandda_dict)
                    self.db.update_site_event_panddaTable(xtal, site_index, event_index, db_pandda_dict)
                    self.Logfile.insert('updating panddaTable for xtal: %s, site: %s' %(xtal,site_index))

#                progress += progress_step
#                self.emit(QtCore.SIGNAL('update_progress_bar'), progress)





class create_png_and_cif_of_compound(QtCore.QThread):
    def __init__(self,external_software,initial_model_directory,compound_list,database_directory,data_source_file,todo,ccp4_scratch_directory,xce_logfile,max_queue_jobs,restraints_program):
        QtCore.QThread.__init__(self)
        self.external_software=external_software
        self.initial_model_directory=initial_model_directory
        self.compound_list=compound_list
        self.database_directory=database_directory
        self.data_source_file=data_source_file
        self.todo=todo
        self.ccp4_scratch_directory=ccp4_scratch_directory
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.max_queue_jobs=max_queue_jobs
        self.db=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
        self.restraints_program=restraints_program

    def run(self):
        # first remove all ACEDRG input scripts in ccp4_scratch directory
        self.Logfile.insert('removing all xce_acedrg scripts from '+self.ccp4_scratch_directory)
        os.chdir(self.ccp4_scratch_directory)
        os.system('/bin/rm -f xce_acedrg*')

        progress_step=100/float(len(self.compound_list))
        progress=0
        counter=1
        for item in self.compound_list:
            sampleID=item[0]
            compoundID=item[1]
            smiles=item[2]
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'creating ACEDRG shell script for '+sampleID)
            if compoundID=='' or compoundID==None:
                compoundID='compound'

            if not os.path.isdir(os.path.join(self.initial_model_directory,sampleID)):
                os.mkdir(os.path.join(self.initial_model_directory,sampleID))

            if self.todo=='ALL' or self.todo=='SELECTED':
                # remove symbolic links if present
                if os.path.isfile(os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.pdb')):
                    os.system('/bin/rm '+os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.pdb'))
#               commented this out since people found the presence of the old.cif file interfering with pandda.inspect
#                if os.path.isfile(os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.cif')):
#                    # copy existing CIF file to old.cif so that it can be used as input in
#                    # restraints generating program
#                    os.chdir(os.path.join(self.initial_model_directory,sampleID))
#                    os.system('/bin/cp %s old.cif' %(compoundID.replace(' ','')+'.cif'))
#                    os.system('/bin/rm '+os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.cif'))
                if os.path.isfile(os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.png')):
                    os.system('/bin/rm '+os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.png'))
                if os.path.isdir(os.path.join(self.initial_model_directory,sampleID,'compound')):
                    os.system('/bin/rm -fr '+os.path.join(self.initial_model_directory,sampleID,'compound'))
                    db_dict={}
                    db_dict['RefinementCIFStatus']='pending'
                    self.Logfile.insert('%s: removed compound directory and all its contents' %sampleID)
                    self.Logfile.insert('%s: setting RefinementCIFStatus flag to started' %sampleID)
                    self.db.update_data_source(sampleID,db_dict)

            # create 'compound' directory if not present
            if not os.path.isdir(os.path.join(self.initial_model_directory,sampleID,'compound')):
                os.mkdir(os.path.join(self.initial_model_directory,sampleID,'compound'))

            # create text file which contains the smiles string
            if not os.path.isfile(os.path.join(self.initial_model_directory,sampleID,'compound',compoundID+'.smiles')):
                os.chdir(os.path.join(self.initial_model_directory,sampleID,'compound'))
                f=open(compoundID+'.smiles','w')
                f.write(smiles)
                f.close()

            if not os.path.isfile(os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.cif')):
                os.chdir(os.path.join(self.initial_model_directory,sampleID,'compound'))

                helpers().make_png( self.initial_model_directory,
                                    sampleID,
                                    compoundID,
                                    smiles,
                                    self.external_software,
                                    self.database_directory,
                                    self.data_source_file,
                                    self.ccp4_scratch_directory,
                                    counter,
                                    self.xce_logfile,
                                    self.restraints_program )
                counter += 1

                db_dict={}
                db_dict['RefinementCIFprogram']=self.restraints_program
                db_dict['RefinementCIFStatus']='started'
                self.Logfile.insert('%s: setting RefinementCIFStatus flag to started' %sampleID)
                self.db.update_data_source(sampleID,db_dict)

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        # submit array job at Diamond
        self.Logfile.insert('created input scripts for '+str(counter)+' ACEDRG jobs in '+self.ccp4_scratch_directory)
        os.chdir(self.ccp4_scratch_directory)
        self.Logfile.insert('changing directory to '+self.ccp4_scratch_directory)
        if counter > 1:
#            if os.getcwd().startswith('/dls'):
            if os.path.isdir('/dls'):
                if self.external_software['qsub_array']:
                    Cmds = (
                            '#PBS -joe -N xce_acedrg_master\n'
                            './xce_acedrg_$SGE_TASK_ID.sh\n'
                            )
                    f = open('acedrg_master.sh','w')
                    f.write(Cmds)
                    f.close()
                    self.Logfile.insert('submitting array job with maximal 100 jobs running on cluster')
                    self.Logfile.insert('using the following command:')
                    self.Logfile.insert('         qsub -P labxchem -t 1:%s -tc %s acedrg_master.sh' %(str(counter),self.max_queue_jobs))
                    os.system('qsub -P labxchem -t 1:%s -tc %s acedrg_master.sh' %(str(counter),self.max_queue_jobs))
                else:
                    self.Logfile.insert("cannot start ARRAY job: make sure that 'module load global/cluster' is in your .bashrc or .cshrc file")
            elif self.external_software['qsub']:
                self.Logfile.insert('submitting %s individual jobs to cluster' %(str(counter)))
                self.Logfile.insert('WARNING: this could potentially lead to a crash...')
                for i in range(counter):
                    self.Logfile.insert('qsub xce_acedrg_%s.sh' %(str(i+1)))
                    os.system('qsub xce_acedrg_%s.sh' %(str(i+1)))
            else:
                self.Logfile.insert('running %s consecutive ACEDRG jobs on your local machine')
                for i in range(counter):
                    self.Logfile.insert('starting xce_acedrg_%s.sh' %(str(i+1)))
                    os.system('./xce_acedrg_%s.sh' %(str(i+1)))

#        self.emit(QtCore.SIGNAL("finished()"))
        self.emit(QtCore.SIGNAL('datasource_menu_reload_samples'))

class run_dimple_on_all_autoprocessing_files(QtCore.QThread):
    def __init__(self,sample_list,initial_model_directory,external_software,ccp4_scratch_directory,database_directory,data_source_file,max_queue_jobs,xce_logfile):
        QtCore.QThread.__init__(self)
        self.sample_list=sample_list
        self.initial_model_directory=initial_model_directory
        self.external_software=external_software
        self.queueing_system_available=external_software['qsub']
        self.ccp4_scratch_directory=ccp4_scratch_directory
        self.database_directory=database_directory
        self.data_source_file=data_source_file
        self.max_queue_jobs=max_queue_jobs
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.pipeline='dimple'


    def run(self):
        progress_step=1
        if len(self.sample_list) != 0:
            progress_step=100/float(len(self.sample_list))
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        os.chdir(self.ccp4_scratch_directory)
        os.system('/bin/rm -f xce_dimple*sh')

        db=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
        database=os.path.join(self.database_directory,self.data_source_file)

        for n,item in enumerate(self.sample_list):

            xtal =                  item[0]
            visit_run_autoproc =    item[1]
            mtzin =                 item[2]
            ref_pdb =               item[3]
            ref_mtz =               item[4]
            ref_cif =               item[5]

            # check if reference mtzfile has an Rfree column; if not, then ignore
            # DIMPLE assumes an Rfree column and barfs if it is not present
            # note: ref_mtz looks like this: ref mtz  -R reference.mtz
#            if os.path.isfile(ref_mtz.split()[len(ref_mtz.split())-1]):
#                mtz_column_dict=mtztools(ref_mtz.split()[len(ref_mtz.split())-1]).get_all_columns_as_dict()
            if os.path.isfile(ref_mtz):
                mtz_column_dict=mtztools(ref_mtz).get_all_columns_as_dict()
                if 'FreeR_flag' not in mtz_column_dict['RFREE']:
                    self.Logfile.insert('cannot find FreeR_flag in reference mtz file: %s -> ignoring reference mtzfile!!!' %ref_mtz)
                    ref_mtz = ''
                    if mtz_column_dict['RFREE'] != []:
                        self.Logfile.insert('found Rfree set with other column name though: %s' %str(mtz_column_dict['RFREE']))
                        self.Logfile.insert('try renaming Rfree column to FreeR_flag with CAD!')


            db_dict={}
            db_dict['DimpleReferencePDB']=ref_pdb
            db.update_data_source(xtal,db_dict)

            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'creating input script for '+xtal+' in '+visit_run_autoproc)

            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal)):
                os.mkdir(os.path.join(self.initial_model_directory,xtal))
            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'dimple')):
                os.mkdir(os.path.join(self.initial_model_directory,xtal,'dimple'))
            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc)):
                os.mkdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc))
            os.chdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc))
            os.system('touch dimple_run_in_progress')
            os.system('/bin/rm final.mtz 2> /dev/null')
            os.system('/bin/rm final.pdb 2> /dev/null')

            if self.queueing_system_available:
                top_line='#PBS -joe -N XCE_dimple\n'
            else:
                top_line='#!'+os.getenv('SHELL')+'\n'

            if 'csh' in os.getenv('SHELL'):
                ccp4_scratch='setenv CCP4_SCR '+self.ccp4_scratch_directory+'\n'
            elif 'bash' in os.getenv('SHELL'):
                ccp4_scratch='export CCP4_SCR='+self.ccp4_scratch_directory+'\n'
            else:
                ccp4_scratch=''

            if 'dimple_rerun_on_selected_file' in visit_run_autoproc:
                additional_cmds = (
                            'cd %s\n' %os.path.join(self.initial_model_directory,xtal) +
                            '/bin/rm dimple.pdb\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/final.pdb dimple.pdb\n'
                            '/bin/rm dimple.mtz\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/final.mtz dimple.mtz\n'
                            '/bin/rm 2fofc.map\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/2fofc.map .\n'
                            '/bin/rm fofc.map\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/fofc.map .\n'
                            '\n'
                            '$CCP4/libexec/python '+os.path.join(os.getenv('XChemExplorer_DIR'),'helpers','update_data_source_for_new_dimple_pdb.py')+
                            ' %s %s %s\n' %(os.path.join(self.database_directory,self.data_source_file),xtal,self.initial_model_directory)  )

            else:
                additional_cmds=''



            Cmds = (
                    '%s\n' %top_line+
                    '\n'
                    'export XChemExplorer_DIR="'+os.getenv('XChemExplorer_DIR')+'"\n'
                    '\n'
                    'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc) +
                    '\n'
                    'source $XChemExplorer_DIR/setup-scripts/xce.setup-sh\n'
                    '\n'
                    +ccp4_scratch+
                    '\n'
                    '$CCP4/bin/ccp4-python $XChemExplorer_DIR/helpers/update_status_flag.py %s %s %s %s\n' %(database,xtal,'DimpleStatus','running') +
                    '\n'
                    'dimple --no-cleanup %s %s %s %s dimple\n' %(mtzin,ref_pdb,ref_mtz,ref_cif) +
                    '\n'
                    'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc,'dimple') +
                    '\n'
                    'fft hklin final.mtz mapout 2fofc.map << EOF\n'
                    ' labin F1=FWT PHI=PHWT\n'
                    'EOF\n'
                    '\n'
                    'fft hklin final.mtz mapout fofc.map << EOF\n'
                    ' labin F1=DELFWT PHI=PHDELWT\n'
                    'EOF\n'
                    '\n'
                    +additional_cmds+
                    '\n'
                    'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc) +
                    '\n'
                    '/bin/rm dimple_run_in_progress\n'
                    '\n'
                    'ln -s dimple/final.pdb .\n'
                    'ln -s dimple/final.mtz .\n'
                    )

            os.chdir(self.ccp4_scratch_directory)
            f = open('xce_dimple_%s.sh' %str(n+1),'w')
            f.write(Cmds)
            f.close()
            os.system('chmod +x xce_dimple_%s.sh' %str(n+1))
            db_dict={}
            db_dict['DimpleStatus']='started'
            self.Logfile.insert('%s: setting DataProcessingStatus flag to started' %xtal)
            db.update_data_source(xtal,db_dict)


            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        # submit job
        self.Logfile.insert('created input scripts for '+str(n+1)+' in '+self.ccp4_scratch_directory)
        os.chdir(self.ccp4_scratch_directory)
#        if os.getcwd().startswith('/dls'):
        if os.path.isdir('/dls'):
            if self.external_software['qsub_array']:
                Cmds = (
                        '#PBS -joe -N xce_dimple_master\n'
                        './xce_dimple_$SGE_TASK_ID.sh\n'
                        )
                f = open('dimple_master.sh','w')
                f.write(Cmds)
                f.close()
                self.Logfile.insert('submitting array job with maximal 100 jobs running on cluster')
                self.Logfile.insert('using the following command:')
                self.Logfile.insert('qsub -P labxchem -t 1:%s -tc %s dimple_master.sh' %(str(n+1),self.max_queue_jobs))
                os.system('qsub -P labxchem -t 1:%s -tc %s dimple_master.sh' %(str(n+1),self.max_queue_jobs))
            else:
                self.Logfile.insert("cannot start ARRAY job: make sure that 'module load global/cluster' is in your .bashrc or .cshrc file")
        elif self.external_software['qsub']:
            self.Logfile.insert('submitting %s individual jobs to cluster' %(str(n+1)))
            self.Logfile.insert('WARNING: this could potentially lead to a crash...')
            for i in range(n+1):
                self.Logfile.insert('qsub xce_dimple_%s.sh' %(str(i+1)))
                os.system('qsub xce_dimple_%s.sh' %(str(i+1)))
        else:
            self.Logfile.insert('running %s consecutive DIMPLE jobs on your local machine')
            for i in range(n+1):
                self.Logfile.insert('starting xce_dimple_%s.sh' %(str(i+1)))
                os.system('./xce_dimple_%s.sh' %(str(i+1)))

        self.emit(QtCore.SIGNAL('datasource_menu_reload_samples'))

class run_dimple_on_all_autoprocessing_files_new(QtCore.QThread):
    def __init__(self,sample_list,initial_model_directory,external_software,ccp4_scratch_directory,database_directory,data_source_file,max_queue_jobs,xce_logfile):
        QtCore.QThread.__init__(self)
        self.sample_list=sample_list
        self.initial_model_directory=initial_model_directory
        self.external_software=external_software
        self.queueing_system_available=external_software['qsub']
        self.ccp4_scratch_directory=ccp4_scratch_directory
        self.database_directory=database_directory
        self.data_source_file=data_source_file

        self.db=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
        self.database=os.path.join(self.database_directory,self.data_source_file)

        self.max_queue_jobs=max_queue_jobs
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)

        self.n=0

        self.pipeline='dimple'


    def run(self):
        progress_step=1
        if len(self.sample_list) != 0:
            progress_step=100/float(len(self.sample_list))
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        os.chdir(self.ccp4_scratch_directory)
        os.system('/bin/rm -f xce_%s*sh' %self.pipeline)


        for item in sorted(self.sample_list):

            xtal =                  item[0]
            visit_run_autoproc =    item[1]
            mtzin =                 item[2]
            ref_pdb =               item[3]
            ref_mtz =               item[4]
            ref_cif =               item[5]

            if 'dimple_rerun_on_selected_file' in visit_run_autoproc:
                if self.pipeline=='dimple':
                    self.prepare_dimple_shell_script(xtal,visit_run_autoproc,mtzin,ref_pdb,ref_mtz,ref_cif)
                elif self.pipeline=='pipedream':
                    self.prepare_pipedream_shell_script()
                elif self.pipeline=='phenix.ligand_pipeline':
                    self.prepare_phenix_ligand_pipeline_shell_script()
            else:
                self.prepare_dimple_shell_script(xtal,visit_run_autoproc,mtzin,ref_pdb,ref_mtz,ref_cif)

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.run_script()

        self.emit(QtCore.SIGNAL('datasource_menu_reload_samples'))






    def prepare_phenix_ligand_pipeline_shell_script(self,xtal,visit_run_autoproc,mtzin,ref_pdb,ref_mtz,ref_cif):

        # check if reference mtzfile has an Rfree column; if not, then ignore
        # DIMPLE assumes an Rfree column and barfs if it is not present
        # note: ref_mtz looks like this: ref mtz  -R reference.mtz
#        if os.path.isfile(ref_mtz):
#            mtz_column_dict=mtztools(ref_mtz).get_all_columns_as_dict()
#            if 'FreeR_flag' not in mtz_column_dict['RFREE']:
#                self.Logfile.insert('cannot find FreeR_flag in reference mtz file: %s -> ignoring reference mtzfile!!!' %ref_mtz)
#                ref_mtz = ''
#                if mtz_column_dict['RFREE'] != []:
#                    self.Logfile.insert('found Rfree set with other column name though: %s' %str(mtz_column_dict['RFREE']))
#                    self.Logfile.insert('try renaming Rfree column to FreeR_flag with CAD!')
#
#        db_dict={}
#        db_dict['DimpleReferencePDB']=ref_pdb
#        self.db.update_data_source(xtal,db_dict)
#
#        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'creating input script for '+xtal+' in '+visit_run_autoproc)
#

        if not os.path.isdir(os.path.join(self.initial_model_directory,xtal)):
            os.mkdir(os.path.join(self.initial_model_directory,xtal))
        if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'phenix.ligand_pipeline')):
            os.mkdir(os.path.join(self.initial_model_directory,xtal,'phenix.ligand_pipeline'))
        os.chdir(os.path.join(self.initial_model_directory,xtal,'phenix.ligand_pipeline'))
        os.system('touch dimple_run_in_progress')
        os.system('/bin/rm final.mtz 2> /dev/null')
        os.system('/bin/rm final.pdb 2> /dev/null')

        if self.queueing_system_available:
            top_line='#PBS -joe -N XCE_%s\n' %self.pipeline
        else:
            top_line='#!'+os.getenv('SHELL')+'\n'

        if 'csh' in os.getenv('SHELL'):
            ccp4_scratch='setenv CCP4_SCR '+self.ccp4_scratch_directory+'\n'
        elif 'bash' in os.getenv('SHELL'):
            ccp4_scratch='export CCP4_SCR='+self.ccp4_scratch_directory+'\n'
        else:
            ccp4_scratch=''

        if os.path.isdir('/dls'):
            ccp4_scratch+='module load phenix\n'

        Cmds = (
                '%s\n' %top_line+
                '\n'
                'export XChemExplorer_DIR="'+os.getenv('XChemExplorer_DIR')+'"\n'
                '\n'
                'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'phenix.ligand_pipeline') +
                '\n'
                'source $XChemExplorer_DIR/setup-scripts/xce.setup-sh\n'
                '\n'
                +ccp4_scratch+
                '\n'
                '$CCP4/bin/ccp4-python $XChemExplorer_DIR/helpers/update_status_flag.py %s %s %s %s\n' %(database,xtal,'DimpleStatus','running') +
                '\n'
                'phenix.ligand_pipeline %s %s' %(ref_pdb,mtzin)+
                ' mr=False'
                ' ligand_copies=0'
                ' build=False'
                ' prune=False'
                ' remove_waters=False'
                ' stop_if_r_free_greater_than=0.4'
                ' update_waters=False'
                ' build_hydrogens=False\n'
                '\n'
                'cd %s\n' %os.path.join(self.initial_model_directory,xtal) +
                '\n'
                'ln -s phenix.ligand_pipeline/pipeline_1/refine_final.pdb dimple.pdb'
                'ln -s phenix.ligand_pipeline/pipeline_1/refine_final.mtz dimple.mtz'
                '\n'
                'fft hklin dimple.mtz mapout 2fofc.map << EOF\n'
                ' labin F1=2FOFCWT PHI=PH2FOFCWT\n'
                'EOF\n'
                '\n'
                'fft hklin dimple.mtz mapout fofc.map << EOF\n'
                ' labin F1=FOFCWT PHI=PHFOFCWT\n'
                'EOF\n'
                '\n'
                '$CCP4/libexec/python '+os.path.join(os.getenv('XChemExplorer_DIR'),'helpers','update_data_source_for_new_dimple_pdb.py')+
                ' %s %s %s\n' %(os.path.join(self.database_directory,self.data_source_file),xtal,self.initial_model_directory)+
                '\n'
                '/bin/rm dimple_run_in_progress\n'
                )

        os.chdir(self.ccp4_scratch_directory)
        f = open('xce_%s_%s.sh' %(self.pipeline,str(n+1)),'w')
        f.write(Cmds)
        f.close()
        self.n+=1
        os.system('chmod +x xce_%s_%s.sh' %(self.pipeline,str(n+1)))
        db_dict={}
        db_dict['DimpleStatus']='started'
        self.Logfile.insert('%s: setting DataProcessingStatus flag to started' %xtal)
        self.db.update_data_source(xtal,db_dict)


    def prepare_pipedream_shell_script(self,xtal,visit_run_autoproc,mtzin,ref_pdb,ref_mtz,ref_cif):

        # check if reference mtzfile has an Rfree column; if not, then ignore
        # DIMPLE assumes an Rfree column and barfs if it is not present
        # note: ref_mtz looks like this: ref mtz  -R reference.mtz
#        if os.path.isfile(ref_mtz):
#            mtz_column_dict=mtztools(ref_mtz).get_all_columns_as_dict()
#            if 'FreeR_flag' not in mtz_column_dict['RFREE']:
#                self.Logfile.insert('cannot find FreeR_flag in reference mtz file: %s -> ignoring reference mtzfile!!!' %ref_mtz)
#                ref_mtz = ''
#                if mtz_column_dict['RFREE'] != []:
#                    self.Logfile.insert('found Rfree set with other column name though: %s' %str(mtz_column_dict['RFREE']))
#                    self.Logfile.insert('try renaming Rfree column to FreeR_flag with CAD!')
#
#        db_dict={}
#        db_dict['DimpleReferencePDB']=ref_pdb
#        self.db.update_data_source(xtal,db_dict)
#
#        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'creating input script for '+xtal+' in '+visit_run_autoproc)
#

        if not os.path.isdir(os.path.join(self.initial_model_directory,xtal)):
            os.mkdir(os.path.join(self.initial_model_directory,xtal))
        if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'pipedream')):
            os.mkdir(os.path.join(self.initial_model_directory,xtal,'pipedream'))
        os.chdir(os.path.join(self.initial_model_directory,xtal,'pipedream'))
        os.system('touch dimple_run_in_progress')
        os.system('/bin/rm final.mtz 2> /dev/null')
        os.system('/bin/rm final.pdb 2> /dev/null')

        if self.queueing_system_available:
            top_line='#PBS -joe -N XCE_%s\n' %self.pipeline
        else:
            top_line='#!'+os.getenv('SHELL')+'\n'

        if 'csh' in os.getenv('SHELL'):
            ccp4_scratch='setenv CCP4_SCR '+self.ccp4_scratch_directory+'\n'
        elif 'bash' in os.getenv('SHELL'):
            ccp4_scratch='export CCP4_SCR='+self.ccp4_scratch_directory+'\n'
        else:
            ccp4_scratch=''

        if os.path.isdir('/dls'):
            ccp4_scratch+='module load buster\n'

        if os.path.isfile(ref_mtz):
            hklref_line=' -hklref %s' %ref_mtz
        else:
            hklref_line=''

        Cmds = (
                '%s\n' %top_line+
                '\n'
                'export XChemExplorer_DIR="'+os.getenv('XChemExplorer_DIR')+'"\n'
                '\n'
                'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'pipedream') +
                '\n'
                'source $XChemExplorer_DIR/setup-scripts/xce.setup-sh\n'
                '\n'
                +ccp4_scratch+
                '\n'
                '$CCP4/bin/ccp4-python $XChemExplorer_DIR/helpers/update_status_flag.py %s %s %s %s\n' %(database,xtal,'DimpleStatus','running') +
                '\n'
                'pipedream '
                ' -d pipedreamDir'
                ' -xyzin %s' %ref_pdb+
                hklref_line+
                ' -hklin %s' %mtzin+
                ' -keepwater\n'
                '\n'
                'cd %s\n' %os.path.join(self.initial_model_directory,xtal) +
                '\n'
                'ln -s phenix.ligand_pipeline/pipeline_1/refine_final.pdb dimple.pdb'
                'ln -s phenix.ligand_pipeline/pipeline_1/refine_final.mtz dimple.mtz'
                '\n'
                'fft hklin dimple.mtz mapout 2fofc.map << EOF\n'
                ' labin F1=FWT PHI=PHWT\n'
                'EOF\n'
                '\n'
                'fft hklin dimple.mtz mapout fofc.map << EOF\n'
                ' labin F1=DELFWT PHI=PHDELWT\n'
                'EOF\n'
                '\n'
                '$CCP4/libexec/python '+os.path.join(os.getenv('XChemExplorer_DIR'),'helpers','update_data_source_for_new_dimple_pdb.py')+
                ' %s %s %s\n' %(os.path.join(self.database_directory,self.data_source_file),xtal,self.initial_model_directory)+
                '\n'
                '/bin/rm dimple_run_in_progress\n'
                )

        os.chdir(self.ccp4_scratch_directory)
        f = open('xce_%s_%s.sh' %(self.pipeline,str(n+1)),'w')
        f.write(Cmds)
        f.close()
        self.n+=1
        os.system('chmod +x xce_%s_%s.sh' %(self.pipeline,str(n+1)))
        db_dict={}
        db_dict['DimpleStatus']='started'
        self.Logfile.insert('%s: setting DataProcessingStatus flag to started' %xtal)
        self.db.update_data_source(xtal,db_dict)




#pipedream -d pipedreamDir -xyzin $reference.pdb -hklref $referenceDataset.mtz -hklin $mtzFileUnique.mtz -rhofit $inhib.cif  -target $reference.pdb -keepwater -nthreads 12 -postquick -allclusters |tee -a pipedream.log


    def prepare_dimple_shell_script(self,xtal,visit_run_autoproc,mtzin,ref_pdb,ref_mtz,ref_cif):

        # check if reference mtzfile has an Rfree column; if not, then ignore
        # DIMPLE assumes an Rfree column and barfs if it is not present
        # note: ref_mtz looks like this: ref mtz  -R reference.mtz
        if os.path.isfile(ref_mtz):
            mtz_column_dict=mtztools(ref_mtz).get_all_columns_as_dict()
            if 'FreeR_flag' not in mtz_column_dict['RFREE']:
                self.Logfile.insert('cannot find FreeR_flag in reference mtz file: %s -> ignoring reference mtzfile!!!' %ref_mtz)
                ref_mtz = ''
                if mtz_column_dict['RFREE'] != []:
                    self.Logfile.insert('found Rfree set with other column name though: %s' %str(mtz_column_dict['RFREE']))
                    self.Logfile.insert('try renaming Rfree column to FreeR_flag with CAD!')

        db_dict={}
        db_dict['DimpleReferencePDB']=ref_pdb
        self.db.update_data_source(xtal,db_dict)

        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'creating input script for '+xtal+' in '+visit_run_autoproc)

        if not os.path.isdir(os.path.join(self.initial_model_directory,xtal)):
            os.mkdir(os.path.join(self.initial_model_directory,xtal))
        if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'dimple')):
            os.mkdir(os.path.join(self.initial_model_directory,xtal,'dimple'))
        if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc)):
            os.mkdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc))
        os.chdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc))
        os.system('touch dimple_run_in_progress')
        os.system('/bin/rm final.mtz 2> /dev/null')
        os.system('/bin/rm final.pdb 2> /dev/null')

        if self.queueing_system_available:
            top_line='#PBS -joe -N XCE_dimple\n'
        else:
            top_line='#!'+os.getenv('SHELL')+'\n'

        if 'csh' in os.getenv('SHELL'):
            ccp4_scratch='setenv CCP4_SCR '+self.ccp4_scratch_directory+'\n'
        elif 'bash' in os.getenv('SHELL'):
            ccp4_scratch='export CCP4_SCR='+self.ccp4_scratch_directory+'\n'
        else:
            ccp4_scratch=''

        if 'dimple_rerun_on_selected_file' in visit_run_autoproc:
            additional_cmds = (
                            'cd %s\n' %os.path.join(self.initial_model_directory,xtal) +
                            '/bin/rm dimple.pdb\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/final.pdb dimple.pdb\n'
                            '/bin/rm dimple.mtz\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/final.mtz dimple.mtz\n'
                            '/bin/rm 2fofc.map\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/2fofc.map .\n'
                            '/bin/rm fofc.map\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/fofc.map .\n'
                            '\n'
                            '$CCP4/libexec/python '+os.path.join(os.getenv('XChemExplorer_DIR'),'helpers','update_data_source_for_new_dimple_pdb.py')+
                            ' %s %s %s\n' %(os.path.join(self.database_directory,self.data_source_file),xtal,self.initial_model_directory)  )

        else:
            additional_cmds=''

        Cmds = (
                '%s\n' %top_line+
                '\n'
                'export XChemExplorer_DIR="'+os.getenv('XChemExplorer_DIR')+'"\n'
                '\n'
                'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc) +
                '\n'
                'source $XChemExplorer_DIR/setup-scripts/xce.setup-sh\n'
                '\n'
                +ccp4_scratch+
                '\n'
                '$CCP4/bin/ccp4-python $XChemExplorer_DIR/helpers/update_status_flag.py %s %s %s %s\n' %(database,xtal,'DimpleStatus','running') +
                '\n'
                'dimple --no-cleanup %s %s %s %s dimple\n' %(mtzin,ref_pdb,ref_mtz,ref_cif) +
                '\n'
                'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc,'dimple') +
                '\n'
                'fft hklin final.mtz mapout 2fofc.map << EOF\n'
                ' labin F1=FWT PHI=PHWT\n'
                'EOF\n'
                '\n'
                'fft hklin final.mtz mapout fofc.map << EOF\n'
                ' labin F1=DELFWT PHI=PHDELWT\n'
                'EOF\n'
                '\n'
                +additional_cmds+
                '\n'
                'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc) +
                '\n'
                '/bin/rm dimple_run_in_progress\n'
                '\n'
                'ln -s dimple/final.pdb .\n'
                'ln -s dimple/final.mtz .\n'
                )

        os.chdir(self.ccp4_scratch_directory)
        f = open('xce_%s_%s.sh' %(self.pipeline,str(n+1)),'w')
        f.write(Cmds)
        f.close()
        self.n+=1
        os.system('chmod +x xce_%s_%s.sh' %(self.pipeline,str(n+1)))
        db_dict={}
        db_dict['DimpleStatus']='started'
        self.Logfile.insert('%s: setting DataProcessingStatus flag to started' %xtal)
        self.db.update_data_source(xtal,db_dict)


    def run_script(self):
        # submit job
        self.Logfile.insert('created input scripts for '+str(self.n)+' in '+self.ccp4_scratch_directory)
        os.chdir(self.ccp4_scratch_directory)
        if os.path.isdir('/dls'):
            if self.external_software['qsub_array']:
                Cmds = (
                        '#PBS -joe -N xce_%s_master\n' %self.pipeline+
                        './xce_%s_$SGE_TASK_ID.sh\n' %self.pipeline
                        )
                f = open('%s_master.sh' %self.pipeline,'w')
                f.write(Cmds)
                f.close()
                self.Logfile.insert('submitting array job with maximal 100 jobs running on cluster')
                self.Logfile.insert('using the following command:')
                self.Logfile.insert('qsub -P labxchem -t 1:%s -tc %s %s_master.sh' %(str(self.n),self.max_queue_jobs,self.pipeline))
                os.system('qsub -P labxchem -t 1:%s -tc %s %s_master.sh' %(str(self.n),self.max_queue_jobs,self.pipeline))
            else:
                self.Logfile.insert("cannot start ARRAY job: make sure that 'module load global/cluster' is in your .bashrc or .cshrc file")
        elif self.external_software['qsub']:
            self.Logfile.insert('submitting %s individual jobs to cluster' %(str(self.n)))
            self.Logfile.insert('WARNING: this could potentially lead to a crash...')
            for i in range(self.n):
                self.Logfile.insert('qsub xce_%s_%s.sh' %(str(i+1),self.pipeline))
                os.system('qsub xce_%s_%s.sh' %(str(i+1),self.pipeline))
        else:
            self.Logfile.insert('running %s consecutive %s jobs on your local machine' %(self.pipeline))
            for i in range(self.n):
                self.Logfile.insert('starting xce_%s_%s.sh' %(str(i+1),self.pipeline))
                os.system('./xce_%s_%s.sh' %(str(i+1),self.pipeline))



class remove_selected_dimple_files(QtCore.QThread):
    def __init__(self,sample_list,initial_model_directory,xce_logfile,database_directory,data_source_file,):
        QtCore.QThread.__init__(self)
        self.sample_list=sample_list
        self.initial_model_directory=initial_model_directory
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.db=XChemDB.data_source(os.path.join(database_directory,data_source_file))

    def run(self):
        progress_step=1
        if len(self.sample_list) != 0:
            progress_step=100/float(len(self.sample_list))
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        for n,xtal in enumerate(self.sample_list):
            db_dict={}
            os.chdir(os.path.join(self.initial_model_directory,xtal))
            self.Logfile.insert('%s: removing dimple.pdb/dimple.mtz' %xtal)
            os.system('/bin/rm dimple.pdb 2> /dev/null')
            os.system('/bin/rm dimple.mtz 2> /dev/null')
            if os.path.isdir(os.path.join(self.initial_model_directory,xtal,'dimple','dimple_rerun_on_selected_file')):
                os.chdir('dimple')
                self.Logfile.insert('%s removing directory dimple/dimple_rerun_on_selected_file' %xtal)
                os.system('/bin/rm -fr dimple_rerun_on_selected_file')

            db_dict['DimpleResolutionHigh']=''
            db_dict['DimpleRcryst']=''
            db_dict['DimpleRfree']=''
            db_dict['DimplePathToPDB']=''
            db_dict['DimplePathToMTZ']=''
            db_dict['DimpleReferencePDB']=''
            db_dict['DimplePANDDAwasRun']='False'
            db_dict['DimplePANDDAhit']='False'
            db_dict['DimplePANDDAreject']='False'
            db_dict['DimplePANDDApath']=''
            db_dict['DimpleStatus']='pending'

            self.Logfile.insert('%s: updating database' %xtal)
            self.db.update_data_source(xtal,db_dict)


            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.emit(QtCore.SIGNAL('datasource_menu_reload_samples'))




class start_COOT(QtCore.QThread):

    def __init__(self,settings,interface):
        QtCore.QThread.__init__(self)
        self.settings=settings
        if interface=='old':
            self.pylib='XChemCoot.py'
        elif interface=='new':
            self.pylib='XChemCootNew.py'

    def run(self):
        cwd=os.getcwd()
        # coot at Diamond always or sometimes at least open in home directory, so then it won't find the .pkl file
        pickle.dump(self.settings,open(os.path.join(os.getenv('HOME'),'.xce_settings.pkl'),'wb'))
        os.system('cd %s\ncoot --no-guano --no-state-script --script %s' %(os.getenv('HOME'),os.path.join(os.getenv('XChemExplorer_DIR'),'lib',self.pylib)))



class start_ICM(QtCore.QThread):

    def __init__(self,html_export_directory):
        QtCore.QThread.__init__(self)
        self.html_export_directory=html_export_directory

    def run(self):
        cwd=os.getcwd()
        if cwd.startswith('/dls'):
#            if not os.path.isfile(os.path.join('/home',getpass.getuser(),'.flexlmrc')):
#                os.system('touch '+os.path.join('/home',getpass.getuser(),'.flexlmrc'))
#                f=open(os.path.join('/home',getpass.getuser(),'.flexlmrc'),'w')
#                f.write('MOLSOFTD_LICENSE_FILE=@diamvicmpro.diamond.ac.uk')
#                f.close()
            os.system('nautilus %s &' %self.html_export_directory)
            os.system('/dls/science/groups/i04-1/software/icm-3.8-5/icm64 -g')

class start_pandda_inspect(QtCore.QThread):

    def __init__(self,settings,xce_logfile):
        QtCore.QThread.__init__(self)
#        self.settings=settings
        self.panddas_directory=settings['panddas_directory']
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)

    def run(self):
        if os.getenv('SHELL') == '/bin/tcsh' or os.getenv('SHELL') == '/bin/csh':
            source_file=os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-csh')
        elif os.getenv('SHELL') == '/bin/bash':
            source_file=os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-sh')
        else:
            source_file=''

        Cmds = (
                '#!'+os.getenv('SHELL')+'\n'
                'unset PYTHONPATH\n'
                'source '+source_file+'\n'
                'cd '+self.panddas_directory+'\n'
                'pandda.inspect\n'
            )

        self.Logfile.insert('starting pandda.inspect with the following command:\n'+Cmds)
        os.system(Cmds)

#        Cmds = (
#                'source '+os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-sh')+'\n'
#                'cd '+self.panddas_directory+'\n'
#                'pandda.inspect\n'
#            )
#
#        self.Logfile.insert('starting pandda.inspect with the following command:\n/bin/bash\n'+Cmds)
#        # need to do this because we're having all sorts of csh bash issues at SGC
#        os.system('/bin/bash\n'+Cmds)





class start_dials_image_viewer(QtCore.QThread):

    def __init__(self,diffraction_image):
        QtCore.QThread.__init__(self)
        self.diffraction_image=diffraction_image

    def run(self):
        os.system('dials.image_viewer '+self.diffraction_image)




class save_autoprocessing_results_to_disc(QtCore.QThread):
    def __init__(self,dataset_outcome_dict,
                      data_collection_table_dict,
                      data_collection_column_three_dict,
                      data_collection_dict,
                      database_directory,data_source_file,
                      initial_model_directory,
                      preferences,
                      data_collection_summary_file):
        QtCore.QThread.__init__(self)
        self.dataset_outcome_dict=dataset_outcome_dict
        self.data_collection_table_dict=data_collection_table_dict
        self.data_collection_column_three_dict=data_collection_column_three_dict
        self.data_collection_dict=data_collection_dict
        self.database_directory=database_directory
        self.data_source_file=data_source_file
        self.initial_model_directory=initial_model_directory
        self.processed_data_to_copy=preferences['processed_data_to_copy']
        self.data_collection_summary_file=data_collection_summary_file

    def run(self):

        if not len(self.dataset_outcome_dict)==0:
            progress_step=100/float(len(self.dataset_outcome_dict))
        progress=0

        data_source=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))

        # previously we did first update the data source so that it reflects the latest state of selections;
        # this is not necessary anymore, because the XCE now updates the data source whenever something changes;
        # the only thing we need to change for a selected sample is the path and name entries;
        # which also needs to be done in the pkl file!
        # another thing that is new is that we copy now files of ALL auto-processing runs;
        # but do not run ctruncate anymore

        ########################################################
        # 1. copy files
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
        for sample in sorted(self.data_collection_dict):
            # find out which row was selected in respective data collection table
            selected_processing_result='n/a'
            indexes=self.data_collection_column_three_dict[sample][0].selectionModel().selectedRows()
            if indexes != []:       # i.e. logfile exists
                for index in sorted(indexes):
                    selected_processing_result=index.row()

            for n,entry in enumerate(self.data_collection_dict[sample]):
                if entry[0]=='logfile':
                    visit=entry[1]
                    run=entry[2]
                    autoproc=entry[4]
                    db_dict=entry[6]
                    outcome=self.dataset_outcome_dict[sample]
                    if outcome.startswith('Failed'):
                        continue
                    self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'writing all files from data processing to project folder -> '+sample+', visit: '+visit+', run: '+run+', program: '+autoproc)
                    path_to_procdir=db_dict['DataProcessingDirectoryOriginal']
                    path_to_logfile=db_dict['DataProcessingPathToLogfile']
                    try:
                        path_to_mtzfile=db_dict['DataProcessingPathToMTZfile']
                    except KeyError:
                        continue
                    if path_to_mtzfile=='':     # in case a log file exists, but the job did not produce an mtz file
                        continue
                    mtz_filename=db_dict['DataProcessingMTZfileName']
                    log_filename=db_dict['DataProcessingLOGfileName']
                    dimple_destination=''
                    try:
                        path_to_dimple_pdbfile=db_dict['DataProcessingPathToDimplePDBfile']
                        path_to_dimple_mtzfile=db_dict['DataProcessingPathToDimpleMTZfile']
                    except KeyError:
                        path_to_dimple_pdbfile=''
                        path_to_dimple_mtzfile=''

                    # create all the directories if necessary
                    if not os.path.isdir(os.path.join(self.initial_model_directory,sample)):
                        os.mkdir(os.path.join(self.initial_model_directory,sample))
                    if not os.path.isdir(os.path.join(self.initial_model_directory,sample,'autoprocessing')):
                        os.mkdir(os.path.join(self.initial_model_directory,sample,'autoprocessing'))
                    if not os.path.isdir(os.path.join(self.initial_model_directory,sample,'autoprocessing',visit+'-'+run+autoproc)):
                        os.mkdir(os.path.join(self.initial_model_directory,sample,'autoprocessing',visit+'-'+run+autoproc))

                    if path_to_dimple_pdbfile != '':
                        if not os.path.isdir(os.path.join(self.initial_model_directory,sample,'dimple')):
                            os.mkdir(os.path.join(self.initial_model_directory,sample,'dimple'))
                        if not os.path.isdir(os.path.join(self.initial_model_directory,sample,'dimple',visit+'-'+run+autoproc)):
                            os.mkdir(os.path.join(self.initial_model_directory,sample,'dimple',visit+'-'+run+autoproc))
                        dimple_destination=os.path.join(self.initial_model_directory,sample,'dimple',visit+'-'+run+autoproc)

                    if self.processed_data_to_copy=='mtz_log_only':
                        path_to_logfile,path_to_mtzfile,mtz_filename=self.copy_mtz_and_logfiles_only(sample,autoproc,run,visit,path_to_procdir,path_to_logfile,path_to_mtzfile,mtz_filename,log_filename)
                        db_dict['DataProcessingPathToLogfile']=os.path.join(self.initial_model_directory,sample,path_to_logfile)
                        db_dict['DataProcessingPathToMTZfile']=os.path.join(self.initial_model_directory,sample,path_to_mtzfile)
                        self.copy_and_link_selected_dimple_files(dimple_destination,sample,path_to_dimple_mtzfile,path_to_dimple_pdbfile)
                        db_dict['DataProcessingPathToDimplePDBfile']=dimple_destination
                        db_dict['DataProcessingPathToDimpleMTZfile']=dimple_destination

                    # update pkl file
                    entry[6]=db_dict
                    self.data_collection_dict[sample][n]=entry

#                    elif self.processed_data_to_copy=='everything':
#                        path_to_logfile,path_to_mtzfile,mtz_filename,log_filename=self.copy_complete_autoprocessing_folder(sample,autoproc,run,visit,path_to_procdir,path_to_logfile,path_to_mtzfile,mtz_filename,log_filename)
#                        self.link_mtz_log_files_to_sample_directory(sample,autoproc,run,visit,path_to_procdir,path_to_logfile,path_to_mtzfile,mtz_filename,log_filename)
#                        self.copy_and_link_selected_dimple_files(dimple_destination,sample,path_to_dimple_mtzfile,path_to_dimple_pdbfile)

                    # update data source if this is the selected file
                    # and make respective links
                    if entry[7]==selected_processing_result:
                        db_dict=entry[6]
                        db_dict['DataCollectionOutcome']=self.dataset_outcome_dict[sample]
                        db_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
                        tmp_dict={}
                        tmp_dict['RefinementMTZfree'],tmp_dict['DimpleRcryst'],tmp_dict['DimpleRfree'],tmp_dict['RefinementOutcome'],tmp_dict['RefinementSpaceGroup'] =self.link_mtz_log_files_to_sample_directory(sample,autoproc,run,visit,path_to_procdir,path_to_logfile,path_to_mtzfile,mtz_filename,log_filename,dimple_destination)
                        if tmp_dict['RefinementOutcome'] != '':
                            # this assumes that no manual DIMPLE run was launched successfully
                            db_dict['RefinementMTZfree'] = tmp_dict['RefinementMTZfree']
                            db_dict['DimpleRcryst'] = tmp_dict['DimpleRcryst']
                            db_dict['DimpleRfree'] = tmp_dict['DimpleRfree']
                            db_dict['RefinementOutcome'] = tmp_dict['RefinementOutcome']
                            db_dict['RefinementSpaceGroup'] = tmp_dict['RefinementSpaceGroup']
                        current_refinement_outcome=data_source.get_value_from_field(sample,'RefinementOutcome')
                        if str(current_refinement_outcome[0]).split()[0].lower().startswith('none'):
                            db_dict['RefinementOutcome']='1 - Analysis Pending'
                        data_source.update_insert_data_source(sample,db_dict)


            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'saving updated pkl file')
        pickle.dump(self.data_collection_dict,open(self.data_collection_summary_file,'wb'))

        self.emit(QtCore.SIGNAL("finished()"))

    def copy_and_link_selected_dimple_files(self,dimple_destination,sample,path_to_dimple_mtzfile,path_to_dimple_pdbfile):
        # dimple files
        if dimple_destination != '':
            os.chdir(dimple_destination)
            if not os.path.isfile('final.pdb'):
                os.system('/bin/cp '+path_to_dimple_pdbfile+' .')
            if not os.path.isfile('final.mtz'):
                os.system('/bin/cp '+path_to_dimple_mtzfile+' .')


    def link_mtz_log_files_to_sample_directory(self,sample,autoproc,run,visit,path_to_procdir,path_to_logfile,path_to_mtzfile,mtz_filename,log_filename,dimple_destination):
        mtzfree=''
        Rcryst=''
        Rfree=''
        refinement_stage=''
        spg=''
        relative_dimple_destination='.'+dimple_destination.replace(os.path.join(self.initial_model_directory,sample),'')
        relative_path_to_mtzfile='./'+path_to_mtzfile.replace(self.initial_model_directory,'')
        relative_path_to_logfile='./'+path_to_logfile.replace(self.initial_model_directory,'')
        # move up to sample directory and link respective files
        # first remove any old symbolic links
        os.chdir(os.path.join(self.initial_model_directory,sample))
        if os.path.islink(os.path.join(self.initial_model_directory,sample,sample+'.mtz')):
            os.system('/bin/rm '+os.path.join(self.initial_model_directory,sample,sample+'.mtz'))
        if os.path.islink(os.path.join(self.initial_model_directory,sample,sample+'.log')):
            os.system('/bin/rm '+os.path.join(self.initial_model_directory,sample,sample+'.log'))
        # then link new files
#        os.symlink(os.path.join(path_to_mtzfile,sample+'.mtz'),sample+'.mtz')
#        os.symlink(os.path.join(path_to_logfile,sample+'.log'),sample+'.log')
        os.symlink(os.path.join(relative_path_to_mtzfile,sample+'.mtz'),sample+'.mtz')
        os.symlink(os.path.join(relative_path_to_logfile,sample+'.log'),sample+'.log')

        # only continue changing DIMPLE links if it was not yet run successfully an selected mtz files
        if not os.path.isfile(os.path.join(self.initial_model_directory,sample,'dimple','dimple_rerun_on_selected_file','dimple','final.pdb')):

            if dimple_destination != '':
                if os.path.isfile(os.path.join(dimple_destination,'final.pdb')):
                    # remove old symbolic links if necessary
                    #if os.path.isfile('dimple.pdb'):
                    # forget about the if statement; it won't catch broken links,
                    # but broken links will still trip os.symlink
                    os.system('/bin/rm dimple.pdb 2> /dev/null')
                    # symlink with absolute path
#                    os.symlink(os.path.join(dimple_destination,'final.pdb'),'dimple.pdb')
                    # symlink with relative paths
                    os.symlink(os.path.join(relative_dimple_destination,'final.pdb'),'dimple.pdb')
                    pdb_info=parse().PDBheader(os.path.join(dimple_destination,'final.pdb'))
                    Rcryst=pdb_info['Rcryst']
                    Rfree=pdb_info['Rfree']
                    spg=pdb_info['SpaceGroup']
                    refinement_stage='1 - Analysis Pending'

                if os.path.isfile(os.path.join(dimple_destination,'final.mtz')):
                    # remove old symbolic links if necessary
                    #if os.path.isfile('dimple.mtz'):
                    os.system('/bin/rm dimple.mtz 2> /dev/null')
#                   os.symlink(os.path.join(dimple_destination,'final.mtz'),'dimple.mtz')
                    os.symlink(os.path.join(relative_dimple_destination,'final.mtz'),'dimple.mtz')
                # if no refinement was carried out yet, then we also want to link the dimple files to refine.pdb/refine.log
                # so that we can look at them with the COOT plugin
                found_previous_refinement=False
                for dirs in glob.glob('*'):
                    if os.path.isdir(dirs) and dirs.startswith('Refine_'):
                        found_previous_refinement=True
                        break
                if not found_previous_refinement:
                    # first delete possible old symbolic links
                    #if os.path.isfile('refine.pdb'):
                    os.system('/bin/rm refine.pdb 2> /dev/null')
#                    os.symlink('dimple.pdb','refine.pdb')
                    #if os.path.isfile('refine.mtz'):
                    os.system('/bin/rm refine.mtz 2> /dev/null')
#                    os.symlink('dimple.mtz','refine.mtz')
                # remove any previous <sample>.free.mtz file, and link new dimple.mtz
                # so if we continue refining, then we do so against the correct file
                # think that REFMAC does not tinker with F,SIGF as long as there is no twinning
                #if os.path.isfile(sample+'.free.mtz'):
                os.system('/bin/rm '+sample+'.free.mtz 2> /dev/null')
#                os.symlink(os.path.join(dimple_destination,'final.mtz'),sample+'.free.mtz')
                os.symlink(os.path.join(relative_dimple_destination,'final.mtz'),sample+'.free.mtz')
                mtzfree=os.path.join(dimple_destination,'final.mtz')
        return mtzfree,Rcryst,Rfree,refinement_stage,spg


    def copy_mtz_and_logfiles_only(self,sample,autoproc,run,visit,path_to_procdir,path_to_logfile,path_to_mtzfile,mtz_filename,log_filename):
        os.chdir(os.path.join(self.initial_model_directory,sample,'autoprocessing',visit+'-'+run+autoproc))
        # don't do anything if file already exists
        if not os.path.isfile(mtz_filename):
            os.system('/bin/cp '+path_to_mtzfile+' .')
        if not os.path.isfile(log_filename):
            os.system('/bin/cp '+path_to_logfile+' .')
        # filenames may be rather long and cryptic, but we link them to <sample>.mtz, <sample>.log;
        # won't have to worry later what they're called; even though info is stored in pkl file
        if not os.path.isfile(sample+'.mtz'):
            os.symlink(mtz_filename,sample+'.mtz')
        if not os.path.isfile(sample+'.log'):
            os.symlink(log_filename,sample+'.log')
        # in case the user copied the results from several data processing pipelines and just wants to
        # set the current one
        path_to_logfile=os.path.join('autoprocessing',visit+'-'+run+autoproc)
        path_to_mtzfile=os.path.join('autoprocessing',visit+'-'+run+autoproc)
        return path_to_logfile,path_to_mtzfile,mtz_filename

#    def copy_complete_autoprocessing_folder(self,sample,autoproc,run,visit,path_to_procdir,path_to_logfile,path_to_mtzfile,mtz_filename,log_filename):
#        os.chdir(os.path.join(self.initial_model_directory,sample,'autoprocessing',visit+'-'+run+autoproc))
#        # in this case, ignore if directory already exists
#        if not os.path.isdir(autoproc):
#            os.system('/bin/cp -Rf '+path_to_procdir+' .')
#            if 'xia2' in path_to_logfile:
##                path_to_logfile=os.path.join('autoprocessing',visit+'-'+run+autoproc)+'/'+ '/'.join(path_to_logfile.split('/')[len(path_to_logfile.split('/'))-3:len(path_to_logfile.split('/'))-1])
##                path_to_mtzfile=os.path.join('autoprocessing',visit+'-'+run+autoproc)+'/'+ '/'.join(path_to_mtzfile.split('/')[len(path_to_mtzfile.split('/'))-3:len(path_to_mtzfile.split('/'))-1])
#                path_to_logfile='./'+'/'.join(path_to_logfile.split('/')[len(path_to_logfile.split('/'))-3:len(path_to_logfile.split('/'))-1])
#                path_to_mtzfile='./'+'/'.join(path_to_mtzfile.split('/')[len(path_to_mtzfile.split('/'))-3:len(path_to_mtzfile.split('/'))-1])
#            elif 'fast_dp' in path_to_logfile:
#                os.chdir('fast_dp')
#                self.run_ctruncate(sample)
#                os.chdir(os.path.join(self.initial_model_directory,sample,'autoprocessing',visit+'-'+run+autoproc))
##                path_to_logfile=os.path.join('autoprocessing',visit+'-'+run+autoproc,'fast_dp')
##                path_to_mtzfile=os.path.join('autoprocessing',visit+'-'+run+autoproc,'fast_dp')
#                path_to_logfile=os.path.join('./','fast_dp')
#                path_to_mtzfile=os.path.join('./','fast_dp')
#                mtz_filename='ctruncate.mtz'
#            elif 'autoPROC' in path_to_logfile:
#                path_to_logfile=os.path.join('./','autoPROC','ap-run')
#                path_to_mtzfile=os.path.join('./','autoPROC','ap-run')
#
#            os.symlink(os.path.join(path_to_mtzfile,mtz_filename),sample+'.mtz')
#            os.symlink(os.path.join(path_to_logfile,log_filename),sample+'.log')
#
#        # since all mtz/log files are already linked  as <sample>.mtz/log in visit+'-'+run+autoproc directory
#        path_to_mtzfile=os.path.join(self.initial_model_directory,sample,'autoprocessing',visit+'-'+run+autoproc)
#        mtz_filename=sample+'.mtz'
#        path_to_logfile=os.path.join(self.initial_model_directory,sample,'autoprocessing',visit+'-'+run+autoproc)
#        log_filename=sample+'.log'
#
#        return path_to_logfile,path_to_mtzfile,mtz_filename,log_filename





class read_autoprocessing_results_from_disc(QtCore.QThread):
    def __init__(self,visit_list,
                 target,
                 reference_file_list,
                 database_directory,
                 data_collection_dict,
                 preferences,
                 data_collection_summary_file,
                 initial_model_directory,
                 rescore_only,
                 acceptable_low_resolution_limit_for_data,
                 data_source_file,
                 xce_logfile):
        QtCore.QThread.__init__(self)
        self.visit_list=visit_list
        self.target=target
        self.reference_file_list=reference_file_list
        self.data_collection_dict=data_collection_dict
        self.database_directory=database_directory
        self.selection_mechanism=preferences['dataset_selection_mechanism']
        self.data_collection_summary_file=data_collection_summary_file
        self.initial_model_directory=initial_model_directory
        self.rescore_only=rescore_only
        self.acceptable_low_resolution_limit_for_data=acceptable_low_resolution_limit_for_data
        self.data_source=XChemDB.data_source(os.path.join(data_source_file))
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)
#        self.gda_log_directories_parsed=gda_log_directories_parsed

        # - open data source if possible
        # - get sampleID, xtbm
        # - search lab36 folder for respective xtal image
        # - convert to string and use in data dict
        # - but images can only be found of XCE is started in the respective labchem directory

    def run(self):
        if self.target=='=== SELECT TARGET ===':
            self.visit_list=[]
#            print '==> XCE: please select a target first'
#            return

        if not self.data_collection_summary_file.endswith('.pkl'):
            print '==> XCE: please assign new Summary File or select an existing one'
            return

        if self.rescore_only:
            self.rescore_and_reset_pkl_file()
        else:
#            self.parse_file_system()
            self.parse_file_system_NEW()
#            self.parse_challenging_file_system()

    def max_IsigI_Completeness_Reflections(self,xtal):
        # before creating the table with the results, try to guess which one to select
        # 1. check if there are reference mtz files
        # 1a. if so: take all logfiles forward that fit to the first one found
        #     'fit means': same lattice and delta Vunitcell < 5%
        # 2. if possible: select all datasets with Rmerge low < 5%
        # 3. finally select the dataset with
        #    max(unique_reflections*completeness*Mn(I/sig<I>)

        ############################################################################################
        # STAGE 1:
        # similarity to reference files
        select_stage_one_list = []
        tmp=[]
        for n,entry in enumerate(self.data_collection_dict[xtal]):
            if xtal=='CAMK1DA-x0313': print entry
            found=False
            if entry[0]=='logfile':
                index=self.data_collection_dict[xtal][n][7]
                if isinstance(entry[6],dict):
                    try:
                        if isinstance(float(entry[6]['DataProcessingUnitCellVolume']),float):
                            for reference_file in self.reference_file_list:
                                if not reference_file[4]==0:
                                    unitcell_difference=round((math.fabs(reference_file[4]-float(entry[6]['DataProcessingUnitCellVolume']))/reference_file[4])*100,1)
                                    if xtal=='CAMK1DA-x0313': print 'uc diff',unitcell_difference
                                    if xtal=='CAMK1DA-x0313': print 'lattice xtal',entry[6]['DataProcessingLattice']
                                    if xtal=='CAMK1DA-x0313': print 'lattice ref',reference_file[3]
                                    if unitcell_difference < 5 and reference_file[3]==entry[6]['DataProcessingLattice']:
                                        select_stage_one_list.append(index)
                                        found=True
                    except ValueError:
                        pass
                if xtal=='CAMK1DA-x0313': print 'FOUND?',found
                if not found:
                    tmp.append(index)               # so that if no file passes criterion above
                                                    # or if no reference is given, we still carry over all existing files

        # if none passed Stage 1, carry them over to Stage 2
        if select_stage_one_list == [] and tmp != []:
            select_stage_one_list=tmp


        ############################################################################################
        # STAGE 2:
        # if possible, select only the ones with Rmerge < 10%
        select_stage_two_list=[]
        tmp=[]
        for index in select_stage_one_list:
            # this may be completely over the top!
            found=False
            for entry in self.data_collection_dict[xtal]:
                if entry[0]=='logfile':
                    if isinstance(entry[6],dict):
                        try:
                            if float(entry[6]['DataProcessingRmergeLow']) < 0.1 and entry[7]==index:
                                select_stage_two_list.append(index)
                                found=True
                        except ValueError:
                            pass
            if not found:
                tmp.append(index)

        # if none passed Stage 2, carry them over to Stage 3
        if select_stage_two_list == [] and tmp != []:
            select_stage_two_list=tmp

        ############################################################################################
        # STAGE 3:
        # finally, select the file with the highest
        # max(unique_reflections*completeness*Mn(I/sig<I>)
#        select_stage_three_list=[]
#        for index in select_stage_two_list:
#            for n,entry in enumerate(self.data_collection_dict[xtal]):
#                if entry[0]=='logfile':
#                    if isinstance(entry[6],dict) and entry[7]==index:
#                        ranking=entry[6]['DataProcessingScore']
#                        select_stage_three_list.append([index,ranking])
#        if not select_stage_three_list==[]:
#            self.set_best_file_to_true(xtal,'max',select_stage_three_list)

        ############################################################################################
        # select the file with the highest DataProcessingScore
        select_stage_three_list=[]
        for index in select_stage_two_list:
            for n,entry in enumerate(self.data_collection_dict[xtal]):
                if entry[0]=='logfile':
                    if isinstance(entry[6],dict) and entry[7]==index:
                        try:
                            ranking=entry[6]['DataProcessingScore']
#                            print xtal
#                        print entry
                            if isinstance(ranking,float):
                                select_stage_three_list.append([index,ranking])
                        except KeyError:
                            try:
                                ranking = (float(entry[6]['DataProcessingUniqueReflectionsOverall'])*\
                                           float(entry[6]['DataProcessingCompletenessOverall'])*\
                                           float(entry[6]['DataProcessingIsigOverall']))/float(entry[6]['DataProcessingUnitCellVolume'])
                                select_stage_three_list.append([index,ranking])
                            except ValueError:
                                continue
        if not select_stage_three_list==[]:
                self.set_best_file_to_true(xtal,'max',select_stage_three_list)



    def max_IsigI_Completeness_Reflections_only(self,xtal):
        # => select the dataset with
        #    max(unique_reflections*completeness*Mn(I/sig<I>)

        ############################################################################################
        # select the file with the highest DataProcessingScore
        select_stage_three_list=[]
        for n,entry in enumerate(self.data_collection_dict[xtal]):
            if entry[0]=='logfile':
                index=self.data_collection_dict[xtal][n][7]
                if isinstance(entry[6],dict):
                    try:
                        ranking=entry[6]['DataProcessingScore']
#                        print xtal
#                        print entry
                        if isinstance(ranking,float):
                            select_stage_three_list.append([index,ranking])
                    except KeyError:
                        try:
                            ranking = (float(entry[6]['DataProcessingUniqueReflectionsOverall'])*\
                                       float(entry[6]['DataProcessingCompletenessOverall'])*\
                                       float(entry[6]['DataProcessingIsigOverall']))/float(entry[6]['DataProcessingUnitCellVolume'])
                            select_stage_three_list.append([index,ranking])
                        except ValueError:
                            continue

        if not select_stage_three_list==[]:
            self.set_best_file_to_true(xtal,'max',select_stage_three_list)


    def set_best_file_to_true(self,xtal,min_max,input_list):
        if min_max=='min':
            best_file_index=min(input_list,key=lambda x: x[1])[0]
        elif min_max=='max':
            best_file_index=max(input_list,key=lambda x: x[1])[0]
        for n,entry in enumerate(self.data_collection_dict[xtal]):
            if entry[0]=='logfile':
                if entry[7]==best_file_index:
                    self.data_collection_dict[xtal][n][8]=True
                    # if this was just a rescoring excersise, the files are already in the project directory
                    # hence we want all the links to be reset immediately
                    visit=entry[1]
                    run=entry[2]
                    autoproc=entry[4]
                    db_dict=entry[6]
                    try:
                        path_to_logfile=db_dict['DataProcessingPathToLogfile']
                        path_to_mtzfile=db_dict['DataProcessingPathToMTZfile']
                        mtz_filename=db_dict['DataProcessingMTZfileName']
                        log_filename=db_dict['DataProcessingLOGfileName']
                        # first check if folders and files exist
                        # since user might do this before data are actually copied over
                        if os.path.isdir(os.path.join(self.initial_model_directory,xtal,'autoprocessing',visit+'-'+run+autoproc)):
                            db_dict['DataProcessingAutoAssigned']='False'
                            os.chdir(os.path.join(self.initial_model_directory,xtal))
                            # first remove old links
                            os.system('/bin/rm '+xtal+'.mtz')
                            os.system('/bin/rm '+xtal+'.log')
                            # make new links
                            os.symlink(os.path.join(path_to_logfile,log_filename),xtal+'.log')
                            os.symlink(os.path.join(path_to_mtzfile,mtz_filename),xtal+'.mtz')
                    # it happened that XCE gave a KeyError because DataProcessingMTZfileName does not exist
                    # i.e. there was an autoPROC run with an aimless logfile and stats in it and also an aimless.mtz file
                    # but truncate-unique.mtz did not exist because the file was really rubbish and so XCE
                    # did not update the respective key
                    except KeyError:
                        break
                else:
                    self.data_collection_dict[xtal][n][8]=False


    def min_Rfree(self,xtal):
        tmp=[]
        for entry in self.data_collection_dict[xtal]:
            if entry[0]=='logfile':
                db_dict=entry[6]
                index=entry[7]
                try:
                    tmp.append([index, float(db_dict['DataProcessingRfree']) ] )
                except ValueError:
                    pass
        if tmp != []:
            self.set_best_file_to_true(xtal,'min',tmp)
        else:
            self.min_resolution(xtal)



    def min_resolution(self,xtal):
        tmp=[]
        for entry in self.data_collection_dict[xtal]:
            if entry[0]=='logfile':
                db_dict=entry[6]
                index=entry[7]
                try:
                    tmp.append([index, float(db_dict['DataProcessingResolutionHigh']) ] )
                except ValueError:
                    pass
        if tmp != []:
            self.set_best_file_to_true(xtal,'min',tmp)



    def select_best_dataset(self):
        if not len(self.data_collection_dict)==0:
            progress_step=100/float(len(self.data_collection_dict))
        else:
            progress_step=1
        progress=0

        for xtal in sorted(self.data_collection_dict):
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'Step 2 of 2: selecting "best" aimless logfile ->'+xtal)
            # overwrite previous selection, only if flag not present
            overwrite_previous_selection=True
            for entry in self.data_collection_dict[xtal]:
                if entry[0]=='user_changed_selection':
                    overwrite_previous_selection=False
                    break
            if overwrite_previous_selection:
                for n,entry in enumerate(self.data_collection_dict[xtal]):
                    if entry[0]=='logfile':
                        self.data_collection_dict[xtal][n][8]=False
                if self.selection_mechanism=='IsigI*Comp*UniqueRefl':
#                    self.max_IsigI_Completeness_Reflections_only(xtal)
                    self.max_IsigI_Completeness_Reflections(xtal)
                elif self.selection_mechanism=='lowest_Rfree':
                    self.min_Rfree(xtal)
                elif self.selection_mechanism=='highest_resolution':
                    self.min_resolution(xtal)
            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

    def update_datasource(self):
        if not len(self.data_collection_dict)==0:
            progress_step=100/float(len(self.data_collection_dict))
        else:
            progress_step=1
        progress=0

        # get all samples that are currently in DB
        existing_samples=self.data_source.get_all_samples_in_data_source_as_list()

        for sample in sorted(self.data_collection_dict):
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'updating data source for '+sample)
            logfile_found=False
            user_changed_selection=False
            db_dict={}
            db_dict['DataCollectionOutcome']='Failed - unknown'
            tmpList=[]
            for entry in self.data_collection_dict[sample]:
#                if sample == 'SV3CP-x0209':
#                    print entry
                if entry[0]=='user_changed_selection':
                    user_changed_selection=True
                if entry[0]=='logfile':
                    if entry[8]:        # the best auto-selected or user selected output
                        db_dict=entry[6]
                        logfile_found=True
#                        if sample == 'SV3CP-x0209': print db_dict['DataProcessingResolutionHigh']
                        try:
                            if float(db_dict['DataProcessingResolutionHigh']) <= float(self.acceptable_low_resolution_limit_for_data):
#                                if sample == 'SV3CP-x0209': print 'here'
                                db_dict['DataCollectionOutcome']='success'
                            else:
                                db_dict['DataCollectionOutcome']='Failed - low resolution'
                        except ValueError:
                            db_dict['DataCollectionOutcome']='Failed - unknown'
                        db_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
#                        db_dict['DataProcessingAutoAssigned']='True'
                if entry[0]=='image':
                    if len(entry) >= 9:     # need this because some older pkl files won't have the beamline added
                        tmpList.append([datetime.strptime(entry[3], '%Y-%m-%d %H:%M:%S'),entry[1],entry[2],entry[8]])
                    else:
                        tmpList.append([datetime.strptime(entry[3], '%Y-%m-%d %H:%M:%S'),entry[1],entry[2],'I04-1'])


            if not logfile_found:
                if tmpList != []:
                    # find latest data collection in tmpList
                    latest_run=max(tmpList,key=lambda x: x[0])
                    db_dict={   'DataCollectionVisit':              latest_run[1],
                                'DataCollectionBeamline':           latest_run[3],
                                'DataCollectionDate':               latest_run[0],
                                'DataCollectionOutcome':            'Failed - unknown',
                                'RefinementOutcome':                '-1 - Data Collection Failed'}
                else:
                    db_dict={   'DataCollectionVisit':              'unknown',
                                'DataCollectionBeamline':           'unknown',
                                'DataCollectionDate':               datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                'DataCollectionOutcome':            'Failed - unknown',
                                'RefinementOutcome':                '-1 - Data Collection Failed'}
            else:
                for n,entry in enumerate(self.data_collection_dict[sample]):
                    if entry[0]=='logfile':
                        db_data_collection_dict=entry[6]
                        if user_changed_selection:
                            db_dict['DataProcessingAutoAssigned']='False'
                            db_data_collection_dict['DataProcessingAutoAssigned']='False'
                        else:
                            db_dict['DataProcessingAutoAssigned']='True'
                            db_data_collection_dict['DataProcessingAutoAssigned']='True'
                        entry[6]=db_data_collection_dict
                        self.data_collection_dict[sample][n]=entry

#            if sample == 'SV3CP-x0209': print db_dict
            if self.rescore_only:
                self.data_source.update_insert_data_source(sample,db_dict)
            elif user_changed_selection==False:     # if user changed the selection, then ignore
                self.data_source.update_insert_data_source(sample,db_dict)
            elif sample not in existing_samples:
#                if sample == 'SV3CP-x0209': print 'hallo'
                self.data_source.update_insert_data_source(sample,db_dict)
            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)


    def rescore_and_reset_pkl_file(self):
        # remove 'user_changed_selection' flag from dictionary
        for xtal in self.data_collection_dict:
            for entry in self.data_collection_dict[xtal]:
                if entry[0]=='user_changed_selection':
                    self.data_collection_dict[xtal].remove(entry)

        # and now select again the best dataset
        self.select_best_dataset()

        # and now update datasource
        self.update_datasource()

        # save everything so that it's quicker to reload and is available outside DLS
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'pickling results')
        pickle.dump(self.data_collection_dict,open(self.data_collection_summary_file,'wb'))

        self.emit(QtCore.SIGNAL('create_widgets_for_autoprocessing_results_only'), self.data_collection_dict)


    def parse_file_system(self):
        # only do once, ignore if just refreshing table
        if self.data_collection_dict=={}:
            if os.path.isfile(self.data_collection_summary_file):
                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'unpickling: '+self.data_collection_summary_file)
#                self.data_collection_dict = pickle.load( open(self.data_collection_summary_file, "rb" ) )
                self.data_collection_dict = cPickle.load( open(self.data_collection_summary_file, "rb" ) )

        number_of_visits_to_search=len(self.visit_list)
        search_cycle=1

        # always check for reprocessed files
        self.visit_list.append(self.initial_model_directory)

        for visit_directory in sorted(self.visit_list):
            if len(glob.glob(os.path.join(visit_directory,'processed',self.target,'*')))==0:
                continue
            progress_step=100/float(len(glob.glob(os.path.join(visit_directory,'processed',self.target,'*'))))
            progress=0

            beamline='n/a'
            if 'attic' in visit_directory:
                visit=visit_directory.split('/')[6]
                beamline=visit_directory.split('/')[3]
            else:
                if visit_directory == self.initial_model_directory:
                    visit='unknown'
                    beamline='unknown'
                else:
                    visit=visit_directory.split('/')[5]
                    beamline=visit_directory.split('/')[2]

            if visit_directory == self.initial_model_directory:
                current_directory=self.initial_model_directory
            else:
                current_directory=os.path.join(visit_directory,'processed',self.target)

            for collected_xtals in sorted(glob.glob(os.path.join(visit_directory,'processed',self.target,'*'))):
                # this step is only relevant when several samples are reviewed in one session
                if 'tmp' in collected_xtals or 'results' in collected_xtals or 'scre' in collected_xtals:
                    continue

                xtal=collected_xtals[collected_xtals.rfind('/')+1:]
                protein_name=collected_xtals.split('/')[len(collected_xtals.split('/'))-2]

                # if crystal is not in the data_collection_dict then add a new one
                if xtal not in self.data_collection_dict:
                    self.data_collection_dict[xtal]=[]

                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'Step 1 of 2: searching visit '+ \
                                                                       str(search_cycle)+' of '+str(number_of_visits_to_search)+ \
                                                                       ' ('+visit+'/'+xtal+')')

                # check if there is already an entry for the current run
                # obviously create if not and fill in basic information
                run_number_list=[]
                aimless_index_list=[]
                for runs in sorted(glob.glob(collected_xtals+'/*')):
                    run=runs[runs.rfind('/')+1:]
                    diffraction_image=''
                    timestamp=datetime.fromtimestamp(os.path.getmtime(runs)).strftime('%Y-%m-%d %H:%M:%S')
                    if os.path.isfile(os.path.join(visit_directory,protein_name,xtal,run+'0001.cbf')):
                        diffraction_image=os.path.join(visit_directory,protein_name,xtal,run+'0001.cbf')

                    ###############################################################
                    # image files
                    image_files_in_list=False
                    for entry in self.data_collection_dict[xtal]:
                        image_files_in_list=False
                        if entry[0]=='image':
                            if entry[0]=='image' and entry[1]==visit and entry[2]==run:
                                image_files_in_list=True
                                break

                    if not image_files_in_list:
                        if run_number_list==[]:
                            run_number=0                                # every run gets and keeps(!) a unique digit assigned
                        else:                                           # the order is arbitrary
                            run_number=max(run_number_list)+1
                        run_number_list.append(run_number)

                    if not image_files_in_list:
                        image_list=[]
                        # we're expecting exactly 5 images: 1 x distl plot; 4 x crystal centring images
                        # for all the ones that are not present, IMAGE_NOT_AVAILABLE.png from
                        # $XChemExplorer_DIR/image will be used instead
                        image_counter=0
                        # first four images are the crystal centring images
                        for image in sorted(glob.glob(os.path.join(visit_directory,'jpegs',self.target,xtal,'*t.png'))):
                            if run in image:
                                image_name=image[image.rfind('/')+1:]
                                image_file=open(image,"rb")
                                image_string=base64.b64encode(image_file.read())
                                image_list.append( [image_name,image_string] )
                                image_counter+=1
                        while image_counter < 4:
                            image_file=open( os.path.join(os.getenv('XChemExplorer_DIR'),'image','IMAGE_NOT_AVAILABLE.png') ,"rb")
                            image_string=base64.b64encode(image_file.read())
                            image_list.append( ['image_'+str(image_counter)+'.png',image_string] )
                            image_counter+=1
                        # now comes the distl plot
                        image_name=run+'.png'
                        if os.path.isfile(os.path.join(visit_directory,'jpegs',self.target,xtal,run+'.png')):
                            if os.stat(os.path.join(visit_directory,'jpegs',self.target,xtal,run+'.png')).st_size != 0:
                                image_file=open(os.path.join(visit_directory,'jpegs',self.target,xtal,run+'.png'),"rb")
                            else:
                                image_file=open( os.path.join(os.getenv('XChemExplorer_DIR'),'image','IMAGE_NOT_AVAILABLE.png') ,"rb")
                            image_string=base64.b64encode(image_file.read())
                            image_list.append( [image_name,image_string] )
                        else:
                            image_file=open( os.path.join(os.getenv('XChemExplorer_DIR'),'image','IMAGE_NOT_AVAILABLE.png') ,"rb")
                            image_string=base64.b64encode(image_file.read())
                            image_list.append( [image_name,image_string] )
                        # and finally comes the html page
                        if os.path.isfile(os.path.join(visit_directory,'jpegs',self.target,xtal,run+'index.html')):
                            html_summary=os.path.join(visit_directory,'jpegs',self.target,xtal,run+'index.html')
                        else:
                            html_summary=''
                        self.data_collection_dict[xtal].append(['image',visit,run,timestamp,image_list,
                                                                diffraction_image,run_number,html_summary,beamline])


                    # before we start, check if there are already entries in the aimless_index_list
                    # this list contains integers which serve a unique identifier for each autoprocessing outcome
                    for entry in self.data_collection_dict[xtal]:
                        if entry[0]=='logfile':
                            aimless_index_list.append(entry[7])
                    if aimless_index_list==[]:
                        aimless_index=0
                    else:
                        aimless_index=max(aimless_index_list)+1

                    ##########################################################################
                    # aimless & Dimple information
                    # first for xia2 runs
                    for file_name in glob.glob(os.path.join(visit_directory,'processed',protein_name,xtal,run,'xia2','*','LogFiles','*aimless.log')):
                        db_dict={   'DataCollectionVisit':              visit,
                                    'DataCollectionBeamline':           beamline,
                                    'DataCollectionDate':               timestamp,
                                    'DataProcessingPathToLogfile':      file_name,
                                    'DataProcessingPathToMTZfile':      '',
                                    'DataProcessingLOGfileName':        file_name[file_name.rfind('/')+1:],
                                    'DataProcessingDirectoryOriginal':  '/'.join(file_name.split('/')[:len(file_name.split('/'))-2])    }
                        # try to find free.mtz file
                        data_path='/'.join(file_name.split('/')[:len(file_name.split('/'))-2])
                        for data_file in glob.glob(os.path.join(data_path,'DataFiles','*')):
                            if 'free' in data_file:
                                db_dict['DataProcessingPathToMTZfile']=data_file
                                db_dict['DataProcessingMTZfileName']=data_file[data_file.rfind('/')+1:]
                                break
                        autoproc=file_name.split('/')[len(file_name.split('/'))-3]
                        found_autoproc=False
                        for n,entry in enumerate(self.data_collection_dict[xtal]):
                            if entry[0]=='logfile' and entry[1]==visit and entry[2]==run and entry[4]==autoproc:
                                found_autoproc=True
                                # need to check this because user may have run this later and otherwise he would need to delete pkl file to pick it up
                                if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')):
                                    dimple_file=os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')
                                    pdb_info=parse().PDBheader(dimple_file)
                                    db_dict_old=self.data_collection_dict[xtal][n][6]
                                    db_dict_old['DataProcessingPathToDimplePDBfile']=dimple_file
                                    db_dict_old['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                    db_dict_old['DataProcessingRcryst']  = pdb_info['Rcryst']
                                    db_dict_old['DataProcessingRfree'] = pdb_info['Rfree']
                                    self.data_collection_dict[xtal][n][6]=db_dict_old
                        if not found_autoproc:  # i.e. this run is not in pkl file yet
                            aimless_results=parse().read_aimless_logfile(file_name)
                            db_dict.update(aimless_results)
                            # first check if user ran dimple already manually
                            if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')):
                                dimple_file=os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')
                                pdb_info=parse().PDBheader(dimple_file)
                                db_dict['DataProcessingPathToDimplePDBfile']=dimple_file
                                db_dict['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                db_dict['DataProcessingRcryst']  = pdb_info['Rcryst']
                                db_dict['DataProcessingRfree'] = pdb_info['Rfree']
                            # only then start looking into the regular processed folder
                            elif os.path.isfile(os.path.join(visit_directory,'processed',protein_name,xtal,run,'xia2',autoproc,'dimple','final.pdb')):
                                dimple_file=os.path.join(visit_directory,'processed',protein_name,xtal,run,'xia2',autoproc,'dimple','final.pdb')
                                pdb_info=parse().PDBheader(dimple_file)
                                db_dict['DataProcessingPathToDimplePDBfile']=dimple_file
                                db_dict['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                db_dict['DataProcessingRcryst']  = pdb_info['Rcryst']
                                db_dict['DataProcessingRfree'] = pdb_info['Rfree']
                            else:
                                db_dict['DataProcessingPathToDimplePDBfile']=''
                                db_dict['DataProcessingPathToDimpleMTZfile']=''
                                db_dict['DataProcessingRcryst']  = '999'
                                db_dict['DataProcessingRfree'] = '999'
                            db_dict['DataProcessingProgram']=autoproc
                            # Note: [8]: best automatically selected file=True
                            #       [9]: the moment the user changes the selection manully this changes to True
                            self.data_collection_dict[xtal].append(['logfile',visit,run,timestamp,autoproc,file_name,db_dict,aimless_index,False,False])
                            aimless_index+=1




                    # then exactly the same for fast_dp
                    if os.path.isfile(os.path.join(runs,'fast_dp','aimless.log')):
                        file_name=os.path.join(runs,'fast_dp','aimless.log')
                        db_dict={   'DataCollectionVisit':              visit,
                                    'DataCollectionBeamline':           beamline,
                                    'DataCollectionDate':               timestamp,
                                    'DataProcessingPathToLogfile':      file_name,
                                    'DataProcessingPathToMTZfile':      '',
                                    'DataProcessingLOGfileName':        'aimless.log',
                                    'DataProcessingDirectoryOriginal':  os.path.join(runs,'fast_dp')    }
                        if os.path.isfile(os.path.join(runs,'fast_dp','fast_dp.mtz')):
                            db_dict['DataProcessingPathToMTZfile']=os.path.join(runs,'fast_dp','fast_dp.mtz')
                            db_dict['DataProcessingMTZfileName']='fast_dp.mtz'
                        autoproc=file_name.split('/')[len(file_name.split('/'))-2]
                        found_autoproc=False
                        for n,entry in enumerate(self.data_collection_dict[xtal]):
#                        for entry in self.data_collection_dict[xtal]:
                            if len(entry)>=9:
                                if entry[0]=='logfile' and entry[1]==visit and entry[2]==run and entry[4]==autoproc:
                                    found_autoproc=True
                                    # need to check this because user may have run this later and otherwise he would need to delete pkl file to pick it up
                                    if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')):
                                        dimple_file=os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')
                                        pdb_info=parse().PDBheader(dimple_file)
                                        db_dict_old=self.data_collection_dict[xtal][n][6]
                                        db_dict_old['DataProcessingPathToDimplePDBfile']=dimple_file
                                        db_dict_old['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                        db_dict_old['DataProcessingRcryst']  = pdb_info['Rcryst']
                                        db_dict_old['DataProcessingRfree'] = pdb_info['Rfree']
                                        self.data_collection_dict[xtal][n][6]=db_dict_old
                        if not found_autoproc:
                            aimless_results=parse().read_aimless_logfile(file_name)
                            db_dict.update(aimless_results)
                            # first check if user ran dimple already manually
                            if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')):
                                dimple_file=os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')
                                pdb_info=parse().PDBheader(dimple_file)
                                db_dict['DataProcessingPathToDimplePDBfile']=dimple_file
                                db_dict['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                db_dict['DataProcessingRcryst']  = pdb_info['Rcryst']
                                db_dict['DataProcessingRfree'] = pdb_info['Rfree']
                            elif os.path.isfile(os.path.join(runs,'fast_dp','dimple','final.pdb')):
                                dimple_file=os.path.join(runs,'fast_dp','dimple','final.pdb')
                                pdb_info=parse().PDBheader(dimple_file)
                                db_dict['DataProcessingPathToDimplePDBfile']=dimple_file
                                db_dict['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                db_dict['DataProcessingRcryst']  = pdb_info['Rcryst']
                                db_dict['DataProcessingRfree'] = pdb_info['Rfree']
                            else:
                                db_dict['DataProcessingPathToDimplePDBfile']=''
                                db_dict['DataProcessingPathToDimpleMTZfile']=''
                                db_dict['DataProcessingRcryst']  = '999'
                                db_dict['DataProcessingRfree'] = '999'
                            db_dict['DataProcessingProgram']=autoproc
                            self.data_collection_dict[xtal].append(['logfile',visit,run,timestamp,autoproc,file_name,db_dict,aimless_index,False,False])
                            aimless_index+=1

                    # then exactly the same for autoPROC
                    if os.path.isfile(os.path.join(runs,'autoPROC','ap-run','aimless.log')):
                        file_name=os.path.join(runs,'autoPROC','ap-run','aimless.log')
                        db_dict={   'DataCollectionVisit':              visit,
                                    'DataCollectionBeamline':           beamline,
                                    'DataCollectionDate':               timestamp,
                                    'DataProcessingPathToLogfile':      file_name,
                                    'DataProcessingPathToMTZfile':      '',
                                    'DataProcessingLOGfileName':        'aimless.log',
                                    'DataProcessingDirectoryOriginal':  os.path.join(runs,'autoPROC')   }
                        if os.path.isfile(os.path.join(runs,'autoPROC','ap-run','truncate-unique.mtz')):
                            db_dict['DataProcessingPathToMTZfile']=os.path.join(runs,'autoPROC','ap-run','truncate-unique.mtz')
                            db_dict['DataProcessingMTZfileName']='truncate-unique.mtz'
                        autoproc=file_name.split('/')[len(file_name.split('/'))-3]
                        found_autoproc=False
                        for n,entry in enumerate(self.data_collection_dict[xtal]):
#                        for entry in self.data_collection_dict[xtal]:
                            if len(entry)>=9:
                                if entry[0]=='logfile' and entry[1]==visit and entry[2]==run and entry[4]==autoproc:
                                    found_autoproc=True
                                    # need to check this because user may have run this later and otherwise he would need to delete pkl file to pick it up
                                    if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')):
                                        dimple_file=os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')
                                        pdb_info=parse().PDBheader(dimple_file)
                                        db_dict_old=self.data_collection_dict[xtal][n][6]
                                        db_dict_old['DataProcessingPathToDimplePDBfile']=dimple_file
                                        db_dict_old['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                        db_dict_old['DataProcessingRcryst']  = pdb_info['Rcryst']
                                        db_dict_old['DataProcessingRfree'] = pdb_info['Rfree']
                                        self.data_collection_dict[xtal][n][6]=db_dict_old
                        if not found_autoproc:
                            aimless_results=parse().read_aimless_logfile(file_name)
                            db_dict.update(aimless_results)
                            # first check if user ran dimple already manually
                            if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')):
                                dimple_file=os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')
                                pdb_info=parse().PDBheader(dimple_file)
                                db_dict['DataProcessingPathToDimplePDBfile']=dimple_file
                                db_dict['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                db_dict['DataProcessingRcryst']  = pdb_info['Rcryst']
                                db_dict['DataProcessingRfree'] = pdb_info['Rfree']
                            elif os.path.isfile(os.path.join(runs,'autoPROC','dimple','final.pdb')):
                                dimple_file=os.path.join(runs,'autoPROC','dimple','final.pdb')
                                pdb_info=parse().PDBheader(dimple_file)
                                db_dict['DataProcessingPathToDimplePDBfile']=dimple_file
                                db_dict['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                db_dict['DataProcessingRcryst']  = pdb_info['Rcryst']
                                db_dict['DataProcessingRfree'] = pdb_info['Rfree']
                            else:
                                db_dict['DataProcessingPathToDimplePDBfile']=''
                                db_dict['DataProcessingPathToDimpleMTZfile']=''
                                db_dict['DataProcessingRcryst']  = '999'
                                db_dict['DataProcessingRfree'] = '999'
                            db_dict['DataProcessingProgram']=autoproc
                            self.data_collection_dict[xtal].append(['logfile',visit,run,timestamp,autoproc,file_name,db_dict,aimless_index,False,False])
                            aimless_index+=1



                progress += progress_step
                self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

            search_cycle+=1

        # finally decide which dataset should be pre-selected
        self.select_best_dataset()

        # and now update datasource
        self.update_datasource()

        # save everything so that it's quicker to reload and is available outside DLS
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'pickling results')
#        pickle.dump(self.data_collection_dict,open(self.data_collection_summary_file,'wb'))
        cPickle.dump(self.data_collection_dict,open(self.data_collection_summary_file,'wb'))

        self.emit(QtCore.SIGNAL('create_widgets_for_autoprocessing_results_only'), self.data_collection_dict)



    def parse_file_system_NEW(self):
        # only do once, ignore if just refreshing table
        if self.data_collection_dict=={}:
            if os.path.isfile(self.data_collection_summary_file):
                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'unpickling: '+self.data_collection_summary_file)
#                self.data_collection_dict = pickle.load( open(self.data_collection_summary_file, "rb" ) )
                self.data_collection_dict = cPickle.load( open(self.data_collection_summary_file, "rb" ) )

        number_of_visits_to_search=len(self.visit_list)
        search_cycle=1

        # always check for reprocessed files
        if self.initial_model_directory not in self.visit_list:
            self.visit_list.append(self.initial_model_directory)


        for visit_directory in sorted(self.visit_list):

            if visit_directory == self.initial_model_directory:
                current_directory=self.initial_model_directory
            else:
                current_directory=os.path.join(visit_directory,'processed',self.target)

            if len(glob.glob(os.path.join(current_directory,'*')))==0:
                continue

            self.Logfile.insert('checking for new data processing results in '+current_directory)

            beamline='n/a'
            if 'attic' in visit_directory:
                visit=visit_directory.split('/')[6]
                beamline=visit_directory.split('/')[3]
            else:
                if visit_directory == self.initial_model_directory:
                    visit='unknown'
                    beamline='unknown'
                else:
                    visit=visit_directory.split('/')[5]
                    beamline=visit_directory.split('/')[2]

            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'parsing '+os.path.join('/dls_sw',beamline,'logs')+' for GDA logiles')
            nr_gda_logfiles=XChemMain.get_nr_files_from_gda_log_folder(beamline)
            if nr_gda_logfiles > 0:
                progress_step=100/float(nr_gda_logfiles)
            else:
                progress_step=1
            progress=0

            gda_pin_dict={}
            for files in glob.glob(os.path.join('/dls_sw',beamline,'logs','gda_server*')):
                self.Logfile.insert('parsing '+files+' for sampleID and pinID')
                gda_pin_dict=XChemMain.append_dict_of_gda_barcodes(gda_pin_dict,files,self.xce_logfile)
                progress += progress_step
                self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

#            gda_pin_dict=XChemMain.get_dict_of_gda_barcodes(beamline)

            progress_step=100/float(len(glob.glob(os.path.join(current_directory,'*'))))
            progress=0

            for collected_xtals in sorted(glob.glob(os.path.join(current_directory,'*'))):
                # this step is only relevant when several samples are reviewed in one session
                if 'tmp' in collected_xtals or 'results' in collected_xtals or 'scre' in collected_xtals:
                    continue

                xtal=collected_xtals[collected_xtals.rfind('/')+1:]
                if visit_directory == self.initial_model_directory:
                    protein_name=''
                else:
                    protein_name=collected_xtals.split('/')[len(collected_xtals.split('/'))-2]

                # if crystal is not in the data_collection_dict then add a new one
                found_processing=False
                if xtal not in self.data_collection_dict:
                    if visit_directory == self.initial_model_directory:
                        for new_run in glob.glob(os.path.join(collected_xtals,'processed','*')):
                            if new_run[new_run.rfind('/')+1:].startswith('run'):
                                found_processing=True
                    else:
                        found_processing=True
                    if found_processing:
                        self.data_collection_dict[xtal]=[]

                if xtal in gda_pin_dict:
                    gda_pin_id=gda_pin_dict[xtal]
                    self.data_source.execute_statement("update mainTable set DataCollectionPinBarcode ='%s' where CrystalName = '%s'"%(gda_pin_id,xtal))


                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'Step 1 of 2: searching visit '+ \
                                                                       str(search_cycle)+' of '+str(number_of_visits_to_search)+ \
                                                                       ' ('+visit+'/'+xtal+')')

                # check if there is already an entry for the current run
                # obviously create if not and fill in basic information
                run_number_list=[]
                aimless_index_list=[]

                if visit_directory == self.initial_model_directory:
                    current_run=os.path.join(collected_xtals,'processed')
                else:
                    current_run=collected_xtals

                for runs in sorted(glob.glob(os.path.join(current_run,'*'))):
                    run=runs[runs.rfind('/')+1:]
#                    print 'runs',run
                    diffraction_image=''
                    timestamp=datetime.fromtimestamp(os.path.getmtime(runs)).strftime('%Y-%m-%d %H:%M:%S')
                    if os.path.isfile(os.path.join(visit_directory,protein_name,xtal,run+'0001.cbf')):
                        diffraction_image=os.path.join(visit_directory,protein_name,xtal,run+'0001.cbf')

                    ###############################################################
                    # image files
                    image_files_in_list=False
                    for entry in self.data_collection_dict[xtal]:
                        image_files_in_list=False
                        if entry[0]=='image':
                            if entry[0]=='image' and entry[1]==visit and entry[2]==run:
                                image_files_in_list=True
                                break

                    if not image_files_in_list:
                        if run_number_list==[]:
                            run_number=0                                # every run gets and keeps(!) a unique digit assigned
                        else:                                           # the order is arbitrary
                            run_number=max(run_number_list)+1
                        run_number_list.append(run_number)

                    if not image_files_in_list:
                        image_list=[]
                        # we're expecting exactly 5 images: 1 x distl plot; 4 x crystal centring images
                        # for all the ones that are not present, IMAGE_NOT_AVAILABLE.png from
                        # $XChemExplorer_DIR/image will be used instead
                        image_counter=0
                        # first four images are the crystal centring images
                        for image in sorted(glob.glob(os.path.join(visit_directory,'jpegs',self.target,xtal,'*t.png'))):
                            if run in image:
                                image_name=image[image.rfind('/')+1:]
                                image_file=open(image,"rb")
                                image_string=base64.b64encode(image_file.read())
                                image_list.append( [image_name,image_string] )
                                image_counter+=1
                        while image_counter < 4:
                            image_file=open( os.path.join(os.getenv('XChemExplorer_DIR'),'image','IMAGE_NOT_AVAILABLE.png') ,"rb")
                            image_string=base64.b64encode(image_file.read())
                            image_list.append( ['image_'+str(image_counter)+'.png',image_string] )
                            image_counter+=1
                        # now comes the distl plot
                        image_name=run+'.png'
                        if os.path.isfile(os.path.join(visit_directory,'jpegs',self.target,xtal,run+'.png')):
                            if os.stat(os.path.join(visit_directory,'jpegs',self.target,xtal,run+'.png')).st_size != 0:
                                image_file=open(os.path.join(visit_directory,'jpegs',self.target,xtal,run+'.png'),"rb")
                            else:
                                image_file=open( os.path.join(os.getenv('XChemExplorer_DIR'),'image','IMAGE_NOT_AVAILABLE.png') ,"rb")
                            image_string=base64.b64encode(image_file.read())
                            image_list.append( [image_name,image_string] )
                        else:
                            image_file=open( os.path.join(os.getenv('XChemExplorer_DIR'),'image','IMAGE_NOT_AVAILABLE.png') ,"rb")
                            image_string=base64.b64encode(image_file.read())
                            image_list.append( [image_name,image_string] )
                        # and finally comes the html page
                        if os.path.isfile(os.path.join(visit_directory,'jpegs',self.target,xtal,run+'index.html')):
                            html_summary=os.path.join(visit_directory,'jpegs',self.target,xtal,run+'index.html')
                        else:
                            html_summary=''
                        self.data_collection_dict[xtal].append(['image',visit,run,timestamp,image_list,
                                                                diffraction_image,run_number,html_summary,beamline])


                    # before we start, check if there are already entries in the aimless_index_list
                    # this list contains integers which serve a unique identifier for each autoprocessing outcome
                    for entry in self.data_collection_dict[xtal]:
                        if entry[0]=='logfile':
                            aimless_index_list.append(entry[7])
                    if aimless_index_list==[]:
                        aimless_index=0
                    else:
                        aimless_index=max(aimless_index_list)+1

                    if visit_directory != self.initial_model_directory:

                        ##########################################################################
                        # aimless & Dimple information
                        # first for xia2 runs
                        for file_name in glob.glob(os.path.join(visit_directory,'processed',protein_name,xtal,run,'xia2','*','LogFiles','*aimless.log')):
                            db_dict={   'DataCollectionVisit':              visit,
                                        'DataCollectionBeamline':           beamline,
                                        'DataCollectionDate':               timestamp,
                                        'DataProcessingPathToLogfile':      file_name,
                                        'DataProcessingPathToMTZfile':      '',
                                        'DataProcessingLOGfileName':        file_name[file_name.rfind('/')+1:],
                                        'DataProcessingDirectoryOriginal':  '/'.join(file_name.split('/')[:len(file_name.split('/'))-2])    }
                            # try to find free.mtz file
                            data_path='/'.join(file_name.split('/')[:len(file_name.split('/'))-2])
                            for data_file in glob.glob(os.path.join(data_path,'DataFiles','*')):
                                if 'free' in data_file:
                                    db_dict['DataProcessingPathToMTZfile']=data_file
                                    db_dict['DataProcessingMTZfileName']=data_file[data_file.rfind('/')+1:]
                                    break
                            autoproc=file_name.split('/')[len(file_name.split('/'))-3]
                            found_autoproc=False
                            for n,entry in enumerate(self.data_collection_dict[xtal]):
                                if entry[0]=='logfile' and entry[1]==visit and entry[2]==run and entry[4]==autoproc:
                                    found_autoproc=True
                                    # need to check this because user may have run this later and otherwise he would need to delete pkl file to pick it up
                                    if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')):
                                        dimple_file=os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')
                                        pdb_info=parse().PDBheader(dimple_file)
                                        db_dict_old=self.data_collection_dict[xtal][n][6]
                                        db_dict_old['DataProcessingPathToDimplePDBfile']=dimple_file
                                        db_dict_old['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                        db_dict_old['DataProcessingRcryst']  = pdb_info['Rcryst']
                                        db_dict_old['DataProcessingRfree'] = pdb_info['Rfree']
                                        self.data_collection_dict[xtal][n][6]=db_dict_old
                            if not found_autoproc:  # i.e. this run is not in pkl file yet
                                aimless_results=parse().read_aimless_logfile(file_name)
                                db_dict.update(aimless_results)
                                # first check if user ran dimple already manually
                                if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')):
                                    dimple_file=os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')
                                    pdb_info=parse().PDBheader(dimple_file)
                                    db_dict['DataProcessingPathToDimplePDBfile']=dimple_file
                                    db_dict['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                    db_dict['DataProcessingRcryst']  = pdb_info['Rcryst']
                                    db_dict['DataProcessingRfree'] = pdb_info['Rfree']
                                # only then start looking into the regular processed folder
                                elif os.path.isfile(os.path.join(visit_directory,'processed',protein_name,xtal,run,'xia2',autoproc,'dimple','final.pdb')):
                                    dimple_file=os.path.join(visit_directory,'processed',protein_name,xtal,run,'xia2',autoproc,'dimple','final.pdb')
                                    pdb_info=parse().PDBheader(dimple_file)
                                    db_dict['DataProcessingPathToDimplePDBfile']=dimple_file
                                    db_dict['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                    db_dict['DataProcessingRcryst']  = pdb_info['Rcryst']
                                    db_dict['DataProcessingRfree'] = pdb_info['Rfree']
                                else:
                                    db_dict['DataProcessingPathToDimplePDBfile']=''
                                    db_dict['DataProcessingPathToDimpleMTZfile']=''
                                    db_dict['DataProcessingRcryst']  = '999'
                                    db_dict['DataProcessingRfree'] = '999'
                                db_dict['DataProcessingProgram']=autoproc
                                # Note: [8]: best automatically selected file=True
                                #       [9]: the moment the user changes the selection manully this changes to True
                                self.data_collection_dict[xtal].append(['logfile',visit,run,timestamp,autoproc,file_name,db_dict,aimless_index,False,False])
                                aimless_index+=1


                        # then exactly the same for fast_dp
                        if os.path.isfile(os.path.join(runs,'fast_dp','aimless.log')):
                            file_name=os.path.join(runs,'fast_dp','aimless.log')
                            db_dict={   'DataCollectionVisit':              visit,
                                        'DataCollectionBeamline':           beamline,
                                        'DataCollectionDate':               timestamp,
                                        'DataProcessingPathToLogfile':      file_name,
                                        'DataProcessingPathToMTZfile':      '',
                                        'DataProcessingLOGfileName':        'aimless.log',
                                        'DataProcessingDirectoryOriginal':  os.path.join(runs,'fast_dp')    }
                            if os.path.isfile(os.path.join(runs,'fast_dp','fast_dp.mtz')):
                                db_dict['DataProcessingPathToMTZfile']=os.path.join(runs,'fast_dp','fast_dp.mtz')
                                db_dict['DataProcessingMTZfileName']='fast_dp.mtz'
                            autoproc=file_name.split('/')[len(file_name.split('/'))-2]
                            found_autoproc=False
                            for n,entry in enumerate(self.data_collection_dict[xtal]):
#                            for entry in self.data_collection_dict[xtal]:
                                if len(entry)>=9:
                                    if entry[0]=='logfile' and entry[1]==visit and entry[2]==run and entry[4]==autoproc:
                                        found_autoproc=True
                                        # need to check this because user may have run this later and otherwise he would need to delete pkl file to pick it up
                                        if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')):
                                            dimple_file=os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')
                                            pdb_info=parse().PDBheader(dimple_file)
                                            db_dict_old=self.data_collection_dict[xtal][n][6]
                                            db_dict_old['DataProcessingPathToDimplePDBfile']=dimple_file
                                            db_dict_old['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                            db_dict_old['DataProcessingRcryst']  = pdb_info['Rcryst']
                                            db_dict_old['DataProcessingRfree'] = pdb_info['Rfree']
                                            self.data_collection_dict[xtal][n][6]=db_dict_old
                            if not found_autoproc:
                                aimless_results=parse().read_aimless_logfile(file_name)
                                db_dict.update(aimless_results)
                                # first check if user ran dimple already manually
                                if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')):
                                    dimple_file=os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')
                                    pdb_info=parse().PDBheader(dimple_file)
                                    db_dict['DataProcessingPathToDimplePDBfile']=dimple_file
                                    db_dict['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                    db_dict['DataProcessingRcryst']  = pdb_info['Rcryst']
                                    db_dict['DataProcessingRfree'] = pdb_info['Rfree']
                                elif os.path.isfile(os.path.join(runs,'fast_dp','dimple','final.pdb')):
                                    dimple_file=os.path.join(runs,'fast_dp','dimple','final.pdb')
                                    pdb_info=parse().PDBheader(dimple_file)
                                    db_dict['DataProcessingPathToDimplePDBfile']=dimple_file
                                    db_dict['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                    db_dict['DataProcessingRcryst']  = pdb_info['Rcryst']
                                    db_dict['DataProcessingRfree'] = pdb_info['Rfree']
                                else:
                                    db_dict['DataProcessingPathToDimplePDBfile']=''
                                    db_dict['DataProcessingPathToDimpleMTZfile']=''
                                    db_dict['DataProcessingRcryst']  = '999'
                                    db_dict['DataProcessingRfree'] = '999'
                                db_dict['DataProcessingProgram']=autoproc
                                self.data_collection_dict[xtal].append(['logfile',visit,run,timestamp,autoproc,file_name,db_dict,aimless_index,False,False])
                                aimless_index+=1

                        # then exactly the same for autoPROC
                        if os.path.isfile(os.path.join(runs,'autoPROC','ap-run','aimless.log')):
                            file_name=os.path.join(runs,'autoPROC','ap-run','aimless.log')
                            db_dict={   'DataCollectionVisit':              visit,
                                        'DataCollectionBeamline':           beamline,
                                        'DataCollectionDate':               timestamp,
                                        'DataProcessingPathToLogfile':      file_name,
                                        'DataProcessingPathToMTZfile':      '',
                                        'DataProcessingLOGfileName':        'aimless.log',
                                        'DataProcessingDirectoryOriginal':  os.path.join(runs,'autoPROC')   }
                            if os.path.isfile(os.path.join(runs,'autoPROC','ap-run','truncate-unique.mtz')):
                                db_dict['DataProcessingPathToMTZfile']=os.path.join(runs,'autoPROC','ap-run','truncate-unique.mtz')
                                db_dict['DataProcessingMTZfileName']='truncate-unique.mtz'
                            autoproc=file_name.split('/')[len(file_name.split('/'))-3]
                            found_autoproc=False
                            for n,entry in enumerate(self.data_collection_dict[xtal]):
#                            for entry in self.data_collection_dict[xtal]:
                                if len(entry)>=9:
                                    if entry[0]=='logfile' and entry[1]==visit and entry[2]==run and entry[4]==autoproc:
                                        found_autoproc=True
                                        # need to check this because user may have run this later and otherwise he would need to delete pkl file to pick it up
                                        if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')):
                                            dimple_file=os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')
                                            pdb_info=parse().PDBheader(dimple_file)
                                            db_dict_old=self.data_collection_dict[xtal][n][6]
                                            db_dict_old['DataProcessingPathToDimplePDBfile']=dimple_file
                                            db_dict_old['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                            db_dict_old['DataProcessingRcryst']  = pdb_info['Rcryst']
                                            db_dict_old['DataProcessingRfree'] = pdb_info['Rfree']
                                            self.data_collection_dict[xtal][n][6]=db_dict_old
                            if not found_autoproc:
                                aimless_results=parse().read_aimless_logfile(file_name)
                                db_dict.update(aimless_results)
                                # first check if user ran dimple already manually
                                if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')):
                                    dimple_file=os.path.join(self.initial_model_directory,xtal,'dimple',visit+'-'+run+autoproc,'dimple','final.pdb')
                                    pdb_info=parse().PDBheader(dimple_file)
                                    db_dict['DataProcessingPathToDimplePDBfile']=dimple_file
                                    db_dict['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                    db_dict['DataProcessingRcryst']  = pdb_info['Rcryst']
                                    db_dict['DataProcessingRfree'] = pdb_info['Rfree']
                                elif os.path.isfile(os.path.join(runs,'autoPROC','dimple','final.pdb')):
                                    dimple_file=os.path.join(runs,'autoPROC','dimple','final.pdb')
                                    pdb_info=parse().PDBheader(dimple_file)
                                    db_dict['DataProcessingPathToDimplePDBfile']=dimple_file
                                    db_dict['DataProcessingPathToDimpleMTZfile']=dimple_file.replace('.pdb','.mtz')
                                    db_dict['DataProcessingRcryst']  = pdb_info['Rcryst']
                                    db_dict['DataProcessingRfree'] = pdb_info['Rfree']
                                else:
                                    db_dict['DataProcessingPathToDimplePDBfile']=''
                                    db_dict['DataProcessingPathToDimpleMTZfile']=''
                                    db_dict['DataProcessingRcryst']  = '999'
                                    db_dict['DataProcessingRfree'] = '999'
                                db_dict['DataProcessingProgram']=autoproc
                                self.data_collection_dict[xtal].append(['logfile',visit,run,timestamp,autoproc,file_name,db_dict,aimless_index,False,False])
                                aimless_index+=1

                    else:
                        ##########################################################################
                        # and finally check for reprocessed files
                        # first for xia2 runs
                        for file_name in glob.glob(os.path.join(current_directory,xtal,'processed',run,'*','LogFiles','*aimless.log')):
                            db_dict={   'DataCollectionVisit':              visit,
                                        'DataCollectionBeamline':           beamline,
                                        'DataCollectionDate':               timestamp,
                                        'DataProcessingPathToLogfile':      file_name,
                                        'DataProcessingPathToMTZfile':      '',
                                        'DataProcessingLOGfileName':        file_name[file_name.rfind('/')+1:],
                                        'DataProcessingDirectoryOriginal':  '/'.join(file_name.split('/')[:len(file_name.split('/'))-2])    }
                            print db_dict
                            # try to find free.mtz file
                            data_path='/'.join(file_name.split('/')[:len(file_name.split('/'))-2])
                            for data_file in glob.glob(os.path.join(data_path,'DataFiles','*')):
                                if 'free' in data_file:
                                    db_dict['DataProcessingPathToMTZfile']=data_file
                                    db_dict['DataProcessingMTZfileName']=data_file[data_file.rfind('/')+1:]
                                    break
                            autoproc=file_name.split('/')[len(file_name.split('/'))-3]
                            found_autoproc=False
                            for n,entry in enumerate(self.data_collection_dict[xtal]):
                                if entry[0]=='logfile' and entry[1]==visit and entry[2]==run and entry[4]==autoproc:
                                    found_autoproc=True
                            if not found_autoproc:  # i.e. this run is not in pkl file yet
                                aimless_results=parse().read_aimless_logfile(file_name)
                                db_dict.update(aimless_results)
                                db_dict['DataProcessingProgram']=autoproc
                                # Note: [8]: best automatically selected file=True
                                #       [9]: the moment the user changes the selection manully this changes to True
                                self.data_collection_dict[xtal].append(['logfile',visit,run,timestamp,autoproc,file_name,db_dict,aimless_index,False,False])
                                aimless_index+=1





                progress += progress_step
                self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

            search_cycle+=1

        # finally decide which dataset should be pre-selected
        self.select_best_dataset()

        # and now update datasource
        self.update_datasource()

        # save everything so that it's quicker to reload and is available outside DLS
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'pickling results')
#        pickle.dump(self.data_collection_dict,open(self.data_collection_summary_file,'wb'))
        cPickle.dump(self.data_collection_dict,open(self.data_collection_summary_file,'wb'))

        self.emit(QtCore.SIGNAL('create_widgets_for_autoprocessing_results_only'), self.data_collection_dict)
