import os, sys, glob
from datetime import datetime
from PyQt4 import QtGui, QtCore
#from PyQt4.QtCore import QThread, SIGNAL

# last commited: 03/12/2015

import time
import pickle
import cPickle
import base64
import math
import subprocess
from datetime import datetime
import time

sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
from XChemUtils import process
from XChemUtils import parse
from XChemUtils import queue
from XChemUtils import mtztools
from XChemUtils import helpers
from XChemUtils import reference
import XChemDB


class update_datasource_from_file_system(QtCore.QThread):
    def __init__(self,initial_model_directory,datasource):
        QtCore.QThread.__init__(self)
        self.initial_model_directory=initial_model_directory
        self.datasource=datasource
        self.db=XChemDB.data_source(self.datasource)

    def run(self):
        progress_step=1
        if len(glob.glob(os.path.join(self.initial_model_directory,'*'))) != 0:
            progress_step=100/float(len(glob.glob(os.path.join(self.initial_model_directory,'*'))))
        else:
            progress_step=1
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        all_samples_in_datasource=self.db.get_all_samples_in_data_source_as_list()

        for directory in sorted(glob.glob(os.path.join(self.initial_model_directory,'*'))):
            xtal=directory[directory.rfind('/')+1:]
            if xtal not in all_samples_in_datasource:
                print '==> XCE inserting '+xtal+' into data source'
                self.db.execute_statement("insert into mainTable (CrystalName) values ('%s');" %xtal)
                all_samples_in_datasource.append(xtal)
            compoundID=str(self.db.get_value_from_field(xtal,'CompoundCode')[0])
            db_dict={}
            os.chdir(directory)
            sample_dict=self.db.get_db_dict_for_sample(xtal)

            if os.path.isfile('dimple.pdb'):
                if not os.path.isfile('refine.pdb'):
                    os.system('/bin/rm refine.pdb')         # this removes broken links that could trip the symlink
                    os.symlink('dimple.pdb', 'refine.pdb')
            if os.path.isfile('dimple.mtz'):
                db_dict['DimplePathToMTZ']=os.path.realpath(os.path.join(directory,'dimple.mtz'))
                dimple_mtz=db_dict['DimplePathToMTZ']
                dimple_path=dimple_mtz[:dimple_mtz.rfind('/')]
                if not os.path.isfile('refine.mtz'):
                    os.system('/bin/rm refine.mtz')
                    os.symlink('dimple.mtz', 'refine.mtz')
            if not os.path.isfile(xtal+'.free.mtz'):
                os.system('/bin/rm '+xtal+'.free.mtz')
                if os.path.isfile(os.path.join(dimple_path,'prepared2.mtz')):
                    os.symlink(os.path.join(dimple_path,'prepared2.mtz'),xtal+'.free.mtz')
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
                        if smiles.endwith('smiles'):
                            for line in open(smiles):
                                if len(line.split()) >= 1:
                                    db_dict['CompoundSMILES']=line.split()[0]
                                    db_dict['CompoundCode']=smiles[smiles.rfind('/')+1:smiles.rfind('.')]
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

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)


class create_png_and_cif_of_compound(QtCore.QThread):
    def __init__(self,external_software,initial_model_directory,compound_list,database_directory,data_source_file,todo,ccp4_scratch_directory):
        QtCore.QThread.__init__(self)
        self.external_software=external_software
        self.initial_model_directory=initial_model_directory
        self.compound_list=compound_list
        self.database_directory=database_directory
        self.data_source_file=data_source_file
        self.todo=todo
        self.ccp4_scratch_directory=ccp4_scratch_directory

    def run(self):
        # first remove all ACEDRG input scripts in ccp4_scratch directory
        print '==> XCE: removing all xce_acedrg scripts from',self.ccp4_scratch_directory
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

            if self.todo=='ALL':
                # remove symbolic links if present
                if os.path.isfile(os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.pdb')):
                    os.system('/bin/rm '+os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.pdb'))
                if os.path.isfile(os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.cif')):
                    os.system('/bin/rm '+os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.cif'))
                if os.path.isfile(os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.png')):
                    os.system('/bin/rm '+os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.png'))
                if os.path.isdir(os.path.join(self.initial_model_directory,sampleID,'compound')):
                    os.system('/bin/rm -fr '+os.path.join(self.initial_model_directory,sampleID,'compound'))

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
                                    self.external_software['qsub'],
                                    self.database_directory,
                                    self.data_source_file,
                                    self.ccp4_scratch_directory,
                                    counter )
                counter += 1

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        # submit array job at Diamond
        print '==> XCE: created input scripts for '+str(counter)+' ACEDRG jobs in '+self.ccp4_scratch_directory
        os.chdir(self.ccp4_scratch_directory)
        if counter > 1:
            Cmds = (
                    '#PBS -joe -N xce_acedrg_master\n'
                    './xce_acedrg_$SGE_TASK_ID.sh\n'
                    )
            f = open('acedrg_master.sh','w')
            f.write(Cmds)
            f.close()
            print '==> XCE: submitting array job with maximal 100 jobs running on cluster'
            print '==> XCE: using the following command:'
            print '         qsub -t 1:%s -tc 100 acedrg_master.sh' %(str(counter))
            os.system('qsub -t 1:%s -tc 100 acedrg_master.sh' %(str(counter)))

        self.emit(QtCore.SIGNAL("finished()"))


class run_dimple_on_all_autoprocessing_files(QtCore.QThread):
    def __init__(self,sample_list,initial_model_directory,external_software,ccp4_scratch_directory,database_directory,data_source_file):
        QtCore.QThread.__init__(self)
        self.sample_list=sample_list
        self.initial_model_directory=initial_model_directory
        self.queueing_system_available=external_software['qsub']
        self.ccp4_scratch_directory=ccp4_scratch_directory
        self.database_directory=database_directory
        self.data_source_file=data_source_file
    def run(self):
        progress_step=1
        if len(self.sample_list) != 0:
            progress_step=100/float(len(self.sample_list))
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        os.chdir(self.ccp4_scratch_directory)
        os.system('/bin/rm -f xce_dimple*sh')

        for n,item in enumerate(self.sample_list):

            xtal =                  item[0]
            visit_run_autoproc =    item[1]
            mtzin =                 item[2]
            ref_pdb =               item[3]
            ref_mtz =               item[4]
            ref_cif =               item[5]

            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'creating input script for '+xtal+' in '+visit_run_autoproc)

            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal)):
                os.mkdir(os.path.join(self.initial_model_directory,xtal))
            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'dimple')):
                os.mkdir(os.path.join(self.initial_model_directory,xtal,'dimple'))
            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc)):
                os.mkdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc))
            os.chdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc))
            os.system('touch dimple_run_in_progress')

            if self.queueing_system_available:
                top_line='#PBS -joe -N XCE_dimple'
            else:
                top_line='#!'+os.getenv('SHELL')

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
                    +ccp4_scratch+
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
                    )

            os.chdir(self.ccp4_scratch_directory)
            f = open('xce_dimple_%s.sh' %str(n+1),'w')
            f.write(Cmds)
            f.close()
            os.system('chmod +x xce_dimple_%s.sh' %str(n+1))
#            if self.queueing_system_available:
#                os.system('qsub xce_dimple.sh')
#            else:
#                os.system('chmod +x xce_dimple.sh')
#                os.system('./xce_dimple.sh')
            if not self.queueing_system_available:
                os.system('./xce_dimple_%s.sh' %str(n+1))

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        # submit array job at Diamond
        print '==> XCE: created input scripts for '+str(n+1)+' in '+self.ccp4_scratch_directory
        os.chdir(self.ccp4_scratch_directory)
        Cmds = (
                '#PBS -joe -N xce_dimple_master\n'
                './xce_dimple_$SGE_TASK_ID.sh\n'
                )
        f = open('dimple_master.sh','w')
        f.write(Cmds)
        f.close()
        print '==> XCE: submitting array job with maximal 100 jobs running on cluster'
        print '==> XCE: using the following command:'
        print '         qsub -t 1:%s -tc 100 dimple_master.sh' %(str(n+1))
        os.system('qsub -t 1:%s -tc 100 dimple_master.sh' %(str(n+1)))



class start_COOT(QtCore.QThread):

    def __init__(self,settings):
        QtCore.QThread.__init__(self)
        self.settings=settings

    def run(self):
        cwd=os.getcwd()
#        pickle.dump(self.settings,open(os.path.join(cwd,'XChemExplorer_settings.pkl'),'wb'))
#        os.system('cd %s\ncoot --no-guano --no-state-script --script %s' %(cwd,os.getenv('XChemExplorer_DIR')+'/lib/XChemCoot.py'))
        # coot at Diamond always or sometimes at least open in home directory, so then it won't find the .pkl file
        pickle.dump(self.settings,open(os.path.join(os.getenv('HOME'),'.xce_settings.pkl'),'wb'))
        os.system('cd %s\ncoot --no-guano --no-state-script --script %s' %(os.getenv('HOME'),os.path.join(os.getenv('XChemExplorer_DIR'),'lib','XChemCoot.py')))

class start_pandda_inspect(QtCore.QThread):

    def __init__(self,settings):
        QtCore.QThread.__init__(self)
#        self.settings=settings
        self.panddas_directory=settings['panddas_directory']

    def run(self):
        if os.getenv('SHELL') == '/bin/tcsh' or os.getenv('SHELL') == '/bin/csh':
            source_file=os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-csh')
        elif os.getenv('SHELL') == '/bin/bash':
            source_file=os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-sh')
        else:
            source_file=''

        Cmds = (
                '#!'+os.getenv('SHELL')+'\n'
                '\n'
                'source '+source_file+'\n'
                '\n'
                'cd '+self.panddas_directory+'\n'
                '\n'
                'pandda.inspect\n'
            )
        print Cmds
        os.system(Cmds)




class start_dials_image_viewer(QtCore.QThread):

    def __init__(self,diffraction_image):
        QtCore.QThread.__init__(self)
        self.diffraction_image=diffraction_image

    def run(self):
        os.system('dials.image_viewer '+self.diffraction_image)




class LATEST_save_autoprocessing_results_to_disc(QtCore.QThread):
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
                    path_to_dimple_pdbfile=db_dict['DataProcessingPathToDimplePDBfile']
                    path_to_dimple_mtzfile=db_dict['DataProcessingPathToDimpleMTZfile']

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
                    if entry[7]==selected_processing_result:
                        db_dict=entry[6]
                        db_dict['DataCollectionOutcome']=self.dataset_outcome_dict[sample]
                        db_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
                        db_dict['RefinementMTZfree'],db_dict['DimpleRcryst'],db_dict['DimpleRfree'],db_dict['RefinementOutcome'],db_dict['RefinementSpaceGroup'] =self.link_mtz_log_files_to_sample_directory(sample,autoproc,run,visit,path_to_procdir,path_to_logfile,path_to_mtzfile,mtz_filename,log_filename,dimple_destination)
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
        # move up to sample directory and link respective files
        # first remove any old symbolic links
        os.chdir(os.path.join(self.initial_model_directory,sample))
        if os.path.islink(os.path.join(self.initial_model_directory,sample,sample+'.mtz')):
            os.system('/bin/rm '+os.path.join(self.initial_model_directory,sample,sample+'.mtz'))
        if os.path.islink(os.path.join(self.initial_model_directory,sample,sample+'.log')):
            os.system('/bin/rm '+os.path.join(self.initial_model_directory,sample,sample+'.log'))
        # then link new files
        os.symlink(os.path.join(path_to_mtzfile,sample+'.mtz'),sample+'.mtz')
        os.symlink(os.path.join(path_to_logfile,sample+'.log'),sample+'.log')
        if dimple_destination != '':
            # check if dimple files exist
            if os.path.isfile(os.path.join(dimple_destination,'final.pdb')):
                # remove old symbolic links if necessary
                if os.path.isfile('dimple.pdb'): os.system('/bin/rm dimple.pdb')
                os.symlink(os.path.join(dimple_destination,'final.pdb'),'dimple.pdb')
                pdb_info=parse().PDBheader(os.path.join(dimple_destination,'final.pdb'))
                Rcryst=pdb_info['Rcryst']
                Rfree=pdb_info['Rfree']
                spg=pdb_info['SpaceGroup']
                refinement_stage='1 - Analysis Pending'

            if os.path.isfile(os.path.join(dimple_destination,'final.mtz')):
                # remove old symbolic links if necessary
                if os.path.isfile('dimple.mtz'): os.system('/bin/rm dimple.mtz')
                os.symlink(os.path.join(dimple_destination,'final.mtz'),'dimple.mtz')
            # if no refinement was carried out yet, then we also want to link the dimple files to refine.pdb/refine.log
            # so that we can look at them with the COOT plugin
            found_previous_refinement=False
            for dirs in glob.glob('*'):
                if os.path.isdir(dirs) and dirs.startswith('Refine_'):
                    found_previous_refinement=True
                    break
            if not found_previous_refinement:
                # first delete possible old symbolic links
                if os.path.isfile('refine.pdb'): os.system('/bin/rm refine.pdb')
                os.symlink('dimple.pdb','refine.pdb')
                if os.path.isfile('refine.mtz'): os.system('/bin/rm refine.mtz')
                os.symlink('dimple.mtz','refine.mtz')
            # remove any previous <sample>.free.mtz file, and link new dimple.mtz
            # so if we continue refining, then we do so against the correct file
            # think that REFMAC does not tinker with F,SIGF as long as there is no twinning
            if os.path.isfile(sample+'.free.mtz'): os.system('/bin/rm '+sample+'.free.mtz')
            os.symlink(os.path.join(dimple_destination,'final.mtz'),sample+'.free.mtz')
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





class NEW_read_autoprocessing_results_from_disc(QtCore.QThread):
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
                 data_source_file):
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

        # - open data source if possible
        # - get sampleID, xtbm
        # - search lab36 folder for respective xtal image
        # - convert to string and use in data dict
        # - but images can only be found of XCE is started in the respective labchem directory

    def run(self):
        if self.target=='=== SELECT TARGET ===':
            print '==> XCE: please select a target first'
            return

        if self.rescore_only:
            self.rescore_and_reset_pkl_file()
        else:
            self.parse_file_system()

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
            found=False
            if entry[0]=='logfile':
                index=self.data_collection_dict[xtal][n][7]
                if isinstance(entry[6],dict):
                    try:
                        if isinstance(float(entry[6]['DataProcessingUnitCellVolume']),float):
                            for reference_file in self.reference_file_list:
                                if not reference_file[4]==0:
                                    unitcell_difference=round((math.fabs(reference_file[4]-float(entry[6]['DataProcessingUnitCellVolume']))/reference_file[4])*100,1)
                                    if unitcell_difference < 5 and reference_file[3]==entry[6]['DataProcessingLattice']:
                                        select_stage_one_list.append(index)
                                        found=True
                    except ValueError:
                        pass
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
        select_stage_three_list=[]
        for index in select_stage_two_list:
            for entry in self.data_collection_dict[xtal]:
                if entry[0]=='logfile':
                    if isinstance(entry[6],dict) and entry[7]==index:
                        try:
                            ranking=float(entry[6]['DataProcessingUniqueReflectionsOverall'])*\
                                    float(entry[6]['DataProcessingCompletenessOverall'])*\
                                    float(entry[6]['DataProcessingIsigOverall'])
                            select_stage_three_list.append([index,ranking])
                        except ValueError:
                            pass
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

        for sample in sorted(self.data_collection_dict):
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'updating data source for '+sample)
            logfile_found=False
            user_changed_selection=False
            db_dict={}
            db_dict['DataCollectionOutcome']='Failed - unknown'
            tmpList=[]
            for entry in self.data_collection_dict[sample]:
                if entry[0]=='user_changed_selection':
                    user_changed_selection=True
                if entry[0]=='logfile':
                    if entry[8]:        # the best auto-selected or user selected output
                        db_dict=entry[6]
                        logfile_found=True
                        try:
                            if float(db_dict['DataProcessingResolutionHigh']) <= float(self.acceptable_low_resolution_limit_for_data):
                                db_dict['DataCollectionOutcome']='success'
                            else:
                                db_dict['DataCollectionOutcome']='Failed - low resolution'
                        except ValueError:
                            db_dict['DataCollectionOutcome']='Failed - unknown'
                        db_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
                if entry[0]=='image':
                    if len(entry) >= 9:     # need this because some older pkl files won't have the beamline added
                        tmpList.append([datetime.strptime(entry[3], '%Y-%m-%d %H:%M:%S'),entry[1],entry[2],entry[8]])
                    else:
                        tmpList.append([datetime.strptime(entry[3], '%Y-%m-%d %H:%M:%S'),entry[1],entry[2],'I04-1'])

            if not logfile_found:
                # find latest data collection in tmpList
                latest_run=max(tmpList,key=lambda x: x[0])
                db_dict={   'DataCollectionVisit':              latest_run[1],
                            'DataCollectionBeamline':           latest_run[3],
                            'DataCollectionDate':               latest_run[0],
                            'DataCollectionOutcome':            'Failed - unknown',
                            'RefinementOutcome':                '-1 - Data Collection Failed'}

            if self.rescore_only:
                self.data_source.update_insert_data_source(sample,db_dict)
            elif user_changed_selection==False:     # if user changed the selection, then ignore
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
                visit=visit_directory.split('/')[5]
                beamline=visit_directory.split('/')[2]


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
                        aimless_index=max(aimless_index_list)

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
