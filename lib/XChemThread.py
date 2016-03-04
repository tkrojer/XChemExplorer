import os, sys, glob
from datetime import datetime
from PyQt4 import QtGui, QtCore
#from PyQt4.QtCore import QThread, SIGNAL

# last commited: 03/12/2015

import time
import pickle
import base64
import math
import subprocess

sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
from XChemUtils import process
from XChemUtils import parse
from XChemUtils import queue
from XChemUtils import mtztools
from XChemUtils import helpers
from XChemUtils import reference
import XChemDB


class crystal_from(QtCore.QThread):
    def __init__(self,initial_model_directory,
                      reference_file_list,
                      data_source,
                      filename_root,
                 xtalform_dict,
                 mode):
        QtCore.QThread.__init__(self)

        self.initial_model_directory=initial_model_directory
        self.reference_file_list=reference_file_list
        self.data_source=data_source
        self.allowed_unitcell_difference_percent=10
        self.filename_root=filename_root
        self.xtalform_dict=xtalform_dict
        self.mode=mode

    def suggest_new_crystalfrom_name(self,xtalform_dict,pointgroup):
        temp = []
        found = 0
        for item in xtalform_dict:
            if str(item)[:str(item).rfind('_')]==pointgroup:
                temp.append(int(item[item.rfind('_')+1:]))
                found = 1
                Serial=item
        if found:
            Serial = max(temp) + 1
        else:
            Serial=1
        return pointgroup+'_'+str(Serial)


    def run(self):

        progress_step=100/float(len(glob.glob(self.initial_model_directory+'/*')))
        progress=0

        # xtalform_dict['name']: [pointgroup,unitcell_volume,[unitcell],spacegroup]
        # xf naming convetion: <pointgroup>_<serial>; e.g. 222_1

        db=XChemDB.data_source(self.data_source)
        db.create_missing_columns()             # always do this because additional columns might be added

        if self.mode=='suggest':
            # 1. check for XF in reference files
            for reference in self.reference_file_list:
                add_crystalform=True
                for xf in self.xtalform_dict:
                    if self.xtalform_dict[xf][0]==reference[5]:      # same pointgroup as reference
                        unitcell_difference=round((math.fabs(reference[4]-float(self.xtalform_dict[xf][1]))/reference[4])*100,1)
                        if unitcell_difference < self.allowed_unitcell_difference_percent:      # suggest new crystal form
                            add_crystalform=False
                            break
                if add_crystalform:
                    name=self.suggest_new_crystalfrom_name(self.xtalform_dict,reference[5])
                    self.xtalform_dict[name]=[reference[5],reference[4],str(reference[2]).split(),reference[1]]

            # 2. check for existing XF in datasource
            columns = ( 'CrystalFormName,'
                        'CrystalFormPointGroup,'
                        'CrystalFormVolume,'
                        'CrystalFormA,'
                        'CrystalFormB,'
                        'CrystalFormC,'
                        'CrystalFormAlpha,'
                        'CrystalFormBeta,'
                        'CrystalFormGamma,'
                        'CrystalFormSpaceGroup' )
            all_xtalforms=db.execute_statement("SELECT "+columns+" FROM mainTable;")
            for item in all_xtalforms:
                name=item[0]
                pg=item[1]
                vol=item[2]
                a=item[3]
                b=item[4]
                c=item[5]
                alpha=item[6]
                beta=item[7]
                gamma=item[8]
                spg=item[9]
                add_crystalform=True
                for xf in self.xtalform_dict:
                    if self.xtalform_dict[xf][0]==pg:      # same pointgroup as reference
                        try:
                            unitcell_difference=round((math.fabs(float(vol)-float(self.xtalform_dict[xf][1]))/float(vol))*100,1)
                        except (ValueError,TypeError):
                            unitcell_difference=100
                        if unitcell_difference < self.allowed_unitcell_difference_percent:      # suggest new crystal form
                            add_crystalform=False
                            break
                if add_crystalform:
                    if str(name)=='None':
                        pass
                    elif str(name)=='':
                        pass
                    else:
                        self.xtalform_dict[name]=[pg,vol,[a,b,c,alpha,beta,gamma],spg]
#                    if not str(name)=='None' or not name=='':
#                        self.xtalform_dict[name]=[pg,vol,[a,b,c,alpha,beta,gamma],spg]


        # 3. parse inital_model directory
        # now assiignment will be random, because first sample read will be first xtalfrom
        # in case nothing else was provided
        for sample_dir in sorted(glob.glob(self.initial_model_directory+'/*')):

            sample=sample_dir[sample_dir.rfind('/')+1:]
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'parsing initial_models folder -> '+sample)

            mtzin=''
            sample_mtz=self.filename_root.replace('${samplename}',sample)+'.mtz'
            if os.path.isfile(os.path.join(self.initial_model_directory,sample,sample_mtz)):
                mtzin=os.path.realpath(os.path.join(self.initial_model_directory,sample,sample_mtz))
            elif os.path.isfile(os.path.join(self.initial_model_directory,sample,sample+'.mtz')):
                mtzin=os.path.realpath(os.path.join(self.initial_model_directory,sample,sample+'.mtz'))
            elif os.path.isfile(os.path.join(self.initial_model_directory,sample,sample+'.free.mtz')):
                mtzin=os.path.realpath(os.path.join(self.initial_model_directory,sample,sample+'.free.mtz'))
            elif os.path.isfile(os.path.join(self.initial_model_directory,sample,'refine.mtz')):
                mtzin=os.path.realpath(os.path.join(self.initial_model_directory,sample,'refine.mtz'))
            if mtzin!='':
                mtz_pg=mtztools(mtzin).get_pointgroup_from_mtz()
                mtz_vol=mtztools(mtzin).calc_unitcell_volume_from_mtz()
                add_crystalform=True
                for xf in self.xtalform_dict:
                    if self.xtalform_dict[xf][0]==mtz_pg:      # same pointgroup as reference
                        unitcell_difference=round((math.fabs(mtz_vol-float(self.xtalform_dict[xf][1]))/mtz_vol)*100,1)
                        if unitcell_difference < self.allowed_unitcell_difference_percent:      # suggest new crystal form
                            add_crystalform=False
                            if self.mode=='assign':
                                db_dict = {
                                    'CrystalFormName':          xf,
                                    'CrystalFormSpaceGroup':    self.xtalform_dict[xf][3],
                                    'CrystalFormPointGroup':    self.xtalform_dict[xf][0],
                                    'CrystalFormA':             self.xtalform_dict[xf][2][0],
                                    'CrystalFormB':             self.xtalform_dict[xf][2][1],
                                    'CrystalFormC':             self.xtalform_dict[xf][2][2],
                                    'CrystalFormAlpha':         self.xtalform_dict[xf][2][3],
                                    'CrystalFormBeta':          self.xtalform_dict[xf][2][4],
                                    'CrystalFormGamma':         self.xtalform_dict[xf][2][5],
                                    'CrystalFormVolume':        self.xtalform_dict[xf][1]       }
                                db.update_data_source(sample,db_dict)
                            break
                if self.mode=='suggest':
                    if add_crystalform:
                        name=self.suggest_new_crystalfrom_name(self.xtalform_dict,mtz_pg)
                        unitcell=mtztools(mtzin).get_unit_cell_from_mtz()
                        spacegroup=mtztools(mtzin).get_spg_from_mtz()
                        self.xtalform_dict[name]=[mtz_pg,mtz_vol,unitcell,spacegroup]

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        if self.mode=='suggest':
            self.emit(QtCore.SIGNAL('update_xtalfrom_table'), self.xtalform_dict)



class create_png_and_cif_of_compound(QtCore.QThread):
    def __init__(self,external_software,initial_model_directory,compound_list):
        QtCore.QThread.__init__(self)
        self.external_software=external_software
        self.initial_model_directory=initial_model_directory
        self.compound_list=compound_list

    def run(self):
        progress_step=100/float(len(self.compound_list))
        progress=0
        for item in self.compound_list:
            sampleID=item[0]
            compoundID=item[1]
            smiles=item[2]
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'creating cif/png -> '+sampleID)
            if compoundID=='' or compoundID==None:
                compoundID='compound'
            helpers().make_png(self.initial_model_directory,sampleID,compoundID,smiles,self.external_software['qsub'])
            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
        self.emit(QtCore.SIGNAL("finished()"))



class run_dimple_on_selected_samples(QtCore.QThread):
    def __init__(self,settings,initial_model_dimple_dict,external_software,ccp4_scratch,filename_root):
        QtCore.QThread.__init__(self)
        self.initial_model_directory=settings['initial_model_directory']
        self.reference_directory=settings['reference_directory']
        self.initial_model_dimple_dict=initial_model_dimple_dict
        self.queueing_system_available=external_software['qsub']
        self.ccp4_scratch_directory=ccp4_scratch
        self.filename_root=filename_root

    def run(self):
        todo=0
        if len(self.initial_model_dimple_dict) != 0:
#            progress_step=100/float(len(self.initial_model_dimple_dict))
            for sample in sorted(self.initial_model_dimple_dict):
                if self.initial_model_dimple_dict[sample][0].isChecked(): todo+=1
            if todo==0:
                todo=1
            progress_step=100/float(todo)
        progress=0

        for sample in sorted(self.initial_model_dimple_dict):
            if self.initial_model_dimple_dict[sample][0].isChecked():
                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'running dimple -> '+sample)
                dimple_commands={   'project_directory': self.initial_model_directory,
                                    'delete_old': self.initial_model_dimple_dict[sample][0].isChecked(),
                                    'xtalID': sample,
                                    'compoundID': '',
                                    'smiles': '',
                                    'reference': self.reference_directory+'/'+
                                                 str(self.initial_model_dimple_dict[sample][1].currentText()),
                                    'queueing_system_available': self.queueing_system_available,
                                    'ccp4_scratch': self.ccp4_scratch_directory,
                                    'fileroot_in':  self.filename_root.replace('${samplename}',sample)  }
                process(dimple_commands).dimple()
            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
        self.emit(QtCore.SIGNAL("finished()"))


class start_COOT(QtCore.QThread):

    def __init__(self,settings):
        QtCore.QThread.__init__(self)
        self.settings=settings

    def run(self):
        pickle.dump(self.settings,open('XChemExplorer_settings.pkl','wb'))
        os.system('coot --no-guano --no-state-script --script %s' %(os.getenv('XChemExplorer_DIR')+'/lib/XChemCoot.py'))

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



class read_intial_refinement_results(QtCore.QThread):

    def __init__(self,initial_model_directory,
                      reference_file_list,
                      data_source,
                      allowed_unitcell_difference_percent,
                      filename_root,
                      update_datasource_only):
        QtCore.QThread.__init__(self)
        self.initial_model_directory=initial_model_directory
        self.reference_file_list=reference_file_list
        self.data_source=data_source
        self.allowed_unitcell_difference_percent=allowed_unitcell_difference_percent
        self.filename_root=filename_root
        self.update_datasource_only=update_datasource_only


    def run(self):

        progress_step=100/float(len(glob.glob(self.initial_model_directory+'/*')))
        progress=0

        db=XChemDB.data_source(self.data_source)
        db.create_missing_columns()

        initial_model_list=[]

        for sample_dir in sorted(glob.glob(self.initial_model_directory+'/*')):

            sample=sample_dir[sample_dir.rfind('/')+1:]
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'parsing initial_models folder -> '+sample)

            run_dimple=False
            resolution_high=''
            Rcryst='pending'
            Rfree='pending'
            spg_autoproc=''
            unitcell_autoproc=''
            spg_reference=''
            unitcell_reference=''
            unitcell_difference=''
            reference_file=''
            alert='#E0E0E0'
            outcome='Analysis Pending'
            rmsdBonds='n/a'
            rmsdAngles='n/a'

            # if XCE is used throughout the process then there will be an inital <sample>.mtz file
            # but if not then there could be a few files which serve the same purpose
            sample_mtz=self.filename_root.replace('${samplename}',sample)+'.mtz'

            spg_reference,unitcell_reference,reference_file,found_suitable_reference,\
                resolution_high,spg_autoproc,unitcell_autoproc,unitcell_difference= \
                reference(os.path.join(self.initial_model_directory,sample,sample_mtz),
                          self.reference_file_list).find_suitable_reference(self.allowed_unitcell_difference_percent)

            if os.path.isfile(os.path.join(self.initial_model_directory,sample,'refine.pdb')):
                pdb=parse().PDBheader(os.path.join(self.initial_model_directory,sample,'refine.pdb'))
                Rcryst=pdb['Rcryst']
                Rfree=pdb['Rfree']
                alert=pdb['Alert']
                rmsdBonds=pdb['rmsdBonds']
                rmsdAngles=pdb['rmsdAngles']

            elif os.path.isfile(os.path.join(self.initial_model_directory,sample,'/dimple_run_in_progress')):
                Rcryst='in progress'
                Rfree='in progress'
                alert=(51,153,255)

            if Rcryst=='pending':
                run_dimple=True

            if found_suitable_reference==True:
                run_dimple=True
            else:
                run_dimple=False

            pdb_latest=''
            path_to_refinement_folder=''
            if os.path.isfile(os.path.join(self.initial_model_directory,sample,'refine.pdb')):
                pdb_latest=os.path.realpath(os.path.join(self.initial_model_directory,sample,'refine.pdb'))
                path_to_refinement_folder=os.path.realpath(os.path.join(self.initial_model_directory,sample))

            mtz_latest=''
            if os.path.isfile(os.path.join(self.initial_model_directory,sample,'refine.mtz')):
                mtz_latest=os.path.realpath(os.path.join(self.initial_model_directory,sample,'refine.mtz'))

            mtz_free=''
            if os.path.isfile(os.path.join(self.initial_model_directory,sample,sample+'.free.mtz')):
                mtz_free=os.path.realpath(os.path.join(self.initial_model_directory,sample,sample+'.free.mtz'))

            cif_file=''
            compoundID=db.execute_statement("SELECT CompoundCode FROM mainTable WHERE CrystalName='"+sample+"';")
            if len(compoundID) >= 1:
                if os.path.isfile(os.path.join(self.initial_model_directory,sample,str(compoundID[0][0])+'.cif')):
                    cif_file=os.path.realpath(os.path.join(self.initial_model_directory,sample,str(compoundID[0][0])+'.cif'))

            if os.path.isfile(os.path.join(self.initial_model_directory,sample,sample+'.log')):
                scaling_logfile=os.path.realpath(os.path.join(self.initial_model_directory,sample,sample+'.log'))
                db_dict={'DataProcessingPathToLogfile': scaling_logfile }
                # update path to logfile: might have changed if data were copied to different location
                db.update_data_source(sample,db_dict)
                db_dict=parse().read_aimless_logfile(scaling_logfile)
                # actual values will only be updated if not present;
                db.update_insert_not_null_fields_only(sample,db_dict)

            # update Dimple stuff (useful for PANDDAs)
            db_dict={}
            if os.path.isfile(os.path.join(self.initial_model_directory,sample,'Dimple','dimple','final.pdb')):
                db_dict['DimplePathToPDB']=os.path.join(self.initial_model_directory,sample,'Dimple','dimple','final.pdb')
                pdb=parse().PDBheader(os.path.join(self.initial_model_directory,sample,'refine.pdb'))
                db_dict['DimpleRcryst']=pdb['Rcryst']
                db_dict['DimpleRfree']=pdb['Rfree']
                db_dict['DimpleResolutionHigh']=pdb['ResolutionHigh']
            if os.path.isfile(os.path.join(self.initial_model_directory,sample,'Dimple','dimple','final.mtz')):
                db_dict['DimplePathToMTZ']=os.path.join(self.initial_model_directory,sample,'Dimple','dimple','final.mtz')
            if db_dict != {}:
                print db_dict
                db.update_insert_not_null_fields_only(sample,db_dict)

            # update data source only if field is null
            db_dict={   'RefinementOutcome':            outcome,
                        'RefinementMTZfree':            mtz_free    }
            db.update_insert_not_null_fields_only(sample,db_dict)

            # update to reflect current state
            db_dict={   'RefinementRcryst':                 Rcryst,
                        'RefinementRfree':                  Rfree,
                        'RefinementRmsdBonds':              rmsdBonds,
                        'RefinementRmsdAngles':             rmsdAngles,
                        'RefinementPDB_latest':             pdb_latest,
                        'RefinementMTZ_latest':             mtz_latest,
                        'RefinementPathToRefinementFolder': path_to_refinement_folder,
                        'RefinementCIF':                    cif_file}
            db.update_data_source(sample,db_dict)

            initial_model_list.append( [ sample,
                                        run_dimple,
                                        resolution_high,
                                        Rcryst,
                                        Rfree,
                                        spg_autoproc,
                                        spg_reference,
                                        unitcell_difference,
                                        reference_file,
                                        unitcell_autoproc,
                                        unitcell_reference,
                                        alert ] )

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
        if not self.update_datasource_only:
            self.emit(QtCore.SIGNAL('create_initial_model_table'), initial_model_list)



class save_autoprocessing_results_to_disc(QtCore.QThread):
    def __init__(self,dataset_outcome_dict,data_collection_table_dict,data_collection_statistics_dict,
                 database_directory,data_source_file,initial_model_directory):
        QtCore.QThread.__init__(self)
        self.dataset_outcome_dict=dataset_outcome_dict
        self.data_collection_table_dict=data_collection_table_dict
        self.data_collection_statistics_dict=data_collection_statistics_dict
        self.database_directory=database_directory
        self.data_source_file=data_source_file
        self.initial_model_directory=initial_model_directory

    def run(self):
        if not len(self.dataset_outcome_dict)==0:
            progress_step=100/float(len(self.dataset_outcome_dict))
        progress=0
        data_source=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
        for sample in sorted(self.dataset_outcome_dict):
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'writing files from data processing to inital_model folder -> '+sample)
            outcome=''
            for button in self.dataset_outcome_dict[sample]:
                if button.isChecked():
                    outcome=button.text()
            indexes=self.data_collection_table_dict[sample].selectionModel().selectedRows()
            if indexes == []:       # i.e. no logfile exists
                logfile=None
            else:
                for index in sorted(indexes):
                    logfile=self.data_collection_statistics_dict[sample][index.row()][1]
            if self.data_source_file != '':
                data_dict=data_source.get_data_dict_to_save_autoprocessing_results_to_data_source(sample,str(outcome),logfile)
                data_source.update_data_source(sample,data_dict)
#                data_source.save_autoprocessing_results_to_data_source(sample,str(outcome),logfile)

            # create all the directories if necessary
            if not os.path.isdir(os.path.join(self.initial_model_directory,sample)):
                os.mkdir(os.path.join(self.initial_model_directory,sample))
            if not os.path.isdir(os.path.join(self.initial_model_directory,sample,'autoprocessing')):
                os.mkdir(os.path.join(self.initial_model_directory,sample,'autoprocessing'))

            if logfile != None:
                path_to_logfile=self.data_collection_statistics_dict[sample][index.row()][1]
                # copy files
                if 'xia2' in path_to_logfile:
                    path_to_procdir=os.path.join('/',*path_to_logfile.split('/')[:len(path_to_logfile.split('/'))-2])
                if 'fast_dp' in path_to_logfile:
                    path_to_procdir=os.path.join('/',*path_to_logfile.split('/')[:len(path_to_logfile.split('/'))-1])
                os.system('/bin/cp -Rf '+path_to_procdir+' '+os.path.join(self.initial_model_directory,sample,'autoprocessing'))

                # link files
                if 'xia2' in path_to_logfile:
                    os.chdir(os.path.join(self.initial_model_directory,sample))
                    for datafile in glob.glob('autoprocessing/*/DataFiles/*'):
                        if datafile.endswith('free.mtz'):
                            if os.path.isfile(sample+'.mtz'):
                                os.system('/bin/rm '+sample+'.mtz')
                            os.symlink(datafile,sample+'.mtz')
                            break
                    for logfile in glob.glob('autoprocessing/*/LogFiles/*'):
                        if logfile.endswith('aimless.log'):
                            if os.path.isfile(sample+'.log'):
                                os.system('/bin/rm '+sample+'.log')
                            os.symlink(logfile,sample+'.log')
                            break
                if 'fast_dp' in path_to_logfile:
                    os.chdir(os.path.join(self.initial_model_directory,sample,'autoprocessing','fast_dp'))
                    os.system("ctruncate -hklin fast_dp.mtz "
                              "-hklout ctruncate.mtz -colin '/*/*/[IMEAN,SIGIMEAN]' "
                              "> ctruncate.log")
                    os.chdir(os.path.join(self.initial_model_directory,sample))
                    if os.path.isfile(sample+'.mtz'):
                        os.system('/bin/rm '+sample+'.mtz')
                    if os.path.isfile(sample+'.log'):
                        os.system('/bin/rm '+sample+'.log')
                    os.symlink('autoprocessing/fast_dp/aimless.log',sample+'.log')
                    os.symlink('autoprocessing/fast_dp/ctruncate.mtz',sample+'.mtz')

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.emit(QtCore.SIGNAL("finished()"))



class read_autoprocessing_results_from_disc(QtCore.QThread):
    def __init__(self,visit_list,target,reference_file_list,database_directory):
        QtCore.QThread.__init__(self)
        self.visit_list=visit_list
        self.target=target
        self.reference_file_list=reference_file_list
        self.data_collection_dict={}
        self.data_collection_statistics_dict={}
        self.database_directory=database_directory
        self.data_collection_dict_collected={}
        self.data_collection_statistics_dict_collected={}

        if os.path.isfile(os.path.join(self.database_directory,'data_collection_summary.pkl')):
            summary = pickle.load( open( os.path.join(self.database_directory,'data_collection_summary.pkl'), "rb" ) )
            self.data_collection_dict_collected=summary[0]
            self.data_collection_statistics_dict_collected=summary[1]

    def run(self):
        number_of_visits_to_search=len(self.visit_list)
        search_cycle=1
        for visit_directory in sorted(self.visit_list):
            if len(glob.glob(os.path.join(visit_directory,'processed',self.target,'*')))==0:
                continue
            progress_step=100/float(len(glob.glob(os.path.join(visit_directory,'processed',self.target,'*'))))
            progress=0

            if 'attic' in visit_directory:
                visit=visit_directory.split('/')[6]
            else:
                visit=visit_directory.split('/')[5]

            ####################################################################################################
            # dewar configuration:
            dewar_configuration=[]
            for xml in glob.glob(os.path.join(visit_directory,'xml','exptTableParams-*.xml')):
                prefix=''
                container_reference=''
                sample_location=''
                for line in open(xml):
                    if 'prefix' in line:
                        prefix=line[line.find('>')+1:line.rfind('<')]
                    if 'container_reference' in line:
                        container_reference=line[line.find('>')+1:line.rfind('<')]
                    if 'sample_location' in line:
                        sample_location=line[line.find('>')+1:line.rfind('<')]
                dewar_configuration.append([prefix,container_reference,sample_location])

            for collected_xtals in sorted(glob.glob(os.path.join(visit_directory,'processed',self.target,'*'))):
                # this step is only relevant when several samples are reviewed in one session
                if 'tmp' in collected_xtals or 'results' in collected_xtals:
                    continue
                xtal=collected_xtals[collected_xtals.rfind('/')+1:]
                protein_name=collected_xtals.split('/')[len(collected_xtals.split('/'))-2]
                self.data_collection_dict[xtal]=[[],[],[],[],[]]
                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'Step 1 of 3: searching visit '+ \
                                                                       str(search_cycle)+' of '+str(number_of_visits_to_search)+ \
                                                                       ' ('+visit+'/'+xtal+')')
                run_list=[]
                logfile_list=[]
                image_list=[]
                image_string_list=[]
                puck_position=[]
                for runs in glob.glob(collected_xtals+'/*'):
                    run_is_in_pickle_file=False
                    run=runs[runs.rfind('/')+1:]
                    for xtals_in_collected_dict in self.data_collection_dict_collected:
                        if xtals_in_collected_dict==xtal:
                            for runs_in_collected_dict in self.data_collection_dict_collected[xtals_in_collected_dict][0]:
                                if runs_in_collected_dict[0]==run:
                                    run_is_in_pickle_file=True
                    if run_is_in_pickle_file:
                        for stuff in self.data_collection_dict_collected[xtal][0]:
                            if stuff[0]==run:
                                self.data_collection_dict[xtal][0].append(stuff)
                        for stuff in self.data_collection_dict_collected[xtal][1]:
                            if run in stuff:
                                image_list.append(stuff)
                        for stuff in self.data_collection_dict_collected[xtal][2]:
                            if run in stuff:
                                logfile_list.append(stuff)
                        for stuff in self.data_collection_dict_collected[xtal][3]:
                            if run in stuff[0]:
                                image_string_list.append(stuff)
                        if len(self.data_collection_dict_collected[xtal])==5:
                            puck_position=self.data_collection_dict_collected[xtal][4]
                        continue

                    timestamp=datetime.fromtimestamp(os.path.getmtime(runs)).strftime('%Y-%m-%d %H:%M:%S')
                    diffraction_image_directory=os.path.join(visit_directory,protein_name,xtal)
                    run_list.append([(run,timestamp,visit,diffraction_image_directory)])
                    self.data_collection_dict[xtal][0].append([run,timestamp,visit,diffraction_image_directory])

                    for file_name in glob.glob(os.path.join(runs,'xia2','*','LogFiles','*')):
                        if file_name.endswith('aimless.log') and (self.target in file_name or self.target=='*'):
                            logfile_list.append(file_name)
                    for file_name in glob.glob(os.path.join(runs,'fast_dp','*')):
                        if file_name.endswith('aimless.log') and (self.target in file_name or self.target=='*'):
                            logfile_list.append(file_name)
                    for file_name in glob.glob(os.path.join(runs,'autoPROC','ap-run','*')):
                        if file_name.endswith('aimless.log') and (self.target in file_name or self.target=='*'):
                            logfile_list.append(file_name)



                    for image in glob.glob(os.path.join(visit_directory,'jpegs',self.target,xtal,'*')):
                        if run in image:
                            if image.endswith('t.png') or image.endswith('_.png'):
                                image_list.append(image)
                                image_file=open(image,"rb")
                                image_string=base64.b64encode(image_file.read())
                                image_string_list.append((image[image.rfind('/')+1:],image_string))

                    for item in dewar_configuration:
                        if item[0]==xtal:
                            puck_position=[item[1],item[2]]
                            break

                self.data_collection_dict[xtal][1]+=image_list
                self.data_collection_dict[xtal][2]+=logfile_list
                self.data_collection_dict[xtal][3]+=image_string_list
                self.data_collection_dict[xtal][4]=puck_position
                progress += progress_step
                self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

            search_cycle+=1


        if not len(self.data_collection_dict)==0:
            progress_step=100/float(len(self.data_collection_dict))
        progress=0
        for sample in sorted(self.data_collection_dict):
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'Step 2 of 3: parsing aimless logfiles -> '+sample)
            self.data_collection_statistics_dict[sample]=[]
#            print sample,self.data_collection_dict[sample][2]
            if not self.data_collection_dict[sample][2]==[]:
                for index,logfile in enumerate(self.data_collection_dict[sample][2]):
                    already_parsed=False
                    if sample in self.data_collection_statistics_dict_collected:
                        for entry in self.data_collection_statistics_dict_collected[sample]:
                            if len(entry) > 2:
                                if logfile==entry[1]:
                                    self.data_collection_statistics_dict[sample].append(entry)
                                    already_parsed=True
                    if already_parsed:
                        continue
                    else:
                        aimless_results=parse().GetAimlessLog(logfile)
                        try:
                            self.data_collection_statistics_dict[sample].append([
                        index,                                                                                      # 0
                        logfile,                                                                                    # 1
                        ['Program',                     aimless_results['AutoProc'],                                                        (200,200,200)],
                        ['Run',                         aimless_results['Run'],                                                             (200,200,200)],
                        ['Space\nGroup',                aimless_results['SpaceGroup'],                                                      (200,200,200)],
                        ['Unit Cell',                   aimless_results['UnitCell'],                                                        (200,200,200)],
                        ['Resolution\nOverall',         aimless_results['ResolutionLow']+'-'+aimless_results['ResolutionHigh'],             (200,200,200)],
                        ['Resolution\nInner Shell',     aimless_results['ResolutionLow']+'-'+aimless_results['ResolutionLowInnerShell'],    (200,200,200)],
                        ['Resolution\nOuter Shell',     aimless_results['ResolutionHighOuterShell']+'-'+aimless_results['ResolutionHigh'],  (200,200,200)],
                        ['Rmerge\nOverall',             aimless_results['RmergeOverall'],                                                   (200,200,200)],
                        ['Rmerge\nInner Shell',         aimless_results['RmergeLow'],                                                       (200,200,200)],
                        ['Rmerge\nOuter Shell',         aimless_results['RmergeHigh'],                                                      (200,200,200)],
                        ['Mn(I/sig(I))\nOverall',       aimless_results['IsigOverall'],                                                     (200,200,200)],
                        ['Mn(I/sig(I))\nInner Shell',   aimless_results['IsigLow'],                                                         (200,200,200)],
                        ['Mn(I/sig(I))\nOuter Shell',   aimless_results['IsigHigh'],                                                        (200,200,200)],
                        ['Completeness\nOverall',       aimless_results['CompletenessOverall'],                                             (200,200,200)],
                        ['Completeness\nInner Shell',   aimless_results['CompletenessLow'],                                                 (200,200,200)],
                        ['Completeness\nOuter Shell',   aimless_results['CompletenessHigh'],                                                (200,200,200)],
                        ['Multiplicity\nOverall',       aimless_results['MultiplicityOverall'],                                             (200,200,200)],
                        ['Multiplicity\nInner Shell',   aimless_results['MultiplicityLow'],                                                 (200,200,200)],
                        ['Multiplicity\nOuter Shell',   aimless_results['MultiplicityHigh'],                                                (200,200,200)],
                        aimless_results['Lattice'],                                                                 # 21
                        float(aimless_results['UniqueReflectionsOverall']),                                         # 22
                        float(aimless_results['CompletenessOverall']),                                              # 23
                        float(aimless_results['IsigOverall']),                                                      # 24
                        float(aimless_results['UnitCellVolume']),                                                   # 25
                        float(aimless_results['RmergeLow']),                                                        # 26
                        ['best file',False]                                                                                       # 27
                                        ])
                        except ValueError:
                            self.data_collection_statistics_dict[sample].append([
                        index,                                                                                      # 0
                        logfile,                                                                                    # 1
                        ['Program',                     aimless_results['AutoProc'],                                                        (200,200,200)],
                        ['Run',                         aimless_results['Run'],                                                             (200,200,200)],
                        ['Space\nGroup',                aimless_results['SpaceGroup'],                                                      (200,200,200)],
                        ['Unit Cell',                   aimless_results['UnitCell'],                                                        (200,200,200)],
                        ['Resolution\nOverall',         aimless_results['ResolutionLow']+'-'+aimless_results['ResolutionHigh'],             (200,200,200)],
                        ['Resolution\nInner Shell',     aimless_results['ResolutionLow']+'-'+aimless_results['ResolutionLowInnerShell'],    (200,200,200)],
                        ['Resolution\nOuter Shell',     aimless_results['ResolutionHighOuterShell']+'-'+aimless_results['ResolutionHigh'],  (200,200,200)],
                        ['Rmerge\nOverall',             aimless_results['RmergeOverall'],                                                   (200,200,200)],
                        ['Rmerge\nInner Shell',         aimless_results['RmergeLow'],                                                       (200,200,200)],
                        ['Rmerge\nOuter Shell',         aimless_results['RmergeHigh'],                                                      (200,200,200)],
                        ['Mn(I/sig(I))\nOverall',       aimless_results['IsigOverall'],                                                     (200,200,200)],
                        ['Mn(I/sig(I))\nInner Shell',   aimless_results['IsigLow'],                                                         (200,200,200)],
                        ['Mn(I/sig(I))\nOuter Shell',   aimless_results['IsigHigh'],                                                        (200,200,200)],
                        ['Completeness\nOverall',       aimless_results['CompletenessOverall'],                                             (200,200,200)],
                        ['Completeness\nInner Shell',   aimless_results['CompletenessLow'],                                                 (200,200,200)],
                        ['Completeness\nOuter Shell',   aimless_results['CompletenessHigh'],                                                (200,200,200)],
                        ['Multiplicity\nOverall',       aimless_results['MultiplicityOverall'],                                             (200,200,200)],
                        ['Multiplicity\nInner Shell',   aimless_results['MultiplicityLow'],                                                 (200,200,200)],
                        ['Multiplicity\nOuter Shell',   aimless_results['MultiplicityHigh'],                                                (200,200,200)],
                        aimless_results['Lattice'],                                                                 # 21
                        0.0,                                         # 22
                        0.0,                                              # 23
                        0.0,                                                      # 24
                        0.0,                                                   # 25
                        100.0,                                                        # 26
                        ['best file',False]                                                                                       # 27
                                        ])

            else:
                self.data_collection_statistics_dict[sample]+='###'*27
            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        # before creating the table with the results, try to guess which one to select
        # 1. check if there are reference mtz files
        # 1a. if so: take all logfiles forward that fit to the first one found
        #     'fit means': same lattice and delta Vunitcell < 5%
        # 2. if possible: select all datasets with Rmerge low < 5%
        # 3. finally select the dataset with
        #    max(unique_reflections*completeness*Mn(I/sig<I>)

        if not len(self.data_collection_statistics_dict)==0:
            progress_step=100/float(len(self.data_collection_statistics_dict))
        progress=0

        # if possible, select only the ones which have the same lattice and
        # a unit cell volume difference of less than 5%
        for sample in sorted(self.data_collection_statistics_dict):
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'Step 3 of 3: selecting "best" aimless logfile ->'+sample)
            if self.data_collection_statistics_dict[sample][0]=='#':
                continue
            select_stage_one_list = []
            found=0
            if self.reference_file_list != []:
                for reference_file in self.reference_file_list:
                    for aimless_file in self.data_collection_statistics_dict[sample]:
                        try:
                            if not reference_file[4]==0:
                                unitcell_difference=round((math.fabs(reference_file[4]-aimless_file[25])/reference_file[4])*100,1)
                                if unitcell_difference < 5 and reference_file[3]==aimless_file[21]:
                                    select_stage_one_list.append(aimless_file)
                                    found=1
                        except IndexError:
                            pass
                if not found:                                                   # in case no file fullfils the criteria
                    for aimless_file in self.data_collection_statistics_dict[sample]:
                        if aimless_file != []:
                            select_stage_one_list.append(aimless_file)
            else:                                                               # in case no reference files are available
                for aimless_file in self.data_collection_statistics_dict[sample]:
                    if aimless_file != []:
                        select_stage_one_list.append(aimless_file)

            # if possible, select only the ones with Rmerge < 5%
            select_stage_two_list=[]
            for aimless_file in select_stage_one_list:
                if aimless_file[26] < 0.05:
                    select_stage_two_list.append(aimless_file)
            if select_stage_two_list==[]:
                select_stage_two_list=select_stage_one_list

            # finally, select the file with the highest
            # max(unique_reflections*completeness*Mn(I/sig<I>)
            select_stage_three_list=[]
            for aimless_file in select_stage_two_list:
                select_stage_three_list.append([aimless_file[0],
                                                aimless_file[22] \
                                                * aimless_file[23] \
                                                * aimless_file[24]])
            if select_stage_three_list != []:
                best_file_index=max(select_stage_three_list,key=lambda x: x[1])[0]
                for index,results in enumerate(self.data_collection_statistics_dict[sample]):
                    if index==best_file_index:
                        self.data_collection_statistics_dict[sample][index][27]=['best file',True]
            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        # save everything so that it's quicker to reload and is available outside DLS
        pickle.dump([self.data_collection_dict,self.data_collection_statistics_dict],
                    open(  os.path.join(self.database_directory,'data_collection_summary.pkl'),'wb'))

        self.emit(QtCore.SIGNAL('create_widgets_for_autoprocessing_results'), [self.data_collection_dict,
                                                                            self.data_collection_statistics_dict])





class NEW_read_autoprocessing_results_from_disc(QtCore.QThread):
    def __init__(self,visit_list,target,reference_file_list,database_directory,data_collection_dict):
        QtCore.QThread.__init__(self)
        self.visit_list=visit_list
        self.target=target
        self.reference_file_list=reference_file_list
        self.data_collection_dict=data_collection_dict
        self.database_directory=database_directory


        # - open data source if possible
        # - get sampleID, xtbm
        # - search lab36 folder for respective xtal image
        # - convert to string and use in data dict
        # - but images can only be found of XCE is started in the respective labchem directory

    def run(self):

        # only do once, ignore if just refreshing table
        if self.data_collection_dict=={}:
            if os.path.isfile(os.path.join(self.database_directory,'test.pkl')):
                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'unpickling: '+os.path.join(self.database_directory,'test.pkl'))
                self.data_collection_dict = pickle.load( open( os.path.join(self.database_directory,'test.pkl'), "rb" ) )

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
                if 'tmp' in collected_xtals or 'results' in collected_xtals:
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
                for runs in sorted(glob.glob(collected_xtals+'/*')):
                    run=runs[runs.rfind('/')+1:]
                    diffraction_image=''
                    timestamp=datetime.fromtimestamp(os.path.getmtime(runs)).strftime('%Y-%m-%d %H:%M:%S')
                    if os.path.isfile(os.path.join(visit_directory,protein_name,xtal,run+'0001.cbf')):
                        diffraction_image=os.path.join(visit_directory,protein_name,xtal,run+'0001.cbf')

                    # image files
                    # note: need one more flag which indicates immediately that images belong together
                    #       this makes it afterwards easier to get them together in the table
                    run_number_list=[]
                    image_files_in_list=False
                    for entry in self.data_collection_dict[xtal]:
                        image_files_in_list=False
                        if entry[0]=='image':
                            print 'ok, there is the dict, should not be there during first run'
                            print entry[0],entry[1],entry[2],entry[3],entry[5],entry[6]
                            print run
                            if entry[0]=='image' and entry[1]==visit and entry[2]==run:
                                print 'already in dict:',run,visit
                                image_files_in_list=True
                            if not image_files_in_list:
                                print 'missing:',run,visit
                                run_number_list.append(int(entry[6])+1)
                                print run_number_list
                    if run_number_list==[]:
                        run_number=1
                    else:
                        run_number=max(run_number_list)
                        print 'run_number:',run_number
                    if not image_files_in_list:
                        image_list=[]
                        for image in glob.glob(os.path.join(visit_directory,'jpegs',self.target,xtal,'*')):
                            if run in image:
                                if image.endswith('t.png') or image.endswith('_.png'):
                                    image_name=image[image.rfind('/')+1:]
                                    image_file=open(image,"rb")
                                    image_string=base64.b64encode(image_file.read())
                                    image_list.append( [image_name,image_string] )
                        self.data_collection_dict[xtal].append(['image',visit,run,timestamp,image_list,
                                                                diffraction_image,run_number])


                    ##########################################################################
                    # aimless & Dimple information
                    # first for xia2 runs
                    for file_name in glob.glob(os.path.join(visit_directory,'processed',protein_name,xtal,run,'xia2','*','LogFiles','*aimless.log')):
                        db_dict={   'DataCollectionVisit':          visit,
                                    'DataCollectionBeamline':       beamline,
                                    'DataCollectionDate':           timestamp,
                                    'DataProcessingPathToLogfile':  file_name   }
                        autoproc=file_name.split('/')[len(file_name.split('/'))-3]
                        found_autoproc=False
                        for entry in self.data_collection_dict[xtal]:
                            if len(entry)>=9:
                                if entry[0]=='logfile' and entry[1]==visit and entry[2]==run and entry[4]==autoproc:
                                    found_autoproc=True
                        if not found_autoproc:
                            aimless_results=parse().read_aimless_logfile(file_name)
                            db_dict.update(aimless_results)
                            if os.path.isfile(os.path.join(visit_directory,'processed',protein_name,xtal,run,'xia2',autoproc,'dimple','final.pdb')):
                                dimple_file=os.path.join(visit_directory,'processed',protein_name,xtal,run,'xia2',autoproc,'dimple','final.pdb')
                                pdb_info=parse().PDBheader(dimple_file)
                                db_dict['DataProcessingRcryst']  = pdb_info['Rcryst']
                                db_dict['DataProcessingRfree'] = pdb_info['Rfree']
                            else:
                                db_dict['DataProcessingRcryst']  = ''
                                db_dict['DataProcessingRfree'] = ''
                            self.data_collection_dict[xtal].append(['logfile',visit,run,timestamp,autoproc,file_name,db_dict,0,False])




                    # then exactly the same for fast_dp
                    if os.path.isfile(os.path.join(runs,'fast_dp','aimless.log')):
                        db_dict={   'DataCollectionVisit':          visit,
                                    'DataCollectionBeamline':       beamline,
                                    'DataCollectionDate':           timestamp,
                                    'DataProcessingPathToLogfile':  file_name   }
                        file_name=os.path.join(runs,'fast_dp','aimless.log')
                        autoproc=file_name.split('/')[len(file_name.split('/'))-2]
                        found_autoproc=False
                        for entry in self.data_collection_dict[xtal]:
                            if len(entry)>=9:
                                if entry[0]=='logfile' and entry[1]==visit and entry[2]==run and entry[4]==autoproc:
                                    found_autoproc=True
                        if not found_autoproc:
                            aimless_results=parse().read_aimless_logfile(file_name)
                            db_dict.update(aimless_results)
                            if os.path.isfile(os.path.join(runs,'fast_dp','dimple','final.pdb')):
                                dimple_file=os.path.join(runs,'fast_dp','dimple','final.pdb')
                                pdb_info=parse().PDBheader(dimple_file)
                                db_dict['DataProcessingRcryst']  = pdb_info['Rcryst']
                                db_dict['DataProcessingRfree'] = pdb_info['Rfree']
                            else:
                                db_dict['DataProcessingRcryst']  = ''
                                db_dict['DataProcessingRfree'] = ''
                            self.data_collection_dict[xtal].append(['logfile',visit,run,timestamp,autoproc,file_name,db_dict,0,False])

                    # then exactly the same for autoPROC
                    if os.path.isfile(os.path.join(runs,'autoPROC','ap-run','aimless.log')):
                        db_dict={   'DataCollectionVisit':          visit,
                                    'DataCollectionBeamline':       beamline,
                                    'DataCollectionDate':           timestamp,
                                    'DataProcessingPathToLogfile':  file_name   }
                        file_name=os.path.join(runs,'autoPROC','ap-run','aimless.log')
                        autoproc=file_name.split('/')[len(file_name.split('/'))-3]
                        found_autoproc=False
                        for entry in self.data_collection_dict[xtal]:
                            if len(entry)>=9:
                                if entry[0]=='logfile' and entry[1]==visit and entry[2]==run and entry[4]==autoproc:
                                    found_autoproc=True
                        if not found_autoproc:
                            aimless_results=parse().read_aimless_logfile(file_name)
                            db_dict.update(aimless_results)
                            if os.path.isfile(os.path.join(runs,'autoPROC','dimple','final.pdb')):
                                dimple_file=os.path.join(runs,'autoPROC','dimple','final.pdb')
                                pdb_info=parse().PDBheader(dimple_file)
                                db_dict['DataProcessingRcryst']  = pdb_info['Rcryst']
                                db_dict['DataProcessingRfree'] = pdb_info['Rfree']
                            else:
                                db_dict['DataProcessingRcryst']  = ''
                                db_dict['DataProcessingRfree'] = ''
                            self.data_collection_dict[xtal].append(['logfile',visit,run,timestamp,autoproc,file_name,db_dict,0,False])



                progress += progress_step
                self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

            search_cycle+=1

        # now we try, somewhat haphazardly, to determine which autoprocessing run is currently the best

        # before creating the table with the results, try to guess which one to select
        # 1. check if there are reference mtz files
        # 1a. if so: take all logfiles forward that fit to the first one found
        #     'fit means': same lattice and delta Vunitcell < 5%
        # 2. if possible: select all datasets with Rmerge low < 5%
        # 3. finally select the dataset with
        #    max(unique_reflections*completeness*Mn(I/sig<I>)

        if not len(self.data_collection_dict)==0:
            progress_step=100/float(len(self.data_collection_dict))
        progress=0


        for xtal in sorted(self.data_collection_dict):
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'Step 2 of 2: selecting "best" aimless logfile ->'+xtal)
#            print 'hello im here now'
            ############################################################################################
            # STAGE 1:
            # similarity to reference files
            select_stage_one_list = []
            tmp=[]
            index=0
            for n,entry in enumerate(self.data_collection_dict[xtal]):
#                print entry[0]
                found=False
                if len(entry)==9 and entry[0]=='logfile':
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
                    self.data_collection_dict[xtal][n][7]=index
                    index+=1

 #           print select_stage_one_list
            # if none passed Stage 1, carry them over to Stage 2
            if select_stage_one_list == [] and tmp != []:
                select_stage_one_list=tmp

            ############################################################################################
            # STAGE 2:
            # if possible, select only the ones with Rmerge < 5%
            select_stage_two_list=[]
            tmp=[]
            for index in select_stage_one_list:
                # this may be completely over the top!
                found=False
                for entry in self.data_collection_dict[xtal]:
                    if len(entry)==9 and entry[0]=='logfile':
                        if isinstance(entry[6],dict):
                            try:
                                if float(entry[6]['DataProcessingRmergeLow']) < 0.05 and entry[7]==index:
                                    select_stage_two_list.append(index)
                                    found=True
                            except ValueError:
                                pass
                if not found:
                    tmp.append(index)

            # if none passed Stage 2, carry them over to Stage 3
            if select_stage_two_list == [] and tmp != []:
                select_stage_two_list=tmp
#            print select_stage_two_list

            ############################################################################################
            # STAGE 3:
            # finally, select the file with the highest
            # max(unique_reflections*completeness*Mn(I/sig<I>)
            select_stage_three_list=[]
            for index in select_stage_two_list:
                for entry in self.data_collection_dict[xtal]:
#                    print len(entry),entry
                    if len(entry)>=9 and entry[0]=='logfile':
#                        print 'ok first requirement fullfilled'
                        if isinstance(entry[6],dict) and entry[7]==index:
#                            print 'ok second requirement fullfilled'
#                            print 'ur',entry[6]['UniqueReflectionsOverall']
#                            print 'comp',entry[6]['CompletenessOverall']
#                            print 'isg',entry[6]['IsigOverall']
                            try:
                                ranking=float(entry[6]['DataProcessingUniqueReflectionsOverall'])*\
                                        float(entry[6]['DataProcessingCompletenessOverall'])*\
                                        float(entry[6]['DataProcessingIsigOverall'])
                                select_stage_three_list.append([index,ranking])
                            except ValueError:
                                pass

            if not select_stage_three_list==[]:
                best_file_index=max(select_stage_three_list,key=lambda x: x[1])[0]
                for n,entry in enumerate(self.data_collection_dict[xtal]):
                    if len(entry)==9 and entry[0]=='logfile':
                        if entry[7]==best_file_index:
                            self.data_collection_dict[xtal][n][8]=True


        # save everything so that it's quicker to reload and is available outside DLS
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'pickling results')
        pickle.dump(self.data_collection_dict,open(  os.path.join(self.database_directory,'test.pkl'),'wb'))

        self.emit(QtCore.SIGNAL('create_widgets_for_autoprocessing_results'), self.data_collection_dict)
