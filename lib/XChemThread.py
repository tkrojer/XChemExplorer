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
            helpers().make_png(self.initial_model_directory,sampleID,compoundID,smiles)
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
        if len(self.initial_model_dimple_dict) != 0:
            progress_step=100/float(len(self.initial_model_dimple_dict))
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


#class cluster_datasets(QtCore.QThread):
#    def __init__(self,settings,initial_model_dimple_dict,external_software,ccp4_scratch,filename_root):
#        QtCore.QThread.__init__(self)
#        self.initial_model_directory=settings['initial_model_directory']
#        self.reference_directory=settings['reference_directory']
#        self.initial_model_dimple_dict=initial_model_dimple_dict
#        self.queueing_system_available=external_software['qsub']
#        self.ccp4_scratch_directory=ccp4_scratch
#        self.filename_root=filename_root
#
#    def run(self):
##       list=[[sampleID,point_group,unitcell_volume],
##             [....................................]]
##         item[0] is member of the first arbitrary cluster
##       name of cluster: 1-<point_group>
##        for item in list:
#
##
##        if len(self.initial_model_dimple_dict) != 0:
##            progress_step=100/float(len(self.initial_model_dimple_dict))
##        progress=0
##
#        for sample in sorted(self.initial_model_dimple_dict):
#            if self.initial_model_dimple_dict[sample][0].isChecked():
##                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'running dimple -> '+sample)
##                dimple_commands={   'project_directory': self.initial_model_directory,
##                                    'delete_old': self.initial_model_dimple_dict[sample][0].isChecked(),
##                                    'xtalID': sample,
##                                    'compoundID': '',
##                                    'smiles': '',
##                                    'reference': self.reference_directory+'/'+
##                                                 str(self.initial_model_dimple_dict[sample][1].currentText()),
##                                    'queueing_system_available': self.queueing_system_available,
##                                    'ccp4_scratch': self.ccp4_scratch_directory,
##                                    'fileroot_in':  self.filename_root.replace('${samplename}',sample)  }
##                process(dimple_commands).dimple()
##            progress += progress_step
##            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
##        self.emit(QtCore.SIGNAL("finished()"))
##---------------------------------------------------------------------------------------------------------------



class start_COOT(QtCore.QThread):

    def __init__(self,settings):
        QtCore.QThread.__init__(self)
        self.settings=settings

    def run(self):
        pickle.dump(self.settings,open('XChemExplorer_settings.pkl','wb'))
        os.system('coot --no-guano --no-state-script --script %s' %(os.getenv('XChemExplorer_DIR')+'/lib/XChemCoot.py'))


class read_intial_refinement_results(QtCore.QThread):

    def __init__(self,initial_model_directory,reference_file_list,data_source,allowed_unitcell_difference_percent,filename_root):
        QtCore.QThread.__init__(self)
        self.initial_model_directory=initial_model_directory
        self.reference_file_list=reference_file_list
        self.data_source=data_source
        self.allowed_unitcell_difference_percent=allowed_unitcell_difference_percent
        self.filename_root=filename_root

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
#            sample_mtz=self.filename_root+'.mtz'
#            if os.path.isfile(os.path.join(self.initial_model_directory,sample,sample+'.free.pdb')):
#                tmp=os.path.join(self.initial_model_directory,sample,sample+'.free.pdb')
#                sample_mtz=tmp[tmp.rfind('/')+1:]
#            elif os.path.isfile(os.path.join(self.initial_model_directory,sample,sample+'-pandda-input.mtz')):
#                tmp=os.path.join(self.initial_model_directory,sample,sample+'-pandda-input.mtz')
#                sample_mtz=tmp[tmp.rfind('/')+1:]
#            elif os.path.isfile(os.path.join(self.initial_model_directory,sample,'refine.mtz')):
#                tmp=os.path.join(self.initial_model_directory,sample,'refine.mtz')
#                sample_mtz=tmp[tmp.rfind('/')+1:]

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
            if os.path.isfile(os.path.join(self.initial_model_directory,sample,'refine.pdb')):
                pdb_latest=os.path.realpath(os.path.join(self.initial_model_directory,sample,'refine.pdb'))

            mtz_latest=''
            if os.path.isfile(os.path.join(self.initial_model_directory,sample,'refine.mtz')):
                mtz_latest=os.path.realpath(os.path.join(self.initial_model_directory,sample,'refine.mtz'))

            mtz_free=''
            if os.path.isfile(os.path.join(self.initial_model_directory,sample,sample+'free.mtz')):
                mtz_free=os.path.realpath(os.path.join(self.initial_model_directory,sample,sample+'free.mtz'))


            # update data source only if field is null
            db_dict={   'RefinementOutcome':        outcome,
                        'RefinementMTZfree':        mtz_free    }
            db.update_insert_not_null_fields_only(sample,db_dict)

            # update to reflect current state
            db_dict={   'RefinementRcryst':         Rcryst,
                        'RefinementRfree':          Rfree,
                        'RefinementRmsdBonds':      rmsdBonds,
                        'RefinementRmsdAngles':     rmsdAngles,
                        'RefinementPDB_latest':     pdb_latest,
                        'RefinementMTZ_latest':     mtz_latest          }
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
            #data_collection_dict = pickle.load( open( os.path.join(self.database_directory,'data_collection_summary.pkl'), "rb" ) )
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
                            if logfile==entry[1]:
                                self.data_collection_statistics_dict[sample].append(entry)
                                already_parsed=True
                    if already_parsed:
#                        print 'hallo'
                        continue
                    else:
#                        print 'hallo ',logfile,self.data_collection_dict[sample][2]
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
