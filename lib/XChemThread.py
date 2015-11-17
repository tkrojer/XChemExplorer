import os, sys, glob
from datetime import datetime
from PyQt4 import QtGui, QtCore
#from PyQt4.QtCore import QThread, SIGNAL

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


class run_dimple_on_selected_samples(QtCore.QThread):
    def __init__(self,settings,initial_model_dimple_dict,queueing_system_available):
        QtCore.QThread.__init__(self)
        self.initial_model_directory=settings['initial_model_directory']
        self.reference_directory=settings['reference_directory']
        self.initial_model_dimple_dict=initial_model_dimple_dict
        self.queueing_system_available=queueing_system_available

    def run(self):
        if len(self.initial_model_dimple_dict) != 0:
            progress_step=100/float(len(self.initial_model_dimple_dict))
        progress=0

        for sample in sorted(self.initial_model_dimple_dict):
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'running dimple -> '+sample)
#            print sample,self.initial_model_dimple_dict[sample][0].isChecked(),\
#                str(self.initial_model_dimple_dict[sample][1].currentText())
            dimple_commands={   'project_directory': self.initial_model_directory,
                                'delete_old': self.initial_model_dimple_dict[sample][0].isChecked(),
                                'xtalID': sample,
                                'compoundID': '',
                                'smiles': '',
                                'reference': self.reference_directory+'/'+
                                             str(self.initial_model_dimple_dict[sample][1].currentText()),
                                'queueing_system_available': self.queueing_system_available }
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


class read_intial_refinement_results(QtCore.QThread):

    def __init__(self,initial_model_directory,reference_file_list):
        QtCore.QThread.__init__(self)
        self.initial_model_directory=initial_model_directory
        self.reference_file_list=reference_file_list

    def run(self):

        progress_step=100/float(len(glob.glob(self.initial_model_directory+'/*')))
        progress=0

        initial_model_list=[]

        for sample_dir in sorted(glob.glob(self.initial_model_directory+'/*')):

            sample=sample_dir[sample_dir.rfind('/')+1:]
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'parsing initial_models folder -> '+sample)

            run_dimple=True
            resolution_high=''
            Rcryst='pending'
            Rfree='pending'
            spg_autoproc=''
            unitcell_autoproc=''
            spg_reference=''
            unitcell_reference=''
            unitcell_difference=''
            reference=''
            alert='#E0E0E0'

            if os.path.isfile(os.path.join(self.initial_model_directory,sample,sample+'.mtz')):
                mtz_autoproc=mtztools(os.path.join(self.initial_model_directory,sample,sample+'.mtz')).get_all_values_as_dict()
                resolution_high=mtz_autoproc['resolution_high']
                spg_autoproc=mtz_autoproc['spacegroup']
                unitcell_autoproc=mtz_autoproc['unitcell']
                lattice_autoproc=mtz_autoproc['bravais_lattice']
                unitcell_volume_autoproc=mtz_autoproc['unitcell_volume']
                # check which reference file is most similar
                for o,reference_file in enumerate(self.reference_file_list):
                    unitcell_difference=round((math.fabs(reference_file[4]-unitcell_volume_autoproc)/reference_file[4])*100,1)
                    # reference file is accepted when different in unitcell volume < 5%
                    # and both files have the same lattice type
                    if unitcell_difference < 5 and lattice_autoproc==reference_file[3]:
                        spg_reference=reference_file[1]
                        unitcell_reference=reference_file[2]
                        reference=reference_file[0]
                        break
            if os.path.isdir(os.path.join(self.initial_model_directory,sample,'Dimple')):
                    if os.path.isfile(os.path.join(self.initial_model_directory,sample,'Dimple','dimple','final.pdb')):
                        pdb=parse().PDBheader(os.path.join(self.initial_model_directory,sample,'Dimple','dimple','final.pdb'))
                        Rcryst=pdb['Rcryst']
                        Rfree=pdb['Rfree']
                        alert=pdb['Alert']
                    elif os.path.isfile(os.path.join(self.initial_model_directory,sample,'/dimple_run_in_progress')):
                        Rcryst='in progress'
                        Rfree='in progress'
                        alert=(51,153,255)

            initial_model_list.append( [ sample,
                                        run_dimple,
                                        resolution_high,
                                        Rcryst,
                                        Rfree,
                                        spg_autoproc,
                                        spg_reference,
                                        unitcell_difference,
                                        reference,
                                        unitcell_autoproc,
                                        unitcell_reference,
                                        alert ] )

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
        self.emit(QtCore.SIGNAL('create_initial_model_table'), initial_model_list)



class save_autoprocessing_results_to_disc(QtCore.QThread):
    def __init__(self,dataset_outcome_dict,data_collection_table_dict,data_collection_statistics_dict):
        QtCore.QThread.__init__(self)
        self.dataset_outcome_dict=dataset_outcome_dict
        self.data_collection_table_dict=data_collection_table_dict
        self.data_collection_statistics_dict=data_collection_statistics_dict

    def run(self):
        if not len(self.dataset_outcome_dict)==0:
            progress_step=100/float(len(self.dataset_outcome_dict))
        progress=0
        csv_out=''
        for key in sorted(self.dataset_outcome_dict):
            outcome=''
            for button in self.dataset_outcome_dict[key]:
                if button.isChecked():
                    print key,button.text()
                    outcome=button.text()
            csv_out+=str(key)+','+outcome+','
            if outcome=='success':
                indexes=self.data_collection_table_dict[key].selectionModel().selectedRows()
                for index in sorted(indexes):
                    # csv out
                    for item in self.data_collection_statistics_dict[key][index.row()]:
                        csv_out+=str(item)+','
                    csv_out+='\n'

#                    print self.data_collection_table_dict[key]
                    print self.data_collection_statistics_dict[key][index.row()]
                    print self.data_collection_statistics_dict[key][index.row()]

#if not os.path.isdir(os.path.join(self.initial_model_directory,key)):
#    os.mkdir(os.path.join(self.initial_model_directory,key))
                    os.symlink('0-coot.state.scm','temp.link')
#if not os.path.isdir(os.path.join(self.initial_model_directory,key,'autoprocessing')):
#    os.mkdir(os.path.join(self.initial_model_directory,key,'autoprocessing'))


                    if 'xia2' in self.data_collection_statistics_dict[key][index.row()][1]:
#                        print self.data_collection_statistics_dict[key][index.row()][1]
                        print os.path.join(*self.data_collection_statistics_dict[key][index.row()][1].split('/')[:13])
#                print('Row %d is selected' % index.row())
                    if 'fast_dp' in self.data_collection_statistics_dict[key][index.row()][1]:
#                        print self.data_collection_statistics_dict[key][index.row()][1]
                        print os.path.join(*self.data_collection_statistics_dict[key][index.row()][1].split('/')[:12])


            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'writing files from data processing to inital_model folder -> '+key)

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
        print csv_out

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
                break
            progress_step=100/float(len(glob.glob(os.path.join(visit_directory,'processed',self.target,'*'))))
            progress=0
            if 'attic' in visit_directory:
                visit=visit_directory.split('/')[6]
            else:
                visit=visit_directory.split('/')[5]
            for collected_xtals in sorted(glob.glob(os.path.join(visit_directory,'processed',self.target,'*'))):
                xtal=collected_xtals[collected_xtals.rfind('/')+1:]
                self.data_collection_dict[xtal]=[[],[],[],[]]
                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'Step 1 of 3: searching visit '+ \
                                                                       str(search_cycle)+' of '+str(number_of_visits_to_search)+ \
                                                                       ' ('+visit+'/'+xtal+')')
                run_list=[]
                logfile_list=[]
                image_list=[]
                image_string_list=[]
                for runs in glob.glob(collected_xtals+'/*'):
                    run_is_in_pickle_file=False
                    run=runs[runs.rfind('/')+1:]
                    for xtals_in_collected_dict in self.data_collection_dict_collected:
                        if xtals_in_collected_dict==xtal:
                            for runs_in_collected_dict in self.data_collection_dict_collected[xtals_in_collected_dict][0]:
                                if runs_in_collected_dict[0]==run:
                                    run_is_in_pickle_file=True
                    if run_is_in_pickle_file:
#                        self.data_collection_dict[xtal][0].append(self.data_collection_dict_collected[xtal][0])
                        self.data_collection_dict[xtal][0]+=self.data_collection_dict_collected[xtal][0]
                        self.data_collection_dict[xtal][1]+=self.data_collection_dict_collected[xtal][1]
                        self.data_collection_dict[xtal][2]+=self.data_collection_dict_collected[xtal][2]
                        self.data_collection_dict[xtal][3]+=self.data_collection_dict_collected[xtal][3]
#                        self.data_collection_dict[xtal][4]=self.data_collection_dict_collected[xtal][4]
#                        print 'already done'
                        continue
                    timestamp=datetime.fromtimestamp(os.path.getmtime(runs)).strftime('%Y-%m-%d %H:%M:%S')
                    run_list.append([(run,timestamp,visit)])
                    self.data_collection_dict[xtal][0].append([run,timestamp,visit])
                    for (path, dirs, files) in os.walk(runs):
                        if 'edna' in dirs:
                            dirs.remove('edna')
                        if 'auto_mrbump' in dirs:
                            dirs.remove('auto_mrbump')
                        if 'fast_ep' in dirs:
                            dirs.remove('fast_ep')
                        if 'multi-xia2' in dirs:
                            dirs.remove('multi-xia2')
                        for file_name in files:
                            if file_name.endswith('aimless.log') and self.target in path:
                                logfile_list.append(path+'/'+file_name)
                                continue
                for image in glob.glob(visit_directory+'/jpegs/'+self.target+'/'+xtal+'/*'):
                    if image.endswith('t.png'):
                        image_list.append(image)
#                    if image.endswith('thumb.jpeg'):
#                        image_list.append(image)
                    if image.endswith('_.png'):
                        image_list.append(image)
                self.data_collection_dict[xtal][1]+=image_list
                self.data_collection_dict[xtal][2]+=logfile_list
                progress += progress_step
                self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

                # convert images to strong and attach to
                for image in self.data_collection_dict[xtal][1]:
                    try:
                        if image in self.data_collection_dict_collected[xtal][1]:
                            continue
                    except KeyError:
                        image_file=open(image,"rb")
                        image_string=base64.b64encode(image_file.read())
                        image_string_list.append((image[image.rfind('/')+1:],image_string))
                self.data_collection_dict[xtal][3]+=image_string_list

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
                    if sample in self.data_collection_statistics_dict_collected:
                        for entry in self.data_collection_statistics_dict_collected[sample]:
                            if logfile==entry[1]:
                                self.data_collection_statistics_dict[sample].append(entry)
                                already_parsed=True
                    if already_parsed:
#                        print 'hallo'
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
#            print sample
#            print self.data_collection_statistics_dict[sample]
            if self.data_collection_statistics_dict[sample][0]=='#':
#                print sample
                continue
            select_stage_one_list = []
            found=0
            if self.reference_file_list != []:
                for reference_file in self.reference_file_list:
                    for aimless_file in self.data_collection_statistics_dict[sample]:
                        try:
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
