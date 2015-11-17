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
        if os.path.isfile(os.path.join(self.database_directory,'data_collection_summary.pkl')):
            #data_collection_dict = pickle.load( open( os.path.join(self.database_directory,'data_collection_summary.pkl'), "rb" ) )
            summary = pickle.load( open( os.path.join(self.database_directory,'data_collection_summary.pkl'), "rb" ) )
            self.data_collection_dict_collected=summary[0]

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
                tmp=[]
                for image in self.data_collection_dict[xtal][1]:
#                    if image in self.data_collection_dict_collected[xtal][1]:
#                        continue
                    image_file=open(image,"rb")
                    image_string=base64.b64encode(image_file.read())
                    image_string_list.append((image[image.rfind('/')+1:],image_string))
                self.data_collection_dict[xtal][3]+=image_string_list

            search_cycle+=1

        if not len(self.data_collection_dict)==0:
            progress_step=100/float(len(self.data_collection_dict))
        progress=0
        for sample in sorted(self.data_collection_dict):
            self.data_collection_statistics_dict[sample]=[]
#            print sample,self.data_collection_dict[sample][2]
            if not self.data_collection_dict[sample][2]==[]:
                for index,logfile in enumerate(self.data_collection_dict[sample][2]):
                    self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'Step 2 of 3: parsing aimless logfiles -> '+sample)
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
                self.data_collection_statistics_dict[sample]+='###'*20
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
                print sample
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

        self.emit(QtCore.SIGNAL('create_widgets_for_autoprocessing_results'), [self.data_collection_dict,
                                                                            self.data_collection_statistics_dict])
