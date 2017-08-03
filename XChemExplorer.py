# last edited: 03/08/2017, 18:00

import os, sys, glob
from datetime import datetime
import time
import sqlite3

import subprocess

from PyQt4 import QtGui, QtCore, QtWebKit

import pickle
import base64
import math
import multiprocessing
import webbrowser
import getpass

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'web'))
from XChemUtils import parse
from XChemUtils import external_software
from XChemUtils import helpers
import XChemThread
import XChemDB
import XChemPANDDA
import XChemToolTips
import XChemMain
import XChemPlots
import XChemLog
import XChemProcess
import XChemDeposit
import XChemWeb

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import numpy as np

class XChemExplorer(QtGui.QApplication):
    def __init__(self,args):
        QtGui.QApplication.__init__(self,args)

        self.xce_version='v1.0-beta.4.1'

        # general settings
        self.allowed_unitcell_difference_percent=12
        self.acceptable_low_resolution_limit_for_data=3.5
        self.filename_root='${samplename}'
        self.data_source_set=False
        self.max_queue_jobs=100

        #
        # directories
        #

        self.current_directory=os.getcwd()
        self.xce_logfile=os.path.join(self.current_directory,'xce.log')
        try:
            XChemLog.startLog(self.xce_logfile).create_logfile(self.xce_version)
        except IOError:
            self.xce_logfile=os.path.join(self.current_directory,'xce_'+getpass.getuser()+'.log')
            XChemLog.startLog(self.xce_logfile).create_logfile(self.xce_version)
        self.update_log=XChemLog.updateLog(self.xce_logfile)
        self.update_log.insert('new session started')
        self.diffraction_data_directory=self.current_directory
        self.diffraction_data_search_info='n/a'
        self.diffraction_data_reference_mtz='ignore'
        self.html_export_directory=os.getcwd()

        if 'labxchem' in self.current_directory:
            self.labxchem_directory='/'+os.path.join(*self.current_directory.split('/')[1:6])    # need splat operator: *
            self.beamline_directory=os.path.join(self.labxchem_directory,'processing','beamline')
            self.initial_model_directory=os.path.join(self.labxchem_directory,'processing','analysis','initial_model')
            self.reference_directory=os.path.join(self.labxchem_directory,'processing','reference')
            self.database_directory=os.path.join(self.labxchem_directory,'processing','database')
            self.panddas_directory=os.path.join(self.labxchem_directory,'processing','analysis','panddas')
            self.data_collection_summary_file=os.path.join(self.database_directory,str(os.getcwd().split('/')[5])+'_summary.pkl')
            self.data_source_file=''
            self.html_export_directory=os.path.join(self.labxchem_directory,'processing','html')
            self.group_deposit_directory=os.path.join(self.labxchem_directory,'processing','group_deposition')
            if os.path.isfile(os.path.join(self.labxchem_directory,'processing','database','soakDBDataFile.sqlite')):
                self.data_source_file='soakDBDataFile.sqlite'
                self.database_directory=os.path.join(self.labxchem_directory,'processing','database')
                self.data_source_set=True
                self.db=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
                self.db.create_missing_columns()

            self.ccp4_scratch_directory=os.path.join(self.labxchem_directory,'processing','tmp')

            if not os.path.isdir(self.beamline_directory):
                os.mkdir(self.beamline_directory)
            if not os.path.isdir(os.path.join(self.labxchem_directory,'processing','analysis')):
                os.mkdir(os.path.join(self.labxchem_directory,'processing','analysis'))
            if not os.path.isdir(self.initial_model_directory):
                os.mkdir(self.initial_model_directory)
            if not os.path.isdir(self.panddas_directory):
                os.mkdir(self.panddas_directory)
            if not os.path.isdir(self.reference_directory):
                os.mkdir(self.reference_directory)
            if not os.path.isdir(self.database_directory):
                os.mkdir(self.database_directory)
            if not os.path.isdir(self.ccp4_scratch_directory):
                os.mkdir(self.ccp4_scratch_directory)
            if not os.path.isdir(self.html_export_directory):
                os.mkdir(self.html_export_directory)
            if not os.path.isdir(self.group_deposit_directory):
                os.mkdir(self.group_deposit_directory)

        else:
            self.beamline_directory=self.current_directory
            self.initial_model_directory=self.current_directory
            self.reference_directory=self.current_directory
            self.database_directory=self.current_directory
            self.data_source_file=''
            self.ccp4_scratch_directory=os.getenv('CCP4_SCR')
            self.panddas_directory=self.current_directory
            self.data_collection_summary_file=''
            self.group_deposit_directory=self.current_directory

        #
        # Preferences
        #

        self.preferences_data_to_copy = [
            ['aimless logiles and merged mtz only',                             'mtz_log_only'],
            #['All Files in the respective auto-processing directory',           'everything'],
                    ]

        self.preferences_selection_mechanism = [    'IsigI*Comp*UniqueRefl',
                                                    'highest_resolution',
                                                    'lowest_Rfree'              ]

        self.preferences =  {   'processed_data_to_copy':       'mtz_log_only',
                                'dataset_selection_mechanism':  'IsigI*Comp*UniqueRefl' }


        #
        # Settings
        #

        self.settings =     {'current_directory':               self.current_directory,
                             'beamline_directory':              self.beamline_directory,
                             'data_collection_summary':         self.data_collection_summary_file,
                             'initial_model_directory':         self.initial_model_directory,
                             'panddas_directory':               self.panddas_directory,
                             'reference_directory':             self.reference_directory,
                             'database_directory':              self.database_directory,
                             'data_source':                     os.path.join(self.database_directory,self.data_source_file),
                             'ccp4_scratch':                    self.ccp4_scratch_directory,
                             'unitcell_difference':             self.allowed_unitcell_difference_percent,
                             'too_low_resolution_data':         self.acceptable_low_resolution_limit_for_data,
                             'filename_root':                   self.filename_root,
                             'preferences':                     self.preferences,
                             'xce_logfile':                     self.xce_logfile,
                             'max_queue_jobs':                  self.max_queue_jobs,
                             'diffraction_data_directory':      self.diffraction_data_directory,
                             'html_export_directory':           self.html_export_directory,
                             'group_deposit_directory':         self.group_deposit_directory,
                             'remote_qsub':                     ''  }


        #
        # Deposition
        #

        self.deposit_dict = {}


        #
        # internal lists and dictionaries
        #

        self.data_collection_list=[]
        self.visit_list=[]
        self.target=''
        self.dataset_outcome_combobox_dict={}
        self.data_collection_dict={}
        self.xtal_db_dict={}
        self.pandda_analyse_input_table_dict={}
        self.dewar_configuration_dict={}
        self.data_collection_statistics_dict={}
        self.initial_model_dimple_dict={}       # contains toggle button if dimple should be run
        self.reference_file_list=[]
        self.all_columns_in_data_source=XChemDB.data_source(os.path.join(self.database_directory,
                                                                         self.data_source_file)).return_column_list()
        self.albula_button_dict={}              # using dials.image_viewer instead of albula, but keep name for dictionary
        self.xtalform_dict={}

        self.dataset_outcome_dict={}            # contains the dataset outcome buttons
        self.data_collection_table_dict={}      # contains the dataset table
        self.data_collection_image_dict={}
        self.data_collection_column_three_dict={}
        self.data_collection_summary_dict={}
        self.diffraction_data_table_dict={}
        self.summary_table_dict={}
        self.main_data_collection_table_exists=False
        self.timer_to_check_for_new_data_collection = QtCore.QTimer()
#        self.timer_to_check_for_new_data_collection.timeout.connect(self.check_for_new_autoprocessing_or_rescore(False))

        self.target_list,self.visit_list=XChemMain.get_target_and_visit_list(self.beamline_directory)
#        self.target_list,self.visit_list=XChemMain.get_target_and_visit_list_for_Pietro(self.beamline_directory)

        self.diffraction_data_dict={}

        #
        # internal switches and flags
        #

        self.explorer_active=0
        self.coot_running=0
        self.progress_bar_start=0
        self.progress_bar_step=0
        self.albula = None
        self.albula_subframes=[]
        self.show_diffraction_image = None
        self.data_collection_details_currently_on_display=None      # can be any widget to be displayed in data collection summary tab

        self.dataset_outcome = [    "success",
                                    "Failed - centring failed",
                                    "Failed - no diffraction",
                                    "Failed - processing",
                                    "Failed - loop empty",
                                    "Failed - loop broken",
                                    "Failed - low resolution",
                                    "Failed - no X-rays",
                                    "Failed - unknown"  ]

        self.refinement_stage = [       '0 - All Datasets',
                                        '1 - Analysis Pending',
                                        '2 - PANDDA model',
                                        '3 - In Refinement',
                                        '4 - CompChem ready',
                                        '5 - Deposition ready',
                                        '6 - Deposited'         ]

        #
        # checking for external software packages
        #

        self.external_software=external_software(self.xce_logfile).check()
        if self.external_software['acedrg']:
            self.restraints_program='acedrg'
            self.update_log.insert('will use ACEDRG for generation of ligand coordinates and restraints')
        elif self.external_software['phenix.elbow']:
            self.restraints_program='phenix.elbow'
            self.update_log.insert('will use PHENIX.ELBOW for generation of ligand coordinates and restraints')
        elif self.external_software['grade']:
            self.restraints_program='grade'
            self.update_log.insert('will use GRADE for generation of ligand coordinates and restraints')
        else:
            self.restraints_program=''
            self.update_log.insert('No program for generation of ligand coordinates and restraints available!')

        self.using_remote_qsub_submission=False
        self.remote_qsub_submission="ssh <dls fed ID>@nx.diamond.ac.uk 'module load global/cluster; qsub'"

        # start GUI
        self.start_GUI()
        self.exec_()

    def start_GUI(self):

        # GUI setup
        self.window=QtGui.QWidget()
        self.window.setWindowTitle("XChemExplorer")

        #size_policy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        #self.window.setSizePolicy(size_policy)
        self.screen = QtGui.QDesktopWidget().screenGeometry()
        #self.window.setFixedSize(self.screen.width(),self.screen.height()-70)


        ######################################################################################
        # Menu Widget
        menu_bar = QtGui.QMenuBar()

        file = menu_bar.addMenu("&File")
        load=QtGui.QAction("Open Config File", self.window)
        load.setShortcut('Ctrl+O')
        load.triggered.connect(self.open_config_file)
        save=QtGui.QAction("Save Config File", self.window)
        save.setShortcut('Ctrl+S')
        save.triggered.connect(self.save_config_file)
        quit=QtGui.QAction("Quit", self.window)
        quit.setShortcut('Ctrl+Q')
        quit.triggered.connect(self.quit_xce)
        file.addAction(load)
        file.addAction(save)
        file.addAction(quit)

        datasource_menu = menu_bar.addMenu("&Data Source")
        reload_samples_from_datasource=QtGui.QAction('Reload Samples from Datasource',self.window)
        reload_samples_from_datasource.triggered.connect(self.datasource_menu_reload_samples)
        save_samples_to_datasource=QtGui.QAction('Save Samples to Datasource',self.window)
        save_samples_to_datasource.triggered.connect(self.datasource_menu_save_samples)
        import_csv_file_into_datasource=QtGui.QAction('Import CSV file into Datasource',self.window)
        import_csv_file_into_datasource.triggered.connect(self.datasource_menu_import_csv_file)
        export_csv_file_into_datasource=QtGui.QAction('Export CSV file from Datasource',self.window)
        export_csv_file_into_datasource.triggered.connect(self.datasource_menu_export_csv_file)
        update_datasource=QtGui.QAction('Update Datasource from file system',self.window)
        update_datasource.triggered.connect(self.datasource_menu_update_datasource)
        select_columns_to_show=QtGui.QAction('Select columns to show',self.window)
        select_columns_to_show.triggered.connect(self.select_datasource_columns_to_display)
        create_new_data_source=QtGui.QAction('Create New Data Source (SQLite)',self.window)
        create_new_data_source.triggered.connect(self.create_new_data_source)
        export_csv_for_WONKA=QtGui.QAction('export CSV for WONKA',self.window)
        export_csv_for_WONKA.triggered.connect(self.export_data_for_WONKA)

        datasource_menu.addAction(reload_samples_from_datasource)
        datasource_menu.addAction(save_samples_to_datasource)
        datasource_menu.addAction(import_csv_file_into_datasource)
        datasource_menu.addAction(export_csv_file_into_datasource)
        datasource_menu.addAction(update_datasource)
        datasource_menu.addAction(select_columns_to_show)
        datasource_menu.addAction(create_new_data_source)
        datasource_menu.addAction(export_csv_for_WONKA)

        preferences_menu = menu_bar.addMenu("&Preferences")
        show_preferences=QtGui.QAction('Edit Preferences',self.window)
        show_preferences.triggered.connect(self.show_preferences)
        preferences_menu.addAction(show_preferences)

        deposition_menu = menu_bar.addMenu("&Deposition")
        edit_deposition_info=QtGui.QAction('Edit Information',self.window)
        edit_deposition_info.triggered.connect(self.deposition_data)
        deposition_menu.addAction(edit_deposition_info)
        export_results_to_html=QtGui.QAction('Export to HTML',self.window)
        export_results_to_html.triggered.connect(self.export_to_html)
        deposition_menu.addAction(export_results_to_html)

        find_apo_structures=QtGui.QAction('find PanDDA apo structures',self.window)
        find_apo_structures.triggered.connect(self.create_missing_apo_records_in_depositTable)
        deposition_menu.addAction(find_apo_structures)

        update_file_information_of_apo_records=QtGui.QAction('update file info of apo structures',self.window)
        update_file_information_of_apo_records.triggered.connect(self.update_file_information_of_apo_records)
        deposition_menu.addAction(update_file_information_of_apo_records)

        self.prepare_mmcif_files_dict={}

        prepare_mmcif_files_for_apo_structures=QtGui.QAction('prepare mmcif for apo structures',self.window)
        prepare_mmcif_files_for_apo_structures.triggered.connect(self.prepare_models_for_deposition)
        deposition_menu.addAction(prepare_mmcif_files_for_apo_structures)
        self.prepare_mmcif_files_dict['apo']=prepare_mmcif_files_for_apo_structures

        prepare_mmcif_files_for_ligand_bound_structures=QtGui.QAction('prepare mmcif for ligand bound structures',self.window)
        prepare_mmcif_files_for_ligand_bound_structures.triggered.connect(self.prepare_models_for_deposition)
        deposition_menu.addAction(prepare_mmcif_files_for_ligand_bound_structures)
        self.prepare_mmcif_files_dict['ligand_bound']=prepare_mmcif_files_for_ligand_bound_structures

        prepare_for_group_deposition_upload=QtGui.QAction('copy files to group deposition directory',self.window)
        prepare_for_group_deposition_upload.triggered.connect(self.prepare_for_group_deposition_upload)
        deposition_menu.addAction(prepare_for_group_deposition_upload)

        enter_pdb_codes=QtGui.QAction('Update DB with PDB codes',self.window)
        enter_pdb_codes.triggered.connect(self.enter_pdb_codes)
        deposition_menu.addAction(enter_pdb_codes)

        check_smiles=QtGui.QAction('Check SMILES',self.window)
        check_smiles.triggered.connect(self.check_smiles_in_db_and_pdb)
        deposition_menu.addAction(check_smiles)

        ### RACHAEL'S PROASIS STUFF ###

        self.proasis_directory = '/dls/science/groups/proasis/'

        # function for adding a new project
        def create_project(name):
            # make relevant project directory in proasis LabXChem folder
            os.system(str('mkdir ' + os.path.join(self.proasis_directory, 'LabXChem', name)))
            perm_string = str('chmod u=rw,g=rx,o=rx ' + os.path.join(self.proasis_directory, 'LabXChem', name))
            os.system(perm_string)
            # make reference file directory in project directory
            os.system(str('mkdir ' + os.path.join(self.proasis_directory, 'LabXChem', name, 'reference')))
            perm_string = str('chmod u=rw,g=rx,o=rx ' + os.path.join(self.proasis_directory, 'LabXChem',
                                                                     name, 'reference'))
            os.system(perm_string)
            # create a temporary job to add the project in proasis schedule
            temp_job = open(os.path.join(self.proasis_directory, 'Scripts/scheduled_jobs/temp_jobs',
                                         str(name + '.sh')), 'w')
            perm_string = str(
                'chmod 770 ' + os.path.join(self.proasis_directory, 'Scripts/scheduled_jobs/temp_jobs', str(name
                                                                                                            + '.sh')))
            os.system(perm_string)
            job_string = str('/usr/local/Proasis2/utils/addnewproject.py -q OtherClasses -p ' + name)
            temp_job.write(str(job_string))
            temp_job.close()

        # add proasis menu to main menu
        self.proasis_menu = menu_bar.addMenu('Proasis')
        # connect to soakDB to get proasis info
        if os.path.isfile(os.path.join(self.database_directory, self.data_source_file)):
            conn = sqlite3.connect(os.path.join(self.database_directory, self.data_source_file))
            c = conn.cursor()

        # Project details or add project in menu
        counter = 0
        try:
            # get protein name from soakDB - this will be the proasis project name
            for row in c.execute('SELECT Protein FROM soakDB;'):
                counter+=1
                # if there is only one protein name in soakDB - all is good - happy days
                if counter==1:
                    self.proasis_name = str(row[0])
                # otherwise - give a warning
                # TODO: If this is actually ever encountered - deal with it. Should be fine.
                if counter > 1:
                    print('WARNING: More than one protein name found (proasis)')
                # If the project directory already exists in proasis dir, project should exist in proasis
                if os.path.isdir(os.path.join('/dls/science/groups/proasis/LabXChem/', self.proasis_name)):
                    # show project name in menu (no action when clicked)
                    self.proasis_project = QtGui.QAction(str('Project Name: ' + self.proasis_name), self.window)
                    self.proasis_menu.addAction(self.proasis_project)
        # should catch if project doesnt exitst
        except AttributeError:
            self.update_log.insert('cannot find %s' %os.path.join(self.database_directory, self.data_source_file))
        except UnboundLocalError:
            self.update_log.insert('cannot find %s' %os.path.join(self.database_directory, self.data_source_file))
        except:
                # option to create project, action = create_project()
                self.proasis_project = QtGui.QAction(str('Create Project for ' + self.proasis_name + '...'),
                                                     self.window)
                self.proasis_project.triggered.connect(lambda: create_project(self.proasis_name))
                self.proasis_menu.addAction(self.proasis_project)


        # Lead details or add lead in menu
        counter = 0
        try:
            # check if there is a lead (from soakDB)
            for row in c.execute('SELECT proasisID from proasisLead'):
                counter+=1
                if counter==1:
                    # If so, display id of lead in menu, no action if clicked
                    self.proasis_lead = QtGui.QAction(str('Lead ID: ' + str(row[0])), self.window)
                    self.proasis_menu.addAction(self.proasis_lead)
                # otherwise, if you can find the pandda_analyse_sites.csv file, allow lead to be added
                elif os.path.isfile(os.path.join(self.panddas_directory, 'analyses/pandda_analyse_sites.csv')):
                    self.proasis_lead = QtGui.QAction(str('Create lead from pandda sites...'), self.window)
                    self.proasis_menu.addAction(self.proasis_lead)
        # If no lead or sites file, error message. No action on click
        except:
                self.proasis_lead = QtGui.QAction(str('Site info not found... '
                                                      'please run pandda analyse before adding lead'), self.window)
                self.proasis_lead.triggered.connect(lambda:self.add_lead())
                self.proasis_menu.addAction(self.proasis_lead)

        # Hit details or add hits (refined) in menu
        counter = 0
        try:
            # count the number of hits in proasis if they exist (from soakDB)
            for row in c.execute('SELECT proasisID from proasis'):
                counter += 1
            no_hits = counter
            # display no of hits (proasis) in menu, no action if clicked
            self.proasis_hits = QtGui.QAction(str('Hits in proasis: ' + str(no_hits)), self.window)
            self.proasis_menu.addAction(self.proasis_hits)
        # otherwise, try to add hits to proasis (if there are no hits, the job will still run and hits will be added as
        # they are refined - i.e. when refine.bound.pdb file is detected for a refinement detailed in soakDB)
        except:
            self.proasis_hits = QtGui.QAction(str('Attempt to add refined hits to proasis...'), self.window)
            self.proasis_hits.triggered.connect(lambda:self.add_hits())
            self.proasis_menu.addAction(self.proasis_hits)

        ##############################

        def openFile(file):
   		if sys.platform == 'linux2':
        		subprocess.call(["xdg-open", file])
    		else:
        		os.startfile(file)

        help_menu = menu_bar.addMenu("&Help")
        load_xce_tutorial = QtGui.QAction('Open XCE tutorial', self.window)
        file = '/dls/science/groups/i04-1/software/docs/XChemExplorer.pdf'
        load_xce_tutorial.triggered.connect(lambda:openFile(file))
        help_menu.addAction(load_xce_tutorial)

        load_xce_troubleshoot = QtGui.QAction('Troubleshooting', self.window)
        file2 = '/dls/science/groups/i04-1/software/xce_troubleshooting.pdf'
        load_xce_troubleshoot.triggered.connect(lambda:openFile(file2))
        help_menu.addAction(load_xce_troubleshoot)

        ######################################################################################:
        #
        # Workflow @ Task Containers
        #

        self.workflow =     [   'Overview',         # 0
                                'Datasets',         # 1
                                'Maps',             # 2
                                'PANDDAs',          # 3
                                'Refinement',       # 4
                                'Deposition',       # 6
                                'Settings'   ]      # 5

        self.workflow_dict = {  self.workflow[0]:       'Overview',
                                self.workflow[1]:       'Datasets',
                                self.workflow[2]:       'Maps',
                                self.workflow[3]:       'PANDDAs',
                                self.workflow[4]:       'Refinement',
                                self.workflow[6]:       'Settings',
                                self.workflow[5]:       'Deposition'        }

        self.workflow_widget_dict = {}

        #
        # check http://doc.qt.io/qt-4.8/stylesheet-customizing.html#the-box-model
        #

        headlineLabelfont = QtGui.QFont("Arial", 20, QtGui.QFont.Bold)

        #
        # @ Update from datasource button ###################################################
        #

        update_from_datasource_button=QtGui.QPushButton("Update Tables\nFrom Datasource")
        update_from_datasource_button.setToolTip(XChemToolTips.update_from_datasource_button_tip())
        update_from_datasource_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px; background: rgb(140,140,140) }")
        update_from_datasource_button.setFont(headlineLabelfont)
        update_from_datasource_button.clicked.connect(self.datasource_menu_reload_samples)

        #
        # @ Datasets ########################################################################
        #

        self.dataset_tasks = [  'Get New Results from Autoprocessing',
                               # 'Save Files from Autoprocessing to Project Folder',
                                'Run DIMPLE on All Autoprocessing MTZ files',
                                'Rescore Datasets',
                                'Read PKL file',
                                'Run xia2 on selected datasets',
                                'Run xia2 on selected datasets - overwrite' ]

        frame_dataset_task=QtGui.QFrame()
        frame_dataset_task.setFrameShape(QtGui.QFrame.StyledPanel)
        frame_dataset_task.setStyleSheet("QFrame { border: 1px solid black; border-radius: 1px; padding: 0px; margin: 0px }")
        vboxTask=QtGui.QVBoxLayout()
        label=QtGui.QLabel('Datasets')
        label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        label.setFont(headlineLabelfont)
        label.setStyleSheet(" QLabel { border: 1px solid black; border-radius: 1px; background: rgb(240,255,140); padding: 0px; margin: 0px }")
        vboxTask.addWidget(label)
        hboxAction=QtGui.QHBoxLayout()
        self.dataset_tasks_combobox = QtGui.QComboBox()
        for task in self.dataset_tasks:
            self.dataset_tasks_combobox.addItem(task)
        self.dataset_tasks_combobox.setToolTip(XChemToolTips.dataset_task_tip())
        self.dataset_tasks_combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")
        hboxAction.addWidget(self.dataset_tasks_combobox)
        vboxButton=QtGui.QVBoxLayout()
        dataset_task_run_button=QtGui.QPushButton("Run")
        dataset_task_run_button.setToolTip(XChemToolTips.dataset_task_run_button_tip())
        dataset_task_run_button.clicked.connect(self.button_clicked)
        dataset_task_run_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
        vboxButton.addWidget(dataset_task_run_button)
        dataset_task_status_button=QtGui.QPushButton("Status")
        dataset_task_status_button.setToolTip(XChemToolTips.dataset_task_status_button_tip())
        dataset_task_status_button.clicked.connect(self.button_clicked)
        dataset_task_status_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
        vboxButton.addWidget(dataset_task_status_button)
        hboxAction.addLayout(vboxButton)
        vboxTask.addLayout(hboxAction)
        vboxTask.setSpacing(0)
        vboxTask.setMargin(0)
        frame_dataset_task.setLayout(vboxTask)

        self.workflow_widget_dict['Datasets']=[self.dataset_tasks_combobox,dataset_task_run_button,dataset_task_status_button]


        #
        # @ MAP & CIF files #######################################################################
        #

        self.map_cif_file_tasks = [ 'Run DIMPLE on selected MTZ files',
                                    'Remove selected DIMPLE PDB/MTZ files',
                                    'Create CIF/PDB/PNG file of ALL compounds',
                                    'Create CIF/PDB/PNG file of NEW compounds',
                                    'Create CIF/PDB/PNG file of SELECTED compounds' ]

        frame_map_cif_file_task=QtGui.QFrame()
        frame_map_cif_file_task.setFrameShape(QtGui.QFrame.StyledPanel)
        frame_map_cif_file_task.setStyleSheet("QFrame { border: 1px solid black; border-radius: 1px; padding: 0px; margin: 0px }")
        vboxTask=QtGui.QVBoxLayout()
        label=QtGui.QLabel('Maps & Restraints')
        label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        label.setFont(headlineLabelfont)
        label.setStyleSheet(" QLabel { border: 1px solid black; border-radius: 1px; background: rgb(140,255,150); padding: 0px; margin: 0px }")
        vboxTask.addWidget(label)
        hboxAction=QtGui.QHBoxLayout()
        self.map_cif_file_tasks_combobox = QtGui.QComboBox()
        for task in self.map_cif_file_tasks:
            self.map_cif_file_tasks_combobox.addItem(task)
        self.map_cif_file_tasks_combobox.setToolTip(XChemToolTips.map_cif_file_task_tip())
        self.map_cif_file_tasks_combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")
        hboxAction.addWidget(self.map_cif_file_tasks_combobox)
        vboxButton=QtGui.QVBoxLayout()
        map_cif_file_task_run_button=QtGui.QPushButton("Run")
        map_cif_file_task_run_button.setToolTip(XChemToolTips.map_cif_file_task_run_button_tip())
        map_cif_file_task_run_button.clicked.connect(self.button_clicked)
        map_cif_file_task_run_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
        vboxButton.addWidget(map_cif_file_task_run_button)
        map_cif_file_task_status_button=QtGui.QPushButton("Status")
        map_cif_file_task_status_button.setToolTip(XChemToolTips.map_cif_file_task_status_button_tip())
        map_cif_file_task_status_button.clicked.connect(self.button_clicked)
        map_cif_file_task_status_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
        vboxButton.addWidget(map_cif_file_task_status_button)
        hboxAction.addLayout(vboxButton)
        vboxTask.addLayout(hboxAction)
        vboxTask.setSpacing(0)
        vboxTask.setMargin(0)
        frame_map_cif_file_task.setLayout(vboxTask)

        self.workflow_widget_dict['Maps']=[self.map_cif_file_tasks_combobox,map_cif_file_task_run_button,map_cif_file_task_status_button]

        #####################################################################################

        #
        # @ PANDDAs #########################################################################
        #

        self.panddas_file_tasks = [ 'pandda.analyse',
                                    'pandda.inspect',
                                    'run pandda.inspect at home',
                                    'Export NEW PANDDA models',
                                    'Export ALL PANDDA models',
                                    'Show HTML summary',
                                    'Update datasource with results from pandda.inspect',
                                    'cluster datasets',
                                    'Event Map -> SF',
                                    'check modelled ligands',
                                    'pre-run for ground state model',
                                    'Build ground state model' ]

        frame_panddas_file_task=QtGui.QFrame()
        frame_panddas_file_task.setFrameShape(QtGui.QFrame.StyledPanel)
        frame_panddas_file_task.setStyleSheet("QFrame { border: 1px solid black; border-radius: 1px; padding: 0px; margin: 0px }")
        vboxTask=QtGui.QVBoxLayout()
        label=QtGui.QLabel('Hit Identification')
        label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        label.setFont(headlineLabelfont)
        label.setStyleSheet(" QLabel { border: 1px solid black; border-radius: 1px; background: rgb(140,200,255); padding: 0px; margin: 0px }")
        vboxTask.addWidget(label)
        hboxAction=QtGui.QHBoxLayout()
        self.panddas_file_tasks_combobox = QtGui.QComboBox()
        for task in self.panddas_file_tasks:
            self.panddas_file_tasks_combobox.addItem(task)
        self.panddas_file_tasks_combobox.setToolTip(XChemToolTips.panddas_file_task_tip())
        self.panddas_file_tasks_combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")
        hboxAction.addWidget(self.panddas_file_tasks_combobox)
        vboxButton=QtGui.QVBoxLayout()
        panddas_file_task_run_button=QtGui.QPushButton("Run")
        panddas_file_task_run_button.setToolTip(XChemToolTips.panddas_file_task_run_button_tip())
        panddas_file_task_run_button.clicked.connect(self.button_clicked)
        panddas_file_task_run_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
        vboxButton.addWidget(panddas_file_task_run_button)
        panddas_file_task_status_button=QtGui.QPushButton("Status")
        panddas_file_task_status_button.setToolTip(XChemToolTips.panddas_file_task_status_button_tip())
        panddas_file_task_status_button.clicked.connect(self.button_clicked)
        panddas_file_task_status_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
        vboxButton.addWidget(panddas_file_task_status_button)
        hboxAction.addLayout(vboxButton)
        vboxTask.addLayout(hboxAction)
        vboxTask.setSpacing(0)
        vboxTask.setMargin(0)
        frame_panddas_file_task.setLayout(vboxTask)

        self.workflow_widget_dict['PANDDAs']=[self.panddas_file_tasks_combobox,panddas_file_task_run_button,panddas_file_task_status_button]

        #####################################################################################

        #
        # @ Refine ##########################################################################
        #

        self.refine_file_tasks = [ 'Open COOT',
                                   'Open COOT - new interface',
                                   'Open COOT for old PanDDA',
                                   'Update Deposition Table',
                                   'Prepare Group Deposition'   ]

        frame_refine_file_task=QtGui.QFrame()
        frame_refine_file_task.setFrameShape(QtGui.QFrame.StyledPanel)
        frame_refine_file_task.setStyleSheet("QFrame { border: 1px solid black; border-radius: 1px; padding: 0px; margin: 0px }")
        vboxTask=QtGui.QVBoxLayout()
        label=QtGui.QLabel('Refinement')
        label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        label.setFont(headlineLabelfont)
        label.setStyleSheet(" QLabel { border: 1px solid black; border-radius: 1px; background: rgb(245,190,255); padding: 0px; margin: 0px }")
        vboxTask.addWidget(label)
        hboxAction=QtGui.QHBoxLayout()
        self.refine_file_tasks_combobox = QtGui.QComboBox()
        for task in self.refine_file_tasks:
            self.refine_file_tasks_combobox.addItem(task)
        self.refine_file_tasks_combobox.setToolTip(XChemToolTips.refine_file_task_tip())
        self.refine_file_tasks_combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")
        hboxAction.addWidget(self.refine_file_tasks_combobox)
        vboxButton=QtGui.QVBoxLayout()
        refine_file_task_run_button=QtGui.QPushButton("Run")
        refine_file_task_run_button.setToolTip(XChemToolTips.refine_file_task_run_button_tip())
        refine_file_task_run_button.clicked.connect(self.button_clicked)
        refine_file_task_run_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
        vboxButton.addWidget(refine_file_task_run_button)
        refine_file_task_status_button=QtGui.QPushButton("Status")
        refine_file_task_status_button.setToolTip(XChemToolTips.refine_file_task_status_button_tip())
        refine_file_task_status_button.clicked.connect(self.button_clicked)
        refine_file_task_status_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
        vboxButton.addWidget(refine_file_task_status_button)
        hboxAction.addLayout(vboxButton)
        vboxTask.addLayout(hboxAction)
        vboxTask.setSpacing(0)
        vboxTask.setMargin(0)
        frame_refine_file_task.setLayout(vboxTask)

        self.workflow_widget_dict['Refinement']=[self.refine_file_tasks_combobox,refine_file_task_run_button,refine_file_task_status_button]

        #####################################################################################

        ######################################################################################
        #
        # Workflow @ Tabs for Tasks
        #

        #
        # @ tab widget #######################################################################
        #

        self.main_tab_widget = QtGui.QTabWidget()
        #self.main_tab_widget.setSizePolicy(size_policy)
        self.tab_dict={}
        for page in self.workflow:
            tab=QtGui.QWidget()
            vbox=QtGui.QVBoxLayout(tab)
            self.main_tab_widget.addTab(tab,page)
            self.tab_dict[page]=[tab,vbox]

        #
        # @ Data Source Tab ###################################################################
        #

        overview_tab_widget = QtGui.QTabWidget()
        #overview_tab_widget.setSizePolicy(size_policy)
        self.tab_dict[self.workflow_dict['Overview']][1].addWidget(overview_tab_widget)
        overview_tab_list = [   'Data Source',
                                'Summary'    ]

        self.overview_tab_dict={}
        for page in overview_tab_list:
            tab=QtGui.QWidget()
            vbox=QtGui.QVBoxLayout(tab)
            overview_tab_widget.addTab(tab,page)
            self.overview_tab_dict[page]=[tab,vbox]

        self.data_source_columns_to_display=[   'Sample ID',
                                                'Compound ID',
                                                'Smiles',
                                                'Visit',
                                                'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                                'Refinement\nRfree',
                                                'Data Collection\nDate',
                                                'Puck',
                                                'PuckPosition',
                                                'Ligand\nConfidence'    ]

        self.mounted_crystal_table=QtGui.QTableWidget()
        self.mounted_crystal_table.setSortingEnabled(True)
        self.mounted_crystal_table.resizeColumnsToContents()
        self.overview_tab_dict['Data Source'][1].addWidget(self.mounted_crystal_table)

        # - Overview Graph ####################################################################

        self.overview_figure, self.overview_axes = plt.subplots()
        self.overview_canvas = FigureCanvas(self.overview_figure)
        self.update_summary_plot()
        self.overview_tab_dict['Summary'][1].addWidget(self.overview_canvas)



        ######################################################################################


        #
        # @ Dataset @ Data Collection Tab ####################################################
        #

        self.dls_data_collection_vbox=QtGui.QVBoxLayout()
        self.tab_dict[self.workflow_dict['Datasets']][1].addLayout(self.dls_data_collection_vbox)

        hbox=QtGui.QHBoxLayout()
        check_for_new_data_collection = QtGui.QCheckBox('Check for new data collection every two minutes')
        check_for_new_data_collection.toggle()
        check_for_new_data_collection.setChecked(False)
        check_for_new_data_collection.stateChanged.connect(self.continously_check_for_new_data_collection)
        hbox.addWidget(check_for_new_data_collection)
        hbox.addWidget(QtGui.QLabel('                                             '))
        hbox.addWidget(QtGui.QLabel('Select Target: '))
        self.target_selection_combobox = QtGui.QComboBox()
        self.populate_target_selection_combobox(self.target_selection_combobox)
        self.target_selection_combobox.activated[str].connect(self.target_selection_combobox_activated)
        hbox.addWidget(self.target_selection_combobox)
        self.target=str(self.target_selection_combobox.currentText())

        self.dls_data_collection_vbox.addLayout(hbox)


        dls_tab_widget = QtGui.QTabWidget()
        dls_tab_list = [ 'Summary',
        #                 'Dewar',
                         'Reprocess'    ]

        self.dls_tab_dict={}
        for page in dls_tab_list:
            tab=QtGui.QWidget()
            vbox=QtGui.QVBoxLayout(tab)
            dls_tab_widget.addTab(tab,page)
            self.dls_tab_dict[page]=[tab,vbox]

        # - Summary Sub-Tab ##################################################################

        data_collection_summary_list=[]
        self.data_collection_summary_column_name=[      'Sample ID',
                                                        #'Date',
                                                        'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                                        'DataProcessing\nSpaceGroup',
                                                        'DataProcessing\nRfree',
                                                        'SoakDB\nBarcode',
                                                        'GDA\nBarcode',
                                                        'Rmerge\nLow',
                                                        'auto-assigned',
                                                        'DataCollection\nOutcome',
                                                        'img1',
                                                        'img2',
                                                        'img3',
                                                        'img4',
                                                        'img5',
                                                        'Show\nDetails',
                                                        'Show Diffraction\nImage'
                                                        ]

        self.data_collection_summary_table=QtGui.QTableWidget()
        self.data_collection_summary_table.setColumnCount(len(self.data_collection_summary_column_name))
        self.data_collection_summary_table.setSortingEnabled(True)
        self.data_collection_summary_table.setHorizontalHeaderLabels(self.data_collection_summary_column_name)

        # table
        self.data_collection_summarys_vbox_for_table=QtGui.QVBoxLayout()
        self.dls_tab_dict['Summary'][1].addLayout(self.data_collection_summarys_vbox_for_table)

        # another vbox for details to be shown
        self.data_collection_summarys_vbox_for_details=QtGui.QVBoxLayout()
        self.data_collection_details_currently_on_display=None
        self.dls_tab_dict['Summary'][1].addLayout(self.data_collection_summarys_vbox_for_details)

        self.data_collection_summarys_vbox_for_table.addWidget(self.data_collection_summary_table)

        self.dls_data_collection_vbox.addWidget(dls_tab_widget)

        # - Dewar Sub-Tab ####################################################################

        self.dewar_configuration_dict={}
        self.dewar_sample_configuration_dict={}
        self.dewar_label_active=''
        self.dewar_configuration_layout = QtGui.QGridLayout()

        # create context menu
        self.popMenu = QtGui.QMenu()
        recollect=QtGui.QAction("recollect", self.window)
        recollect.triggered.connect(self.flag_sample_for_recollection)
        undo_recollect=QtGui.QAction("undo", self.window)
        undo_recollect.triggered.connect(self.undo_flag_sample_for_recollection)
        self.popMenu.addAction(recollect)
        self.popMenu.addAction(undo_recollect)

        for puck in range(38):
            for position in range(17):
                frame=QtGui.QFrame()
                frame.setFrameShape(QtGui.QFrame.StyledPanel)
                vbox_for_frame=QtGui.QVBoxLayout()
                if puck==0 and position == 0:
                    label=QtGui.QLabel('')
                    vbox_for_frame.addWidget(label)
                    frame.setLayout(vbox_for_frame)
                elif puck==0 and position != 0:
                    label=QtGui.QLabel(str(position))
                    vbox_for_frame.addWidget(label)
                    frame.setLayout(vbox_for_frame)
                elif position==0 and puck != 0:
                    label=QtGui.QLabel(str(puck))
                    vbox_for_frame.addWidget(label)
                    frame.setLayout(vbox_for_frame)
                else:
                    frame=QtGui.QPushButton('')
                    frame.setStyleSheet("font-size:5px;border-width: 0px")
                    frame.clicked.connect(self.show_html_summary_in_firefox)
                    # how to right click on button
                    self.dewar_configuration_dict[str(puck)+'-'+str(position)]=frame
                    self.dewar_sample_configuration_dict[str(puck)+'-'+str(position)]=[]
                    frame.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
                    frame.customContextMenuRequested.connect(self.on_context_menu)

                self.dewar_configuration_layout.addWidget(frame, position, puck)

#        self.dls_tab_dict['Dewar'][1].addLayout(self.dewar_configuration_layout)


        # - Reprocessing Sub-Tab ####################################################################
        reprocess_vbox=QtGui.QVBoxLayout()

        frame=QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)
        hbox=QtGui.QHBoxLayout()


        frame_select=QtGui.QFrame()
        frame_select.setFrameShape(QtGui.QFrame.StyledPanel)
        hbox_select=QtGui.QHBoxLayout()
        label=QtGui.QLabel('Data collection directory:')
        hbox_select.addWidget(label)
        dir_frame=QtGui.QFrame()
        dir_frame.setFrameShape(QtGui.QFrame.StyledPanel)
        dir_label_box=QtGui.QVBoxLayout()
        self.diffraction_data_dir_label=QtGui.QLabel(self.diffraction_data_directory)
        dir_label_box.addWidget(self.diffraction_data_dir_label)
        dir_frame.setLayout(dir_label_box)
        hbox_select.addWidget(dir_frame)
        button=QtGui.QPushButton("Select")
        button.clicked.connect(self.select_diffraction_data_directory)
        hbox_select.addWidget(button)
        frame_select.setLayout(hbox_select)
        hbox.addWidget(frame_select)
        frame_search=QtGui.QFrame()
        frame_search.setFrameShape(QtGui.QFrame.StyledPanel)
        hbox_search=QtGui.QHBoxLayout()
        button=QtGui.QPushButton("Search Datasets")
        button.clicked.connect(self.search_for_datasets)
        hbox_search.addWidget(button)
        frame_search_info=QtGui.QFrame()
        frame_search_info.setFrameShape(QtGui.QFrame.StyledPanel)
        hbox_search_info=QtGui.QHBoxLayout()
        self.diffraction_data_search_label=QtGui.QLabel(self.diffraction_data_search_info)
        hbox_search_info.addWidget(self.diffraction_data_search_label)
        frame_search_info.setLayout(hbox_search_info)
        hbox_search.addWidget(frame_search_info)
        frame_search.setLayout(hbox_search)
        hbox.addWidget(frame_search)

        frame_translate=QtGui.QFrame()
        frame_translate.setFrameShape(QtGui.QFrame.StyledPanel)
        vbox_translate=QtGui.QVBoxLayout()
        label=QtGui.QLabel('translate:\ndatasetID -> sampleID')
        label.setAlignment(QtCore.Qt.AlignCenter)
        vbox_translate.addWidget(label)
        button=QtGui.QPushButton('Open CSV')
        button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
        button.clicked.connect(self.translate_datasetID_to_sampleID)
        vbox_translate.addWidget(button)
        frame_translate.setLayout(vbox_translate)
        hbox.addWidget(frame_translate)

        hbox.addStretch(0)

        frame.setLayout(hbox)
        reprocess_vbox.addWidget(frame)


        self.reprocess_datasets_column_list=[   'Dataset ID',
                                                'Sample ID',
                                                'Run\nxia2',
                                                'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                                'Rmerge\nLow',
                                                'Dimple\nRfree',
                                                'DataProcessing\nSpaceGroup',
                                                'DataProcessing\nUnitCell',
                                                'DataProcessing\nStatus'    ]

        self.reprocess_datasets_table=QtGui.QTableWidget()
        self.reprocess_datasets_table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.reprocess_datasets_table.setSortingEnabled(True)
        self.reprocess_datasets_table.setColumnCount(len(self.reprocess_datasets_column_list))
        self.reprocess_datasets_table.setHorizontalHeaderLabels(self.reprocess_datasets_column_list)
        reprocess_vbox.addWidget(self.reprocess_datasets_table)

        # create context menu
        self.popMenu_for_reprocess_datasets_table = QtGui.QMenu()
        run_xia2_on_selected=QtGui.QAction("mark selected for reprocessing", self.window)
        run_xia2_on_selected.triggered.connect(self.select_sample_for_xia2)
        self.popMenu_for_reprocess_datasets_table.addAction(run_xia2_on_selected)
        self.reprocess_datasets_table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.reprocess_datasets_table.customContextMenuRequested.connect(self.on_context_menu_reprocess_data)

        frame=QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)
        hbox=QtGui.QHBoxLayout()

        frame_options=QtGui.QFrame()
        frame_options.setFrameShape(QtGui.QFrame.StyledPanel)
        hbox_options=QtGui.QHBoxLayout()
        label=QtGui.QLabel('Data processing protocol:')
        hbox_options.addWidget(label)
        self.xia2_3d_checkbox = QtGui.QCheckBox(' xia2 3d')
        hbox_options.addWidget(self.xia2_3d_checkbox)
        self.xia2_3dii_checkbox = QtGui.QCheckBox('xia2 3dii')
        hbox_options.addWidget(self.xia2_3dii_checkbox)
        self.xia2_dials_checkbox = QtGui.QCheckBox('Dials')
        hbox_options.addWidget(self.xia2_dials_checkbox)
        frame_options.setLayout(hbox_options)
        hbox.addWidget(frame_options)

        frame_sg=QtGui.QFrame()
        frame_sg.setFrameShape(QtGui.QFrame.StyledPanel)
        hbox_sg=QtGui.QHBoxLayout()
        label=QtGui.QLabel('Space Group:')
        hbox_sg.addWidget(label)
        self.reprocess_space_group_comboxbox=QtGui.QComboBox()
        self.reprocess_space_group_comboxbox.addItem('ignore')
        for sg in XChemMain.space_group_list():
            self.reprocess_space_group_comboxbox.addItem(sg)
        hbox_sg.addWidget(self.reprocess_space_group_comboxbox)
        frame_sg.setLayout(hbox_sg)
        hbox.addWidget(frame_sg)

        frame_ref=QtGui.QFrame()
        frame_ref.setFrameShape(QtGui.QFrame.StyledPanel)
        hbox_ref=QtGui.QHBoxLayout()
        label=QtGui.QLabel('Reference MTZ:')
        hbox_ref.addWidget(label)
        frame_ref_info=QtGui.QFrame()
        frame_ref_info.setFrameShape(QtGui.QFrame.StyledPanel)
        hbox_frame_ref_info=QtGui.QHBoxLayout()
        self.reprocess_reference_mtz_file_label=QtGui.QLabel(self.diffraction_data_reference_mtz)
        hbox_frame_ref_info.addWidget(self.reprocess_reference_mtz_file_label)
        frame_ref_info.setLayout(hbox_frame_ref_info)
        hbox_ref.addWidget(frame_ref_info)
        button=QtGui.QPushButton("Select")
        button.clicked.connect(self.select_reprocess_reference_mtz)
        hbox_ref.addWidget(button)
        frame_ref.setLayout(hbox_ref)
        hbox.addWidget(frame_ref)

        frame_isigma=QtGui.QFrame()
        frame_isigma.setFrameShape(QtGui.QFrame.StyledPanel)
        vbox_isigma=QtGui.QVBoxLayout()
        label=QtGui.QLabel('Resolution\nLimit:\nMn<I/sig(I)>')
        label.setAlignment(QtCore.Qt.AlignCenter)
        vbox_isigma.addWidget(label)
        self.reprocess_isigma_combobox=QtGui.QComboBox()
        misigma = ['default','4','3','2.5','2','1.5','1','0.5']
        for item in misigma:
            self.reprocess_isigma_combobox.addItem(item)
        self.reprocess_isigma_combobox.setCurrentIndex(0)
        self.reprocess_isigma_combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")
        vbox_isigma.addWidget(self.reprocess_isigma_combobox)
        frame_isigma.setLayout(vbox_isigma)
        hbox.addWidget(frame_isigma)

        frame_cc_half=QtGui.QFrame()
        frame_cc_half.setFrameShape(QtGui.QFrame.StyledPanel)
        vbox_cc_half=QtGui.QVBoxLayout()
        label=QtGui.QLabel('Resolution\nLimit:\nCC 1/2')
        label.setAlignment(QtCore.Qt.AlignCenter)
        vbox_cc_half.addWidget(label)
        self.reprocess_cc_half_combobox=QtGui.QComboBox()
        cc_half = ['default','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1']
        for item in cc_half:
            self.reprocess_cc_half_combobox.addItem(item)
        self.reprocess_cc_half_combobox.setCurrentIndex(0)
        self.reprocess_cc_half_combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")
        vbox_cc_half.addWidget(self.reprocess_cc_half_combobox)
        frame_cc_half.setLayout(vbox_cc_half)
        hbox.addWidget(frame_cc_half)

        hbox.addStretch(0)

        frame.setLayout(hbox)
        reprocess_vbox.addWidget(frame)

        self.dls_tab_dict['Reprocess'][1].addLayout(reprocess_vbox)

        #
        # @ MAP files Tab ####################################################################
        #

        initial_model_checkbutton_hbox=QtGui.QHBoxLayout()
        select_sample_for_dimple = QtGui.QCheckBox('(de-)select all samples for DIMPLE')
        select_sample_for_dimple.toggle()
        select_sample_for_dimple.setChecked(False)
        select_sample_for_dimple.stateChanged.connect(self.set_run_dimple_flag)
        initial_model_checkbutton_hbox.addWidget(select_sample_for_dimple)

        set_new_reference_button=QtGui.QPushButton("Set New Reference (if applicable)")
        set_new_reference_button.clicked.connect(self.set_new_reference_if_applicable)
        initial_model_checkbutton_hbox.addWidget(set_new_reference_button)

        refresh_reference_file_list_button=QtGui.QPushButton("Set New Reference (if applicable)")
        refresh_reference_file_list_button.clicked.connect(self.set_new_refresh_reference_file_list)
        initial_model_checkbutton_hbox.addWidget(refresh_reference_file_list_button)

        self.reference_file_list=self.get_reference_file_list(' ')
        self.reference_file_selection_combobox = QtGui.QComboBox()
        self.populate_reference_combobox(self.reference_file_selection_combobox)
        initial_model_checkbutton_hbox.addWidget(self.reference_file_selection_combobox)

        self.tab_dict[self.workflow_dict['Maps']][1].addLayout(initial_model_checkbutton_hbox)
        self.initial_model_vbox_for_table=QtGui.QVBoxLayout()
        self.inital_model_column_list=[     'Sample ID',
                                            'Select',
                                            'Compound ID',
                                            'Smiles',
                                            'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                            'Dimple\nRcryst',
                                            'Dimple\nRfree',
                                            'DataProcessing\nSpaceGroup',
                                            'Reference\nSpaceGroup',
                                            'Difference\nUC Volume (%)',
                                            'Reference File',
                                            'DataProcessing\nUnitCell',
                                            'Dimple\nStatus',
                                            'Compound\nStatus',
                                            'LastUpdated'                                ]

        self.initial_model_table=QtGui.QTableWidget()
        self.initial_model_table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.initial_model_table.setSortingEnabled(True)
        self.initial_model_table.setColumnCount(len(self.inital_model_column_list))
        self.initial_model_table.setHorizontalHeaderLabels(self.inital_model_column_list)
        self.initial_model_vbox_for_table.addWidget(self.initial_model_table)
        self.tab_dict[self.workflow_dict['Maps']][1].addLayout(self.initial_model_vbox_for_table)

        # create context menu
        self.popMenu_for_initial_model_table = QtGui.QMenu()
        run_dimple=QtGui.QAction("mark selected for dimple run", self.window)
        run_dimple.triggered.connect(self.select_sample_for_dimple)
        self.popMenu_for_initial_model_table.addAction(run_dimple)
        self.initial_model_table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.initial_model_table.customContextMenuRequested.connect(self.on_context_menu_initial_model)

        #
        # @ Refine Tab #######################################################################
        #

        self.summary_vbox_for_table=QtGui.QVBoxLayout()
        self.summary_column_name=[ 'Sample ID',
                                    'Compound ID',
                                    'Refinement\nSpace Group',
                                    'Refinement\nResolution',
                                    'Refinement\nRcryst',
                                    'Refinement\nRfree',
                                    'Refinement\nOutcome',
                                    'PanDDA site details',
                                    'Refinement\nStatus'    ]
        self.summary_table=QtGui.QTableWidget()
        self.summary_table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.summary_table.setSortingEnabled(True)
        self.summary_table.setColumnCount(len(self.summary_column_name))
        self.summary_table.setHorizontalHeaderLabels(self.summary_column_name)
        self.summary_vbox_for_table.addWidget(self.summary_table)
        self.tab_dict[self.workflow_dict['Refinement']][1].addLayout(self.summary_vbox_for_table)

        ######################################################################################

        #
        # @ DEPOSITION Tab ###################################################################
        #

        self.deposition_vbox=QtGui.QVBoxLayout()

        scroll = QtGui.QScrollArea()
        self.deposition_vbox.addWidget(scroll)
        scrollContent = QtGui.QWidget(scroll)

        scrollLayout = QtGui.QVBoxLayout(scrollContent)
        scrollContent.setLayout(scrollLayout)

        label_title=QtGui.QLabel('HTML export & ZENODO upload')
        label_title.setStyleSheet("font: 30pt Comic Sans MS")
        scrollLayout.addWidget(label_title)
        scrollLayout.addWidget(QtGui.QLabel(''))
        label_text=QtGui.QLabel(XChemToolTips.html_summary_introduction())
        label_text.setStyleSheet("font: 17pt Arial")
        scrollLayout.addWidget(label_text)
        scrollLayout.addWidget(QtGui.QLabel(''))
        image = QtGui.QLabel()
        pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'),'image','html_summary_page.png'))
        image.setPixmap(pixmap)
        scrollLayout.addWidget(image)
        scrollLayout.addWidget(QtGui.QLabel(''))
        label_heading=QtGui.QLabel('1. Specify HTML export directory in the settings tab')
        label_heading.setStyleSheet("font: bold 20pt Arial")
        scrollLayout.addWidget(label_heading)
        label_text=QtGui.QLabel(XChemToolTips.html_export_directory_background())
        label_text.setStyleSheet("font: 17pt Arial")
        scrollLayout.addWidget(label_text)
        if os.getcwd().startswith('/dls'):
            label_text=QtGui.QLabel('Note: default for labxchem project at DLS is <labxchem_directory>/processing/html.')
            label_text.setStyleSheet("font: 17pt Arial")
            scrollLayout.addWidget(label_text)
            label_text=QtGui.QLabel('In your case: '+self.html_export_directory)
            label_text.setStyleSheet("font: 17pt Arial")
            scrollLayout.addWidget(label_text)
        scrollLayout.addWidget(QtGui.QLabel('\n'))
        label_heading=QtGui.QLabel("2. Prepare files for HTML export")
        label_heading.setStyleSheet("font: bold 20pt Arial")
        scrollLayout.addWidget(label_heading)
        label_text=QtGui.QLabel(XChemToolTips.html_export_step())
        label_text.setStyleSheet("font: 17pt Arial")
        scrollLayout.addWidget(label_text)
        button=QtGui.QPushButton('Export to HTML')
        button.clicked.connect(self.export_to_html)
        button.setMaximumWidth(200)
        scrollLayout.addWidget(button)
        scrollLayout.addWidget(QtGui.QLabel('\n'))
        label_heading=QtGui.QLabel("3. Prepare ICB files")
        label_heading.setStyleSheet("font: bold 20pt Arial")
        scrollLayout.addWidget(label_heading)

        label_text=QtGui.QLabel(XChemToolTips.icb_file_background())
        label_text.setStyleSheet("font: 17pt Arial")
        scrollLayout.addWidget(label_text)


        label_text=QtGui.QLabel(XChemToolTips.prepare_ICB_files())
        label_text.setStyleSheet("font: 17pt Arial")
        scrollLayout.addWidget(label_text)
        button=QtGui.QPushButton('Open ICM-pro')
        button.clicked.connect(self.open_icm)
        button.setMaximumWidth(200)
        scrollLayout.addWidget(button)
        scrollLayout.addWidget(QtGui.QLabel(''))
        image=QtGui.QLabel()
        pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'),'image','drag_and_drop_icb_file.png'))
        image.setPixmap(pixmap)
        scrollLayout.addWidget(image)
        scrollLayout.addWidget(QtGui.QLabel(''))
        image=QtGui.QLabel()
        pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'),'image','run_icm_script.png'))
        image.setPixmap(pixmap)
        scrollLayout.addWidget(image)
        scrollLayout.addWidget(QtGui.QLabel('\n'))
        label_heading=QtGui.QLabel("4. Prepare files for ZENODO upload")
        label_heading.setStyleSheet("font: bold 20pt Arial")
        scrollLayout.addWidget(label_heading)
        label_text=QtGui.QLabel(XChemToolTips.zenodo_upload_start(self.html_export_directory))
        label_text.setStyleSheet("font: 17pt Arial")
        scrollLayout.addWidget(label_text)
        button=QtGui.QPushButton('prepare files')
        button.clicked.connect(self.prepare_files_for_zenodo_upload)
        button.setMaximumWidth(200)
        scrollLayout.addWidget(button)
        scrollLayout.addWidget(QtGui.QLabel('\n'))
        label_heading=QtGui.QLabel("5. ZENODO")
        label_heading.setStyleSheet("font: bold 20pt Arial")
        scrollLayout.addWidget(label_heading)
        label_text=QtGui.QLabel(XChemToolTips.zenodo_upload_part_one(self.html_export_directory))
        label_text.setStyleSheet("font: 17pt Arial")
        scrollLayout.addWidget(label_text)
        image=QtGui.QLabel()
        pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'),'image','new_zenodo_upload.png'))
        image.setPixmap(pixmap)
        scrollLayout.addWidget(image)
        scrollLayout.addWidget(QtGui.QLabel('\n'))
        label_heading=QtGui.QLabel("5. ZENODO upload ID")
        label_heading.setStyleSheet("font: bold 20pt Arial")
        scrollLayout.addWidget(label_heading)
        label_text=QtGui.QLabel(XChemToolTips.zenodo_upload_part_two())
        label_text.setStyleSheet("font: 17pt Arial")
        scrollLayout.addWidget(label_text)
        image=QtGui.QLabel()
        pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'),'image','zenodo_upload_id.png'))
        image.setPixmap(pixmap)
        scrollLayout.addWidget(image)
        label_text=QtGui.QLabel(XChemToolTips.zenodo_upload_part_three())
        label_text.setStyleSheet("font: 17pt Arial")
        scrollLayout.addWidget(label_text)
        hbox_zenodo_upload_id=QtGui.QHBoxLayout()
        hbox_zenodo_upload_id.addWidget(QtGui.QLabel('upload ID:'))
        self.zenodo_upload_id_entry = QtGui.QLineEdit()
        self.zenodo_upload_id_entry.setFixedWidth(200)
        hbox_zenodo_upload_id.addWidget(self.zenodo_upload_id_entry)
        hbox_zenodo_upload_id.addStretch(1)
        scrollLayout.addLayout(hbox_zenodo_upload_id)
        button=QtGui.QPushButton('update html files with upload ID')
        button.clicked.connect(self.update_html_for_zenodo_upload)
        button.setMaximumWidth(300)
        scrollLayout.addWidget(button)
        scrollLayout.addWidget(QtGui.QLabel('\n'))
        label_heading=QtGui.QLabel("6. ZENODO upload HTML files")
        label_heading.setStyleSheet("font: bold 20pt Arial")
        scrollLayout.addWidget(label_heading)
        label_text=QtGui.QLabel(XChemToolTips.zenodo_upload_part_four(self.html_export_directory))
        label_text.setStyleSheet("font: 17pt Arial")
        scrollLayout.addWidget(label_text)
        scrollLayout.addWidget(QtGui.QLabel('\n'))
        scrollLayout.addStretch(1)
        scroll.setWidget(scrollContent)
        self.tab_dict[self.workflow_dict['Deposition']][1].addLayout(self.deposition_vbox)

        ######################################################################################

        #
        # @ PANDDAs Tab ######################################################################
        #

        self.panddas_results_vbox=QtGui.QVBoxLayout()
        self.tab_dict[self.workflow_dict['PANDDAs']][1].addLayout(self.panddas_results_vbox)

        pandda_tab_widget = QtGui.QTabWidget()
        #pandda_tab_widget.setSizePolicy(size_policy)
        pandda_tab_list = [ 'pandda.analyse',
                            'Dataset Summary',
                            'Processing Output',
                            'pandda.inspect',
                            'Statistical Map Summaries']

        self.pandda_tab_dict={}
        for page in pandda_tab_list:
            tab=QtGui.QWidget()
            vbox=QtGui.QVBoxLayout(tab)
            pandda_tab_widget.addTab(tab,page)
            self.pandda_tab_dict[page]=[tab,vbox]

        self.pandda_analyse_hbox=QtGui.QHBoxLayout()
        self.pandda_tab_dict['pandda.analyse'][1].addLayout(self.pandda_analyse_hbox)
        self.pandda_map_layout = QtGui.QVBoxLayout()
        self.pandda_map_list = QtGui.QComboBox()
        self.pandda_maps_html = QtWebKit.QWebView()
        self.pandda_map_layout.addWidget(self.pandda_map_list)
        self.pandda_map_layout.addWidget(self.pandda_maps_html)

        self.pandda_tab_dict['Statistical Map Summaries'][1].addLayout(self.pandda_map_layout)
        self.pandda_maps_html.show()

        grid_pandda = QtGui.QGridLayout()
        grid_pandda.setColumnStretch(0,20)
        grid_pandda.setRowStretch(0,20)
        # left hand side: table with information about available datasets
        self.pandda_column_name = [ 'Sample ID',
                                    'Refinement\nSpace Group',
                                    'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                    'Dimple\nRcryst',
                                    'Dimple\nRfree',
                                    'Crystal Form\nName',
                                    #'PanDDA\nlaunched?',
                                    #'PanDDA\nhit?',
                                    #'PanDDA\nreject?'
                                    ]
                                    #'PanDDA\nStatus'    ]

        self.pandda_analyse_data_table=QtGui.QTableWidget()
        self.pandda_analyse_data_table.setSortingEnabled(True)
        self.pandda_analyse_data_table.resizeColumnsToContents()
        self.pandda_analyse_data_table.setColumnCount(len(self.pandda_column_name))
        self.pandda_analyse_data_table.setHorizontalHeaderLabels(self.pandda_column_name)

        frame_pandda=QtGui.QFrame()
        grid_pandda.addWidget(self.pandda_analyse_data_table,0,0)

        self.pandda_status = 'UNKNOWN'
        self.pandda_status_label = QtGui.QLabel()
        if os.path.exists(str(self.panddas_directory + '/pandda.done')):
            self.pandda_status = 'Finished!'
            self.pandda_status_label.setStyleSheet('color: green')
        if os.path.exists(str(self.panddas_directory + '/pandda.running')):
            self.pandda_status = 'Running...'
            self.pandda_status_label.setStyleSheet('color: orange')
        if os.path.exists(str(self.panddas_directory + '/pandda.errored')):
            self.pandda_status = 'Error encountered... please check the log files for pandda!'
            self.pandda_status_label.setStyleSheet('color: red')
        self.pandda_status_label.setText(str('STATUS: ' + self.pandda_status))
        self.pandda_status_label.setFont(QtGui.QFont("Arial",25, QtGui.QFont.Bold))
        grid_pandda.addWidget(self.pandda_status_label,3,0)

        # right hand side: input parameters for PANDDAs run
        frame_right=QtGui.QFrame()
        frame_right.setFrameShape(QtGui.QFrame.StyledPanel)

        self.pandda_analyse_input_params_vbox=QtGui.QVBoxLayout()

        pandda_input_dir_hbox=QtGui.QHBoxLayout()
        label=QtGui.QLabel('data directory')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        self.pandda_input_data_dir_entry = QtGui.QLineEdit()
        self.pandda_input_data_dir_entry.setText(os.path.join(self.initial_model_directory,'*'))
        self.pandda_input_data_dir_entry.setFixedWidth(300)
        pandda_input_dir_hbox.addWidget(self.pandda_input_data_dir_entry)
        self.select_pandda_input_dir_button=QtGui.QPushButton("Select Input Template")
        self.select_pandda_input_dir_button.setMaximumWidth(200)
        self.select_pandda_input_dir_button.clicked.connect(self.select_pandda_input_template)
        pandda_input_dir_hbox.addWidget(self.select_pandda_input_dir_button)
        self.pandda_analyse_input_params_vbox.addLayout(pandda_input_dir_hbox)

        pandda_pdb_style_hbox=QtGui.QHBoxLayout()
        label=QtGui.QLabel('pdb style')
        pandda_pdb_style_hbox.addWidget(label)
        self.pandda_pdb_style_entry=QtGui.QLineEdit()
        self.pandda_pdb_style_entry.setText('dimple.pdb')
        self.pandda_pdb_style_entry.setFixedWidth(200)
        pandda_pdb_style_hbox.addWidget(self.pandda_pdb_style_entry)
        self.pandda_analyse_input_params_vbox.addLayout(pandda_pdb_style_hbox)

        pandda_mtz_style_hbox=QtGui.QHBoxLayout()
        label=QtGui.QLabel('mtz style')
        pandda_mtz_style_hbox.addWidget(label)
        self.pandda_mtz_style_entry=QtGui.QLineEdit()
        self.pandda_mtz_style_entry.setText('dimple.mtz')
        self.pandda_mtz_style_entry.setFixedWidth(200)
        pandda_mtz_style_hbox.addWidget(self.pandda_mtz_style_entry)
        self.pandda_analyse_input_params_vbox.addLayout(pandda_mtz_style_hbox)

        pandda_output_dir_hbox=QtGui.QHBoxLayout()
        label=QtGui.QLabel('output directory')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        self.pandda_output_data_dir_entry = QtGui.QLineEdit()
        self.pandda_output_data_dir_entry.setText(self.panddas_directory)
        self.pandda_output_data_dir_entry.setFixedWidth(300)
        pandda_output_dir_hbox.addWidget(self.pandda_output_data_dir_entry)
        self.select_pandda_output_dir_button=QtGui.QPushButton("Select PANNDAs Directory")
        self.select_pandda_output_dir_button.setMaximumWidth(200)
        self.select_pandda_output_dir_button.clicked.connect(self.settings_button_clicked)
        pandda_output_dir_hbox.addWidget(self.select_pandda_output_dir_button)
        self.pandda_analyse_input_params_vbox.addLayout(pandda_output_dir_hbox)

        # qstat or local machine
        label=QtGui.QLabel('submit')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        self.pandda_submission_mode_selection_combobox = QtGui.QComboBox()
        if self.external_software['qsub']:
            self.pandda_submission_mode_selection_combobox.addItem('qsub')
        self.pandda_submission_mode_selection_combobox.addItem('local machine')
        self.pandda_submission_mode_selection_combobox.setMaximumWidth(200)
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_submission_mode_selection_combobox)

        label=QtGui.QLabel('number of processors')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        self.pandda_nproc=multiprocessing.cpu_count()-1
        self.pandda_nproc_entry = QtGui.QLineEdit()
        self.pandda_nproc_entry.setText(str(self.pandda_nproc).replace(' ',''))
        self.pandda_nproc_entry.setFixedWidth(200)
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_nproc_entry)

        label=QtGui.QLabel('order events by:')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        self.pandda_sort_event_combobox = QtGui.QComboBox()
        self.pandda_sort_event_combobox.addItem('cluster_size')
        self.pandda_sort_event_combobox.addItem('z_peak')
        self.pandda_sort_event_combobox.setMaximumWidth(200)
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_sort_event_combobox)

        # crystal form option
        label=QtGui.QLabel('Use space group of reference file as filter:')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        # reference file combobox, label with spg display
        hbox=QtGui.QHBoxLayout()
        self.reference_file_list=self.get_reference_file_list(' ')
        self.pandda_reference_file_selection_combobox = QtGui.QComboBox()
        self.populate_reference_combobox(self.pandda_reference_file_selection_combobox)
        self.pandda_reference_file_selection_combobox.activated[str].connect(self.change_pandda_spg_label)
        hbox.addWidget(self.pandda_reference_file_selection_combobox)
        self.pandda_reference_file_spg_label=QtGui.QLabel()
        hbox.addWidget(self.pandda_reference_file_spg_label)
        self.pandda_analyse_input_params_vbox.addLayout(hbox)

        label=QtGui.QLabel('\nExpert Parameters (only change if you know what you are doing!):')
        self.pandda_analyse_input_params_vbox.addWidget(label)

        self.wilson_checkbox = QtGui.QCheckBox('Wilson B-factor Scaling')
        self.wilson_checkbox.toggle()
        self.wilson_checkbox.setChecked(False)
        self.pandda_analyse_input_params_vbox.addWidget(self.wilson_checkbox)

        # minimum number of datasets
        label=QtGui.QLabel('min_build_datasets')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        self.pandda_min_build_dataset_entry = QtGui.QLineEdit()
        self.pandda_min_build_dataset_entry.setText('40')
        self.pandda_min_build_dataset_entry.setFixedWidth(200)
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_min_build_dataset_entry)

        label=QtGui.QLabel('max_new_datasets')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        self.pandda_max_new_datasets_entry = QtGui.QLineEdit()
        self.pandda_max_new_datasets_entry.setText('500')
        self.pandda_max_new_datasets_entry.setFixedWidth(200)
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_max_new_datasets_entry)

        label=QtGui.QLabel('grid_spacing (default=0.5)\nNote: higher values speed up calculations, but maps might be less pretty)')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        self.pandda_grid_spacing_entry = QtGui.QLineEdit()
        self.pandda_grid_spacing_entry.setText('0.5')
        self.pandda_grid_spacing_entry.setFixedWidth(200)
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_grid_spacing_entry)

        self.pandda_analyse_input_params_vbox.addStretch(1)

        frame_right.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Minimum)
        frame_right.setLayout(self.pandda_analyse_input_params_vbox)

        grid_pandda.addWidget(frame_right,0,1,5,5)
        frame_pandda.setLayout(grid_pandda)
        self.pandda_analyse_hbox.addWidget(frame_pandda)

        #######################################################
        # next three blocks display html documents created by pandda.analyse
        if os.path.exists(str(self.panddas_directory+'/interesting_datasets')):
            print('WARNING: USING RESULTS FROM OLD PANDDA ANALYSE! THIS IS NOT FULLY SUPPORTED IN XCE2')
            print('PLEASE CHANGE YOUR PANDDA DIRECTORY TO A NEW RUN, OR USE THE OLD VERSION OF XCE!')
            self.pandda_initial_html_file=str(self.panddas_directory+'/results_summareis/pandda_initial.html')
            self.pandda_analyse_html_file = str(self.panddas_directory + '/results_summaries/pandda_analyse.html')
        self.pandda_initial_html_file=str(self.panddas_directory+'/analyses/html_summaries/'+'pandda_initial.html')
        self.pandda_analyse_html_file=str(self.panddas_directory+'/analyses/html_summaries/'+'pandda_analyse.html')
        self.pandda_inspect_html_file=str(self.panddas_directory+'/analyses/html_summaries/'+'pandda_inspect.html')

        self.pandda_initial_html = QtWebKit.QWebView()
        self.pandda_tab_dict['Dataset Summary'][1].addWidget(self.pandda_initial_html)
        self.pandda_initial_html.load(QtCore.QUrl(self.pandda_initial_html_file))
        self.pandda_initial_html.show()

        self.pandda_analyse_html = QtWebKit.QWebView()
        self.pandda_tab_dict['Processing Output'][1].addWidget(self.pandda_analyse_html)
        self.pandda_analyse_html.load(QtCore.QUrl(self.pandda_analyse_html_file))
        self.pandda_analyse_html.show()

        self.pandda_inspect_html = QtWebKit.QWebView()
        self.pandda_tab_dict['pandda.inspect'][1].addWidget(self.pandda_inspect_html)
        self.pandda_analyse_html.load(QtCore.QUrl(self.pandda_inspect_html_file))
        self.pandda_analyse_html.show()

        self.panddas_results_vbox.addWidget(pandda_tab_widget)
        self.show_pandda_html_summary()

        ######################################################################################


        #
        # @ Settings Tab #####################################################################
        #

        ######################################################################################
        self.settings_container=QtGui.QWidget()
        self.buttons_etc = QtGui.QWidget()
        self.settings_vbox=QtGui.QVBoxLayout()

        self.scroll = QtGui.QScrollArea(self.settings_container)
        self.settings_vbox.addWidget(self.scroll)
        #scroll.setSizePolicy(size_policy)  #setWidgetResizable(True)
        scrollContent_settings = QtGui.QWidget(scroll)
        #scrollContent_settings.setSizePolicy(size_policy)

        scrollLayout_settings = QtGui.QVBoxLayout(scrollContent_settings)
        scrollContent_settings.setLayout(scrollLayout_settings)


        # Settings Tab
        self.data_collection_vbox_for_settings=QtGui.QVBoxLayout()

        self.buttons_etc.setLayout(self.data_collection_vbox_for_settings)
        self.scroll.setWidget(self.buttons_etc)

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nProject Directory: - REQUIRED -'))
        settings_hbox_initial_model_directory=QtGui.QHBoxLayout()
        self.initial_model_directory_label=QtGui.QLabel(self.initial_model_directory)
        settings_hbox_initial_model_directory.addWidget(self.initial_model_directory_label)
        settings_buttoon_initial_model_directory=QtGui.QPushButton('Select Project Directory')
        settings_buttoon_initial_model_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_initial_model_directory.addWidget(settings_buttoon_initial_model_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_initial_model_directory)


        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nReference Structure Directory: - OPTIONAL -'))
        settings_hbox_reference_directory=QtGui.QHBoxLayout()
        self.reference_directory_label=QtGui.QLabel(self.reference_directory)
        settings_hbox_reference_directory.addWidget(self.reference_directory_label)
        settings_buttoon_reference_directory=QtGui.QPushButton('Select Reference Structure Directory')
        settings_buttoon_reference_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_reference_directory.addWidget(settings_buttoon_reference_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_reference_directory)

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nData Source: - REQUIRED -'))
        settings_hbox_data_source_file=QtGui.QHBoxLayout()
        if self.data_source_file != '':
            self.data_source_file_label=QtGui.QLabel(os.path.join(self.database_directory,self.data_source_file))
        else:
            self.data_source_file_label=QtGui.QLabel('')
        settings_hbox_data_source_file.addWidget(self.data_source_file_label)
        settings_buttoon_data_source_file=QtGui.QPushButton('Select Data Source File')
        settings_buttoon_data_source_file.clicked.connect(self.settings_button_clicked)
        settings_hbox_data_source_file.addWidget(settings_buttoon_data_source_file)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_data_source_file)

        #################
        # Data Collection
        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nData Collection Directory: - OPTIONAL -'))

        settings_beamline_frame = QtGui.QFrame()
        settings_beamline_frame.setFrameShape(QtGui.QFrame.StyledPanel)
        settings_beamline_vbox = QtGui.QVBoxLayout()

        settings_hbox_beamline_directory=QtGui.QHBoxLayout()
        self.beamline_directory_label=QtGui.QLabel(self.beamline_directory)
        settings_hbox_beamline_directory.addWidget(self.beamline_directory_label)
        settings_buttoon_beamline_directory=QtGui.QPushButton('Select Data Collection Directory')
        settings_buttoon_beamline_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_beamline_directory.addWidget(settings_buttoon_beamline_directory)
        settings_beamline_vbox.addLayout(settings_hbox_beamline_directory)

        settings_hbox_data_collection_summary_file=QtGui.QHBoxLayout()
        self.data_collection_summary_file_label=QtGui.QLabel(self.data_collection_summary_file)
        settings_hbox_data_collection_summary_file.addWidget(self.data_collection_summary_file_label)
        settings_button_data_collection_summary_file=QtGui.QPushButton('Select Existing\nCollection Summary File')
        settings_button_data_collection_summary_file.clicked.connect(self.settings_button_clicked)
        settings_hbox_data_collection_summary_file.addWidget(settings_button_data_collection_summary_file)

        settings_button_new_data_collection_summary_file=QtGui.QPushButton('Assign New\nCollection Summary File')
        settings_button_new_data_collection_summary_file.clicked.connect(self.settings_button_clicked)
        settings_hbox_data_collection_summary_file.addWidget(settings_button_new_data_collection_summary_file)


        settings_beamline_vbox.addLayout(settings_hbox_data_collection_summary_file)

        settings_beamline_frame.setLayout(settings_beamline_vbox)
        self.data_collection_vbox_for_settings.addWidget(settings_beamline_frame)
        #################

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nCCP4_SCR Directory: - OPTIONAL -'))
        settings_hbox_ccp4_scratch_directory=QtGui.QHBoxLayout()
        self.ccp4_scratch_directory_label=QtGui.QLabel(self.ccp4_scratch_directory)
        settings_hbox_ccp4_scratch_directory.addWidget(self.ccp4_scratch_directory_label)
        settings_buttoon_ccp4_scratch_directory=QtGui.QPushButton('Select CCP4_SCR Directory')
        settings_buttoon_ccp4_scratch_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_ccp4_scratch_directory.addWidget(settings_buttoon_ccp4_scratch_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_ccp4_scratch_directory)

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nPANDDAs directory: - OPTIONAL -'))
        settings_hbox_panddas_directory=QtGui.QHBoxLayout()
        self.panddas_directory_label=QtGui.QLabel(self.panddas_directory)
        settings_hbox_panddas_directory.addWidget(self.panddas_directory_label)
        settings_button_panddas_directory=QtGui.QPushButton('Select PANNDAs Directory')
        settings_button_panddas_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_panddas_directory.addWidget(settings_button_panddas_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_panddas_directory)

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nHTML export directory: - OPTIONAL -'))
        settings_hbox_html_export_directory=QtGui.QHBoxLayout()
        self.html_export_directory_label=QtGui.QLabel(self.html_export_directory)
        settings_hbox_html_export_directory.addWidget(self.html_export_directory_label)
        settings_button_html_export_directory=QtGui.QPushButton('Select HTML Export Directory')
        settings_button_html_export_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_html_export_directory.addWidget(settings_button_html_export_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_html_export_directory)

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nGroup deposition directory: - OPTIONAL -'))
        settings_hbox_group_deposition_directory=QtGui.QHBoxLayout()
        self.group_deposition_directory_label=QtGui.QLabel(self.group_deposit_directory)
        settings_hbox_group_deposition_directory.addWidget(self.group_deposition_directory_label)
        settings_button_group_deposition_directory=QtGui.QPushButton('Select Group deposition Directory')
        settings_button_group_deposition_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_group_deposition_directory.addWidget(settings_button_group_deposition_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_group_deposition_directory)

        self.data_collection_vbox_for_settings.setSpacing(0)
        self.data_collection_vbox_for_settings.setContentsMargins(30,0,0,0)

	self.buttons_etc.resize(self.screen.width()-100, self.buttons_etc.sizeHint().height())
	self.tab_dict[self.workflow_dict['Settings']][1].addLayout(self.settings_vbox)

        ######################################################################################



        self.status_bar=QtGui.QStatusBar()
        self.progress_bar=QtGui.QProgressBar()
        self.progress_bar.setMaximum(100)
        self.status_bar.setMaximumWidth(self.screen.width())
	self.progress_bar.setMaximumWidth(self.screen.width())
	hbox_status=QtGui.QHBoxLayout()
        hbox_status.addWidget(self.status_bar)
        hbox_status.addWidget(self.progress_bar)

        vbox_main = QtGui.QVBoxLayout()
	menu_bar.setMaximumWidth(self.screen.width())
        vbox_main.addWidget(menu_bar)
        self.main_tab_widget.setMaximumSize(self.screen.width(),self.screen.height()-245)
	vbox_main.addWidget(self.main_tab_widget)

        hboxTaskFrames=QtGui.QHBoxLayout()
        update_from_datasource_button.setMaximumWidth((self.screen.width()-20)/5)
	frame_dataset_task.setMaximumWidth((self.screen.width()-20)/5)
	frame_map_cif_file_task.setMaximumWidth((self.screen.width()-20)/5)
	frame_panddas_file_task.setMaximumWidth((self.screen.width()-20)/5)
	frame_refine_file_task.setMaximumWidth((self.screen.width())-20/5)

	hboxTaskFrames.addWidget(update_from_datasource_button)
        hboxTaskFrames.addWidget(frame_dataset_task)
        hboxTaskFrames.addWidget(frame_map_cif_file_task)
        hboxTaskFrames.addWidget(frame_panddas_file_task)
        hboxTaskFrames.addWidget(frame_refine_file_task)
#        hboxTaskFrames.setSpacing(0)
#        hboxTaskFrames.setMargin(0)

        vbox_main.addLayout(hboxTaskFrames)

        vbox_main.addLayout(hbox_status)

        self.window.setLayout(vbox_main)

        self.status_bar.showMessage('Ready')
       # print self.window.minimumSize()
        self.window.show()

        if self.data_source_file != '':
            write_enabled=self.check_write_permissions_of_data_source()
            if not write_enabled:
                self.data_source_set=False

    def select_sample_for_dimple(self):
        indexes = self.initial_model_table.selectionModel().selectedRows()
        for index in sorted(indexes):
            xtal=str(self.initial_model_table.item(index.row(), 0).text())
            self.update_log.insert('{0!s} is marked for DIMPLE'.format(index.row()))
            self.initial_model_dimple_dict[xtal][0].setChecked(True)

    def select_sample_for_xia2(self):
        indexes = self.reprocess_datasets_table.selectionModel().selectedRows()
        for index in sorted(indexes):
            xtal=str(self.reprocess_datasets_table.item(index.row(), 1).text())
            print xtal,self.diffraction_data_table_dict[xtal][0]
            self.update_log.insert('{0!s} marked for reprocessing'.format(index.row()))
            self.diffraction_data_table_dict[xtal][0].setChecked(True)


    def update_summary_plot(self):
        if self.data_source_set:
            XChemPlots.summary_plot(os.path.join(self.database_directory,self.data_source_file),self.overview_axes).update_overview()
            self.overview_canvas.draw()


    def show_preferences(self):
        preferences = QtGui.QMessageBox()
        preferencesLayout = preferences.layout()

        vbox = QtGui.QVBoxLayout()
        settings_hbox_filename_root=QtGui.QHBoxLayout()
        filename_root_label=QtGui.QLabel('filename root:')
        settings_hbox_filename_root.addWidget(filename_root_label)
        filename_root_input = QtGui.QLineEdit()
        filename_root_input.setFixedWidth(400)
        filename_root_input.setText(str(self.filename_root))
        filename_root_input.textChanged[str].connect(self.change_filename_root)
        settings_hbox_filename_root.addWidget(filename_root_input)
        vbox.addLayout(settings_hbox_filename_root)

        settings_hbox_adjust_allowed_unit_cell_difference=QtGui.QHBoxLayout()
        adjust_allowed_unit_cell_difference_label=QtGui.QLabel('Max. Allowed Unit Cell Difference between Reference and Target (%):')
        settings_hbox_adjust_allowed_unit_cell_difference.addWidget(adjust_allowed_unit_cell_difference_label)
        adjust_allowed_unit_cell_difference = QtGui.QLineEdit()
        adjust_allowed_unit_cell_difference.setFixedWidth(200)
        adjust_allowed_unit_cell_difference.setText(str(self.allowed_unitcell_difference_percent))
        adjust_allowed_unit_cell_difference.textChanged[str].connect(self.change_allowed_unitcell_difference_percent)
        settings_hbox_adjust_allowed_unit_cell_difference.addWidget(adjust_allowed_unit_cell_difference)
        vbox.addLayout(settings_hbox_adjust_allowed_unit_cell_difference)

        settings_hbox_acceptable_low_resolution_limit=QtGui.QHBoxLayout()
        adjust_acceptable_low_resolution_limit_label=QtGui.QLabel('Acceptable low resolution limit for datasets (in Angstrom):')
        settings_hbox_acceptable_low_resolution_limit.addWidget(adjust_acceptable_low_resolution_limit_label)
        adjust_acceptable_low_resolution_limit = QtGui.QLineEdit()
        adjust_acceptable_low_resolution_limit.setFixedWidth(200)
        adjust_acceptable_low_resolution_limit.setText(str(self.acceptable_low_resolution_limit_for_data))
        adjust_acceptable_low_resolution_limit.textChanged[str].connect(self.change_acceptable_low_resolution_limit)
        settings_hbox_acceptable_low_resolution_limit.addWidget(adjust_acceptable_low_resolution_limit)
        vbox.addLayout(settings_hbox_acceptable_low_resolution_limit)

        vbox_data=QtGui.QVBoxLayout()
        vbox_data.addWidget(QtGui.QLabel('Select amount of processed data you wish to copy to initial_model directory:'))
        self.preferences_data_to_copy_combobox = QtGui.QComboBox()
        for item in self.preferences_data_to_copy:
            self.preferences_data_to_copy_combobox.addItem(item[0])
        self.preferences_data_to_copy_combobox.currentIndexChanged.connect(self.preferences_data_to_copy_combobox_changed)
        vbox_data.addWidget(self.preferences_data_to_copy_combobox)
        vbox.addLayout(vbox_data)

        vbox_select=QtGui.QVBoxLayout()
        vbox_select.addWidget(QtGui.QLabel('Dataset Selection Mechanism:'))
        self.preferences_selection_mechanism_combobox = QtGui.QComboBox()
        for item in self.preferences_selection_mechanism:
            self.preferences_selection_mechanism_combobox.addItem(item)
        self.preferences_selection_mechanism_combobox.currentIndexChanged.connect(self.preferences_selection_mechanism_combobox_changed)
        vbox_select.addWidget(self.preferences_selection_mechanism_combobox)
        vbox.addLayout(vbox_select)

        vbox_restraints=QtGui.QVBoxLayout()
        vbox_restraints.addWidget(QtGui.QLabel('Restraints generation program:'))
        self.preferences_restraints_generation_combobox = QtGui.QComboBox()
        program_list=[]
        if self.external_software['acedrg']:       program_list.append('acedrg')
        if self.external_software['phenix.elbow']: program_list.append('phenix.elbow')
        if self.external_software['grade']:        program_list.append('grade')
        for item in program_list:
            self.preferences_restraints_generation_combobox.addItem(item)
        self.preferences_restraints_generation_combobox.currentIndexChanged.connect(self.preferences_restraints_generation_combobox_changed)
        index = self.preferences_restraints_generation_combobox.findText(self.restraints_program, QtCore.Qt.MatchFixedString)
        self.preferences_restraints_generation_combobox.setCurrentIndex(index)
        vbox_restraints.addWidget(self.preferences_restraints_generation_combobox)
        vbox.addLayout(vbox_restraints)

        hbox=QtGui.QHBoxLayout()
        hbox.addWidget(QtGui.QLabel('XCE logfile:'))
        self.xce_logfile_label=QtGui.QLabel(self.xce_logfile)
        hbox.addWidget(self.xce_logfile_label)
        button=QtGui.QPushButton("Change")
        button.clicked.connect(self.set_xce_logfile)
        hbox.addWidget(button)
        vbox.addLayout(hbox)

        settings_hbox_max_queue_jobs=QtGui.QHBoxLayout()
        adjust_max_queue_jobs_label=QtGui.QLabel('Max. number of jobs running at once on DLS cluster:')
        settings_hbox_max_queue_jobs.addWidget(adjust_max_queue_jobs_label)
        adjust_max_queue_jobs = QtGui.QLineEdit()
        adjust_max_queue_jobs.setFixedWidth(200)
        adjust_max_queue_jobs.setText(str(self.max_queue_jobs))
        adjust_max_queue_jobs.textChanged[str].connect(self.change_max_queue_jobs)
        settings_hbox_max_queue_jobs.addWidget(adjust_max_queue_jobs)
        vbox.addLayout(settings_hbox_max_queue_jobs)

        settings_hbox_remote_qsub=QtGui.QHBoxLayout()
        remote_qsub_label=QtGui.QLabel('remote qsub:')
        settings_hbox_remote_qsub.addWidget(remote_qsub_label)
        self.remote_qsub_checkbox = QtGui.QCheckBox('use')
        self.remote_qsub_checkbox.toggled.connect(self.run_qsub_remotely)
        if self.using_remote_qsub_submission:
            self.remote_qsub_checkbox.setChecked(True)
        settings_hbox_remote_qsub.addWidget(self.remote_qsub_checkbox)
        self.remote_qsub_command = QtGui.QLineEdit()
        self.remote_qsub_command.setFixedWidth(550)
        self.remote_qsub_command.setText(self.remote_qsub_submission)
#        remote_qsub.textChanged[str].connect(self.change_max_queue_jobs)
        settings_hbox_remote_qsub.addWidget(self.remote_qsub_command)
        vbox.addLayout(settings_hbox_remote_qsub)


        preferencesLayout.addLayout(vbox,0,0)

        preferences.exec_();

    def run_qsub_remotely(self):
        self.remote_qsub_submission=str(self.remote_qsub_command.text())
        if self.remote_qsub_checkbox.isChecked():
            self.update_log.insert('submitting jobs to remote machine with: %s' %self.remote_qsub_submission)
            self.external_software['qsub_remote']=self.remote_qsub_submission
            self.using_remote_qsub_submission=True
            self.settings['remote_qsub']=self.remote_qsub_submission
        else:
            self.update_log.insert('switching off remote job submission')
            self.external_software['qsub_remote']=''
            self.settings['remote_qsub']=''
            self.using_remote_qsub_submission=False

    def enter_pdb_codes(self):
        pdbID_entry = QtGui.QMessageBox()
        pdbID_entryLayout = pdbID_entry.layout()

        vbox = QtGui.QVBoxLayout()

        frame=QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel('Text from PDB email'), 0,0)
        self.pdb_code_entry = QtGui.QTextEdit()
        self.pdb_code_entry.setText('')
        self.pdb_code_entry.setFixedWidth(500)
        grid.addWidget(self.pdb_code_entry, 1,0,20,1)

        frame.setLayout(grid)
        vbox.addWidget(frame)

        hbox=QtGui.QHBoxLayout()
        button=QtGui.QPushButton('Update Database')
        button.clicked.connect(self.update_database_with_pdb_codes)
        hbox.addWidget(button)

        vbox.addLayout(hbox)
        pdbID_entryLayout.addLayout(vbox,0,0)
        pdbID_entry.exec_();

    def export_to_html(self):
        self.update_log.insert('exporting contents of SQLite database into '+self.html_export_directory)
        os.system('ccp4-python '+os.getenv('XChemExplorer_DIR')+'/web/process_sqlite.py -t Summary -s '+os.path.join(self.database_directory,self.data_source_file)+' -d '+self.html_export_directory)
        XChemWeb.create_ICM_input_file(self.html_export_directory,os.path.join(self.database_directory,self.data_source_file))
        self.update_log.insert('open ICMpro:')
        self.update_log.insert('/dls/science/groups/i04-1/software/icm-3.8-5/icm64 -g')
        self.update_log.insert('open file browser and navigate to '+self.html_export_directory)
        self.update_log.insert('drag and drop dsEvent_sqlite.icm into the main window')
        self.update_log.insert('the script will appear in the Workspace Panel')
        self.update_log.insert('right click on the script and select RUN')
        self.update_log.insert('be patient, this may take a while, depending on the number of events')
        self.status_bar.showMessage('please check terminal window for further information')

    def open_icm(self):
        self.update_log.insert('starting ICM...')
        self.work_thread=XChemThread.start_ICM(self.html_export_directory)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def prepare_files_for_zenodo_upload(self):
        self.update_log.insert('preparing files for ZENODO upload...')
        os.system('ccp4-python '+os.getenv('XChemExplorer_DIR')+'/helpers/prepare_for_zenodo_upload.py '+self.html_export_directory)

    def update_html_for_zenodo_upload(self):
        try:
            uploadID=int(self.zenodo_upload_id_entry.text())
            self.update_log.insert('updating html files for ZENODO upload,...')
            self.update_log.insert('ZENODO upload = '+str(uploadID))
            os.system('ccp4-python '+os.getenv('XChemExplorer_DIR')+'/helpers/prepare_for_zenodo_upload.py {0!s} {1!s}'.format(self.html_export_directory, uploadID))
        except ValueError:
            self.update_log.insert('zenodo upload ID must be an integer!')

    def create_missing_apo_records_in_depositTable(self):
        self.db.create_missing_apo_records_for_all_structures_in_depositTable(self.initial_model_directory,self.xce_logfile)

    def update_file_information_of_apo_records(self):
        XChemDeposit.update_file_locations_of_apo_structuresin_DB(os.path.join(self.database_directory,self.data_source_file),self.initial_model_directory,self.xce_logfile)

    def prepare_models_for_deposition(self):

        for key in self.prepare_mmcif_files_dict:
            if self.sender() == self.prepare_mmcif_files_dict[key]:
                structureType=key

        overwrite_existing_mmcif=True
        self.work_thread=XChemDeposit.prepare_mmcif_files_for_deposition(   os.path.join(self.database_directory,self.data_source_file),
                                                                            self.xce_logfile,
                                                                            overwrite_existing_mmcif,
                                                                            self.initial_model_directory,
                                                                            structureType   )
        self.explorer_active=1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def prepare_for_group_deposition_upload(self):

        self.work_thread=XChemDeposit.prepare_for_group_deposition_upload(  os.path.join(self.database_directory,self.data_source_file),
                                                                            self.xce_logfile,
                                                                            self.group_deposit_directory   )
        self.explorer_active=1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def check_smiles_in_db_and_pdb(self):

        self.work_thread=XChemDeposit.compare_smiles_in_db_with_ligand_in_pdb(  self.initial_model_directory,
                                                                                os.path.join(self.database_directory,self.data_source_file),
                                                                                self.xce_logfile   )
        self.explorer_active=1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.connect(self.work_thread, QtCore.SIGNAL("show_error_dict"), self.show_error_dict)
        self.work_thread.start()

    def deposition_data(self):

        depositData = QtGui.QMessageBox()
        depositDataLayout = depositData.layout()


        vbox = QtGui.QVBoxLayout()

        deposit_tab_widget = QtGui.QTabWidget()
        deposit_tab_list = [ 'Contact',
                             'General',
                             'Authors',
                             'Citation',
                             'Molecule',
                             'Misc',
                             'Methods',
                             'Software' ]

        deposit_tab_dict={}
        for page in deposit_tab_list:
            tab=QtGui.QWidget()
            vb=QtGui.QVBoxLayout(tab)
            deposit_tab_widget.addTab(tab,page)
            deposit_tab_dict[page]=[tab,vb]


        #
        # PI & Scientist information
        #

        vb=QtGui.QVBoxLayout()
        hbox = QtGui.QHBoxLayout()

        frame=QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel('Principal Investigator'), 0,0)

        grid.addWidget(QtGui.QLabel('Salutation'), 1,0)
        self.contact_author_PI_salutation = QtGui.QLineEdit()
        self.contact_author_PI_salutation.setText('Dr.')
        self.contact_author_PI_salutation.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_salutation, 1,1)

        grid.addWidget(QtGui.QLabel('First name'), 2,0)
        self.contact_author_PI_first_name = QtGui.QLineEdit()
        self.contact_author_PI_first_name.setText('')
        self.contact_author_PI_first_name.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_first_name, 2,1)

        grid.addWidget(QtGui.QLabel('Last name'), 3,0)
        self.contact_author_PI_last_name = QtGui.QLineEdit()
        self.contact_author_PI_last_name.setText('')
        self.contact_author_PI_last_name.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_last_name, 3,1)

        grid.addWidget(QtGui.QLabel('Middle name'), 4,0)
        self.contact_author_PI_middle_name = QtGui.QLineEdit()
        self.contact_author_PI_middle_name.setText('')
        self.contact_author_PI_middle_name.setFixedWidth(200)
        self.contact_author_PI_middle_name.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.contact_author_PI_middle_name, 4,1)

        grid.addWidget(QtGui.QLabel('PI role'), 5,0)
        self.contact_author_PI_role = QtGui.QComboBox()
        PIroles = ['group leader','principal investigator','investigator']
        for item in PIroles: self.contact_author_PI_role.addItem(item)
        grid.addWidget(self.contact_author_PI_role, 5,1)

        grid.addWidget(QtGui.QLabel('Organization type'), 6,0)
        self.contact_author_PI_organization_type = QtGui.QComboBox()
        Organizations = ['academic','commercial','government']
        for item in Organizations: self.contact_author_PI_organization_type.addItem(item)
        grid.addWidget(self.contact_author_PI_organization_type, 6,1)

        grid.addWidget(QtGui.QLabel('Organization Name'), 7,0)
        self.contact_author_PI_organization_name = QtGui.QLineEdit()
        self.contact_author_PI_organization_name.setText('')
        self.contact_author_PI_organization_name.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_organization_name, 7,1)


        grid.addWidget(QtGui.QLabel('Email'), 8,0)
        self.contact_author_PI_email = QtGui.QLineEdit()
        self.contact_author_PI_email.setText('')
        self.contact_author_PI_email.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_email, 8,1)

        grid.addWidget(QtGui.QLabel('Street'), 9,0)
        self.contact_author_PI_address = QtGui.QLineEdit()
        self.contact_author_PI_address.setText('')
        self.contact_author_PI_address.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_address, 9,1)

        grid.addWidget(QtGui.QLabel('City'), 10,0)
        self.contact_author_PI_city = QtGui.QLineEdit()
        self.contact_author_PI_city.setText('')
        self.contact_author_PI_city.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_city, 10,1)

        grid.addWidget(QtGui.QLabel('State'), 11,0)
        self.contact_author_PI_State_or_Province = QtGui.QLineEdit()
        self.contact_author_PI_State_or_Province.setText('')
        self.contact_author_PI_State_or_Province.setFixedWidth(200)
        self.contact_author_PI_State_or_Province.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.contact_author_PI_State_or_Province, 11,1)

        grid.addWidget(QtGui.QLabel('ZIP code'), 12,0)
        self.contact_author_PI_Zip_Code = QtGui.QLineEdit()
        self.contact_author_PI_Zip_Code.setText('')
        self.contact_author_PI_Zip_Code.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_Zip_Code, 12,1)

        grid.addWidget(QtGui.QLabel('Country'), 13,0)
        self.contact_author_PI_Country = QtGui.QLineEdit()
        self.contact_author_PI_Country.setText('')
        self.contact_author_PI_Country.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_Country, 13,1)

        grid.addWidget(QtGui.QLabel('Phone'), 14,0)
        self.contact_author_PI_phone_number = QtGui.QLineEdit()
        self.contact_author_PI_phone_number.setText('')
        self.contact_author_PI_phone_number.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_phone_number, 14,1)

        frame.setLayout(grid)
        hbox.addWidget(frame)

        frame=QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)
        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel('Responsible Scientist'), 0,0)

        grid.addWidget(QtGui.QLabel('Salutation'), 1,0)
        self.contact_author_salutation = QtGui.QLineEdit()
        self.contact_author_salutation.setText('Dr.')
        self.contact_author_salutation.setFixedWidth(200)
        grid.addWidget(self.contact_author_salutation, 1,1)

        grid.addWidget(QtGui.QLabel('First name'), 2,0)
        self.contact_author_first_name = QtGui.QLineEdit()
        self.contact_author_first_name.setText('')
        self.contact_author_first_name.setFixedWidth(200)
        grid.addWidget(self.contact_author_first_name, 2,1)

        grid.addWidget(QtGui.QLabel('Last name'), 3,0)
        self.contact_author_last_name = QtGui.QLineEdit()
        self.contact_author_last_name.setText('')
        self.contact_author_last_name.setFixedWidth(200)
        grid.addWidget(self.contact_author_last_name, 3,1)

        grid.addWidget(QtGui.QLabel('Middle name'), 4,0)
        self.contact_author_middle_name = QtGui.QLineEdit()
        self.contact_author_middle_name.setText('')
        self.contact_author_middle_name.setFixedWidth(200)
        self.contact_author_middle_name.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.contact_author_middle_name, 4,1)

        grid.addWidget(QtGui.QLabel('Role'), 5,0)

        self.contact_author_role = QtGui.QComboBox()
        ScientistRoles = ['responsible scientist','investigator']
        for item in ScientistRoles: self.contact_author_role.addItem(item)
        grid.addWidget(self.contact_author_role, 5,1)

        grid.addWidget(QtGui.QLabel('Organization type'), 6,0)

        self.contact_author_organization_type = QtGui.QComboBox()
        for item in Organizations: self.contact_author_organization_type.addItem(item)
        grid.addWidget(self.contact_author_organization_type, 6,1)

        grid.addWidget(QtGui.QLabel('Organization Name'), 7,0)
        self.contact_author_organization_name = QtGui.QLineEdit()
        self.contact_author_organization_name.setText('')
        self.contact_author_organization_name.setFixedWidth(200)
        grid.addWidget(self.contact_author_organization_name, 7,1)

        grid.addWidget(QtGui.QLabel('Email'), 8,0)
        self.contact_author_email = QtGui.QLineEdit()
        self.contact_author_email.setText('')
        self.contact_author_email.setFixedWidth(200)
        grid.addWidget(self.contact_author_email, 8,1)

        grid.addWidget(QtGui.QLabel('Street'), 9,0)
        self.contact_author_address = QtGui.QLineEdit()
        self.contact_author_address.setText('')
        self.contact_author_address.setFixedWidth(200)
        grid.addWidget(self.contact_author_address, 9,1)

        grid.addWidget(QtGui.QLabel('City'), 10,0)
        self.contact_author_city = QtGui.QLineEdit()
        self.contact_author_city.setText('')
        self.contact_author_city.setFixedWidth(200)
        grid.addWidget(self.contact_author_city, 10,1)

        grid.addWidget(QtGui.QLabel('State'), 11,0)
        self.contact_author_State_or_Province = QtGui.QLineEdit()
        self.contact_author_State_or_Province.setText('')
        self.contact_author_State_or_Province.setFixedWidth(200)
        self.contact_author_State_or_Province.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.contact_author_State_or_Province, 11,1)

        grid.addWidget(QtGui.QLabel('ZIP code'), 12,0)
        self.contact_author_Zip_Code = QtGui.QLineEdit()
        self.contact_author_Zip_Code.setText('')
        self.contact_author_Zip_Code.setFixedWidth(200)
        grid.addWidget(self.contact_author_Zip_Code, 12,1)

        grid.addWidget(QtGui.QLabel('Country'), 13,0)
        self.contact_author_Country = QtGui.QLineEdit()
        self.contact_author_Country.setText('')
        self.contact_author_Country.setFixedWidth(200)
        grid.addWidget(self.contact_author_Country, 13,1)

        grid.addWidget(QtGui.QLabel('Phone'), 14,0)
        self.contact_author_phone_number = QtGui.QLineEdit()
        self.contact_author_phone_number.setText('')
        self.contact_author_phone_number.setFixedWidth(200)
        grid.addWidget(self.contact_author_phone_number, 14,1)

        frame.setLayout(grid)
        hbox.addWidget(frame)

        vb.addLayout(hbox)
        vb.addWidget(QtGui.QLabel(XChemToolTips.deposition_interface_note()))
        vb.addStretch(1)


        deposit_tab_dict['Contact'][1].addLayout(vb)

        #
        # Release status
        #

        vb=QtGui.QVBoxLayout()

        frame=QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel('Release status'), 0,0)

        grid.addWidget(QtGui.QLabel('Release Status for sequence'), 4,0)

        self.Release_status_for_sequence = QtGui.QComboBox()
        codeStatus = ['RELEASE NOW','HOLD FOR RELEASE']
        for item in codeStatus: self.Release_status_for_sequence.addItem(item)
        grid.addWidget(self.Release_status_for_sequence, 4,1)

        grid.addWidget(QtGui.QLabel('Release Status for coordinates/ SF'), 8,0)
        self.Release_status_for_coordinates = QtGui.QComboBox()
        coordStatus = ['RELEASE NOW','HOLD FOR PUBLICATION','HOLD FOR 4 WEEKS','HOLD FOR 6 MONTHS','HOLD FOR 1 YEAR']
        for item in coordStatus: self.Release_status_for_coordinates.addItem(item)
        grid.addWidget(self.Release_status_for_coordinates,                                        8,1)

        frame.setLayout(grid)
        vb.addWidget(frame)

        #
        # Release status
        #

        frame=QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel('Title & Details'), 0,0)
        note = ( 'Note: supported wildcards: $ProteinName,$CompoundName; e.g. "Crystal Structure of human JMJD2D in complex with N2317a"' )
        grid.addWidget(QtGui.QLabel(note), 1,0)

        grid.addWidget(QtGui.QLabel('Group deposition title'), 2,0)
        self.group_deposition_title = QtGui.QLineEdit()
        self.group_deposition_title.setText('PanDDA analysis group deposition')
        self.group_deposition_title.setFixedWidth(600)
#        self.group_deposition_title.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.group_deposition_title, 2,1)

        grid.addWidget(QtGui.QLabel('Description'), 3,0)
        self.group_description = QtGui.QLineEdit()
        self.group_description.setText('XDomainX of XOrganismX $ProteinName screened against the XXX Fragment Library by X-ray Crystallography at the XChem facility of Diamond Light Source beamline I04-1')
        self.group_description.setFixedWidth(600)
        grid.addWidget(self.group_description, 3,1)

        grid.addWidget(QtGui.QLabel('Structure Title (ligand bound)'), 4,0)
        self.structure_title = QtGui.QLineEdit()
        self.structure_title.setText('Crystal Structure of $ProteinName in complex with $CompoundName')
        self.structure_title.setFixedWidth(600)
        grid.addWidget(self.structure_title, 4,1)

        note = ( '\n\nApo Structure:\nonly use if you want to deposit PanDDA models!'        )
        grid.addWidget(QtGui.QLabel(note), 6,0)


        grid.addWidget(QtGui.QLabel('Structure Title (apo)'), 7,0)
        self.structure_title_apo = QtGui.QLineEdit()
        self.structure_title_apo.setText('Crystal Structure of $ProteinName after initial refinement with no ligand modelled (structure $n)')
        self.structure_title_apo.setFixedWidth(600)
        grid.addWidget(self.structure_title_apo, 7,1)


        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict['General'][1].addLayout(vb)

        #
        # Authors
        #

        vb=QtGui.QVBoxLayout()

        frame=QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel('Deposition authors (e.g. Surname, F.M.)'), 0,0)

        self.structure_author_name_List = []

        for column in range(0,2):
            for row in range(1,15):
                structure_author_name = QtGui.QLineEdit()
                structure_author_name.setText('')
                structure_author_name.setFixedWidth(300)
                grid.addWidget(structure_author_name, row,column)
                self.structure_author_name_List.append(structure_author_name)


        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict['Authors'][1].addLayout(vb)


        #
        # Primary citation
        #

        vb=QtGui.QVBoxLayout()

        frame=QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel('Primary Citation'), 0,0)

        grid.addWidget(QtGui.QLabel('ID'), 1,0)
        self.primary_citation_id = QtGui.QLineEdit()
        self.primary_citation_id.setText('primary')
        self.primary_citation_id.setFixedWidth(500)
        grid.addWidget(self.primary_citation_id, 1,1)

        grid.addWidget(QtGui.QLabel('Journal'), 2,0)
        self.primary_citation_journal_abbrev = QtGui.QLineEdit()
        self.primary_citation_journal_abbrev.setText('To be published')
        self.primary_citation_journal_abbrev.setFixedWidth(500)
        grid.addWidget(self.primary_citation_journal_abbrev, 2,1)

        grid.addWidget(QtGui.QLabel('Title'), 3,0)
        self.primary_citation_title = QtGui.QLineEdit()
        self.primary_citation_title.setText('')
        self.primary_citation_title.setFixedWidth(500)
        self.primary_citation_title.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.primary_citation_title, 3,1)

        grid.addWidget(QtGui.QLabel('Year'), 4,0)
        self.primary_citation_year = QtGui.QLineEdit()
        self.primary_citation_year.setText('')
        self.primary_citation_year.setFixedWidth(500)
        self.primary_citation_year.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.primary_citation_year, 4,1)

        grid.addWidget(QtGui.QLabel('Volume'), 5,0)
        self.primary_citation_journal_volume = QtGui.QLineEdit()
        self.primary_citation_journal_volume.setText('')
        self.primary_citation_journal_volume.setFixedWidth(500)
        self.primary_citation_journal_volume.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.primary_citation_journal_volume, 5,1)

        grid.addWidget(QtGui.QLabel('Page, first'), 6,0)
        self.primary_citation_page_first = QtGui.QLineEdit()
        self.primary_citation_page_first.setText('')
        self.primary_citation_page_first.setFixedWidth(500)
        self.primary_citation_page_first.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.primary_citation_page_first, 6,1)

        grid.addWidget(QtGui.QLabel('Page, last'), 7,0)
        self.primary_citation_page_last = QtGui.QLineEdit()
        self.primary_citation_page_last.setText('')
        self.primary_citation_page_last.setFixedWidth(500)
        self.primary_citation_page_last.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.primary_citation_page_last, 7,1)

        frame.setLayout(grid)
        vb.addWidget(frame)


        #
        # citation authors
        #

        frame=QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        set_primary_citation_authors = QtGui.QCheckBox('same as deposition authors')
        set_primary_citation_authors.toggle()
        set_primary_citation_authors.setChecked(False)
        set_primary_citation_authors.stateChanged.connect(self.set_primary_citation_as_structure_authors)
        grid.addWidget(set_primary_citation_authors, 0,0)

        self.primary_citation_author_name_List=[]

        for column in range(0,2):
            for row in range(1,15):
                primary_citation_author_name = QtGui.QLineEdit()
                primary_citation_author_name.setText('')
                primary_citation_author_name.setFixedWidth(300)
                grid.addWidget(primary_citation_author_name, row,column)
                self.primary_citation_author_name_List.append(primary_citation_author_name)

        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict['Citation'][1].addLayout(vb)

        #
        # Molecule Information
        #

        vb=QtGui.QVBoxLayout()

        frame=QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel('Entity 1'), 1,0)

        grid.addWidget(QtGui.QLabel('Molecule Name'), 2,0)
        self.molecule_name = QtGui.QLineEdit()
        self.molecule_name.setText('')
        self.molecule_name.setFixedWidth(300)
        self.molecule_name.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.molecule_name, 2,1)
        grid.addWidget(QtGui.QLabel('(e.g. RNA Hammerhead Ribozyme)'), 2,2)

        grid.addWidget(QtGui.QLabel('Fragment Name'), 3,0)
        self.fragment_name_one = QtGui.QLineEdit()
        self.fragment_name_one.setText('')
        self.fragment_name_one.setFixedWidth(300)
        self.fragment_name_one.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.fragment_name_one, 3,1)
        grid.addWidget(QtGui.QLabel('(e.g. ligand binding domain, hairpin)'), 3,2)

        grid.addWidget(QtGui.QLabel('Specific Mutation'), 4,0)
        self.fragment_name_one_specific_mutation = QtGui.QLineEdit()
        self.fragment_name_one_specific_mutation.setText('')
        self.fragment_name_one_specific_mutation.setFixedWidth(300)
        self.fragment_name_one_specific_mutation.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.fragment_name_one_specific_mutation, 4,1)
        grid.addWidget(QtGui.QLabel('(e.g. C280S)'), 4,2)

        grid.addWidget(QtGui.QLabel('Enzyme Comission Number'), 5,0)
        self.fragment_name_one_enzyme_comission_number = QtGui.QLineEdit()
        self.fragment_name_one_enzyme_comission_number.setText('')
        self.fragment_name_one_enzyme_comission_number.setFixedWidth(300)
        self.fragment_name_one_enzyme_comission_number.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.fragment_name_one_enzyme_comission_number, 5,1)
        grid.addWidget(QtGui.QLabel('(if known: e.g. 2.7.7.7)'), 5,2)

        grid.addWidget(QtGui.QLabel('Genetically Manipulated Source'), 6,0)

        grid.addWidget(QtGui.QLabel('Source organism scientific name'), 7,0)

        self.Source_organism_scientific_name = QtGui.QComboBox()
        taxonomy_dict=XChemMain.NCBI_taxonomy_ID()
        for item in taxonomy_dict:
            self.Source_organism_scientific_name.addItem(taxonomy_dict[item])
        grid.addWidget(self.Source_organism_scientific_name, 7,1)

        grid.addWidget(QtGui.QLabel('Source organism gene'), 8,0)
        self.Source_organism_gene = QtGui.QLineEdit()
        self.Source_organism_gene.setText('')
        self.Source_organism_gene.setFixedWidth(300)
        grid.addWidget(self.Source_organism_gene, 8,1)
        grid.addWidget(QtGui.QLabel('(e.g. RPOD, ALKA...)'), 8,2)

        grid.addWidget(QtGui.QLabel('Source organism strain'), 9,0)
        self.Source_organism_strain = QtGui.QLineEdit()
        self.Source_organism_strain.setText('')
        self.Source_organism_strain.setFixedWidth(300)
        self.Source_organism_strain.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.Source_organism_strain, 9,1)
        grid.addWidget(QtGui.QLabel('(e.g. BH10 ISOLATE, K-12...)'), 9,2)

        grid.addWidget(QtGui.QLabel('Expression system scientific name'), 10,0)

        self.Expression_system_scientific_name = QtGui.QComboBox()
        for item in taxonomy_dict:
            self.Expression_system_scientific_name.addItem(taxonomy_dict[item])
        grid.addWidget(self.Expression_system_scientific_name, 10,1)


        grid.addWidget(QtGui.QLabel('Expression system strain'), 11,0)
        self.Expression_system_strain = QtGui.QLineEdit()
        self.Expression_system_strain.setText('')
        self.Expression_system_strain.setFixedWidth(300)
        self.Expression_system_strain.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.Expression_system_strain, 11,1)
        grid.addWidget(QtGui.QLabel('(e.g. BL21(DE3))'), 11,2)

        grid.addWidget(QtGui.QLabel('Expression system vector type'), 12,0)
        self.Expression_system_vector_type = QtGui.QLineEdit()
        self.Expression_system_vector_type.setText('')
        self.Expression_system_vector_type.setFixedWidth(300)
        self.Expression_system_vector_type.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.Expression_system_vector_type, 12,1)
        grid.addWidget(QtGui.QLabel('(e.g. plasmid)'), 12,2)

        grid.addWidget(QtGui.QLabel('Expression_system_plasmid_name'), 13,0)
        self.Expression_system_plasmid_name = QtGui.QLineEdit()
        self.Expression_system_plasmid_name.setText('')
        self.Expression_system_plasmid_name.setFixedWidth(300)
        self.Expression_system_plasmid_name.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.Expression_system_plasmid_name, 13,1)
        grid.addWidget(QtGui.QLabel('(e.g. pET26)'), 13,2)

        grid.addWidget(QtGui.QLabel('Manipulated_source_details'), 14,0)
        self.Manipulated_source_details = QtGui.QLineEdit()
        self.Manipulated_source_details.setText('')
        self.Manipulated_source_details.setFixedWidth(300)
        self.Manipulated_source_details.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.Manipulated_source_details, 14,1)
        grid.addWidget(QtGui.QLabel('(any other relevant information)'), 14,2)

        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict['Molecule'][1].addLayout(vb)


        #
        # Misc
        #

        vb=QtGui.QVBoxLayout()

        frame=QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel('Keywords'), 1,0)
        self.structure_keywords = QtGui.QLineEdit()
        self.structure_keywords.setText('SGC - Diamond I04-1 fragment screening, PanDDA, XChemExplorer')
        self.structure_keywords.setFixedWidth(300)
        grid.addWidget(self.structure_keywords, 1,1)
        grid.addWidget(QtGui.QLabel('(e.g. beta barrel, protein-DNA complex)'), 1,2)

        grid.addWidget(QtGui.QLabel('Biological Assembly'), 2,0)
        self.biological_assembly_chain_number = QtGui.QLineEdit()
        self.biological_assembly_chain_number.setText('')
        self.biological_assembly_chain_number.setFixedWidth(300)
        grid.addWidget(self.biological_assembly_chain_number, 2,1)
        grid.addWidget(QtGui.QLabel('(e.g.  1 for monomer, 2 for dimer ..)'), 2,2)

        grid.addWidget(QtGui.QLabel('Sequence UNIPROT ID'), 3,0)
        self.molecule_one_letter_sequence_uniprot_id = QtGui.QLineEdit()
        self.molecule_one_letter_sequence_uniprot_id.setText('')
        self.molecule_one_letter_sequence_uniprot_id.setFixedWidth(300)
        grid.addWidget(self.molecule_one_letter_sequence_uniprot_id, 3,1)
        grid.addWidget(QtGui.QLabel('(e.g.  Q6B0I6)'), 3,2)

        grid.addWidget(QtGui.QLabel('Sequence'), 4,0)
        self.molecule_one_letter_sequence = QtGui.QTextEdit()
        self.molecule_one_letter_sequence.setText('')
        self.molecule_one_letter_sequence.setFixedWidth(300)
        grid.addWidget(self.molecule_one_letter_sequence, 4,1,7,2)

        grid.addWidget(QtGui.QLabel('Structural Genomic (optional)'), 8,0)

        grid.addWidget(QtGui.QLabel('Project Name'), 9,0)
        self.SG_project_name = QtGui.QLineEdit()
        self.SG_project_name.setText('')
        self.SG_project_name.setFixedWidth(300)
        grid.addWidget(self.SG_project_name, 9,1)
        grid.addWidget(QtGui.QLabel('(e.g. PSI, Protein Structure Initiative)'), 9,2)

        grid.addWidget(QtGui.QLabel('Full Name'), 10,0)
        self.full_name_of_SG_center = QtGui.QLineEdit()
        self.full_name_of_SG_center.setText('')
        self.full_name_of_SG_center.setFixedWidth(300)
        grid.addWidget(self.full_name_of_SG_center, 10,1)
        grid.addWidget(QtGui.QLabel('(e.g. Berkeley Structural Genomic Center)'), 10,2)


        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict['Misc'][1].addLayout(vb)


        #
        # Methods
        #

        vb=QtGui.QVBoxLayout()

        frame=QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel('Crystallization'), 1,0)

        grid.addWidget(QtGui.QLabel('Method'), 2,0)

        self.crystallization_method = QtGui.QComboBox()
        for item in XChemMain.crystal_growth_methods(): self.crystallization_method.addItem(item)
        grid.addWidget(self.crystallization_method, 2,1)

        grid.addWidget(QtGui.QLabel('pH'), 3,0)
        self.crystallization_pH = QtGui.QLineEdit()
        self.crystallization_pH.setText('')
        self.crystallization_pH.setFixedWidth(300)
        grid.addWidget(self.crystallization_pH, 3,1)
        grid.addWidget(QtGui.QLabel('(e.g. 7.5 ...)'), 3,2)

        grid.addWidget(QtGui.QLabel('Temperature'), 4,0)
        self.crystallization_temperature = QtGui.QLineEdit()
        self.crystallization_temperature.setText('')
        self.crystallization_temperature.setFixedWidth(300)
        grid.addWidget(self.crystallization_temperature, 4,1)
        grid.addWidget(QtGui.QLabel('(e.g. 298) (in Kelvin)'), 4,2)

        grid.addWidget(QtGui.QLabel('Condition'), 5,0)
        self.crystallization_details = QtGui.QLineEdit()
        self.crystallization_details.setText('')
        self.crystallization_details.setFixedWidth(300)
        grid.addWidget(self.crystallization_details, 5,1)
        grid.addWidget(QtGui.QLabel('(e.g. PEG 4000, NaCl etc.)'), 5,2)

        grid.addWidget(QtGui.QLabel('Diffraction Experiment'), 6,0)
        note = ( 'Note: this information will only be used if it is\n'
                 'not already available in the mainTable!\n'
                 'Ignore if data were collected at DLS' )
        grid.addWidget(QtGui.QLabel(note), 7,0)

        grid.addWidget(QtGui.QLabel('Source'), 8,0)

        self.radiation_source = QtGui.QComboBox()
        for item in XChemMain.radiationSource(): self.radiation_source.addItem(item)
        grid.addWidget(self.radiation_source, 8,1)

        grid.addWidget(QtGui.QLabel('Source Type'), 9,0)

        self.radiation_source_type = QtGui.QComboBox()
        for item in XChemMain.wwBeamlines(): self.radiation_source_type.addItem(item)
        grid.addWidget(self.radiation_source_type, 9,1)


        grid.addWidget(QtGui.QLabel('Wavelength'), 10,0)
        self.radiation_wavelengths = QtGui.QLineEdit()
        self.radiation_wavelengths.setText('')
        self.radiation_wavelengths.setFixedWidth(300)
        grid.addWidget(self.radiation_wavelengths, 10,1)
        grid.addWidget(QtGui.QLabel('(e.g. 1.502)'), 10,2)

        grid.addWidget(QtGui.QLabel('Detector'), 11,0)

        self.radiation_detector = QtGui.QComboBox()
        for item in XChemMain.detector(): self.radiation_detector.addItem(item)
        grid.addWidget(self.radiation_detector, 11,1)


        grid.addWidget(QtGui.QLabel('Detector Type'), 12,0)

        self.radiation_detector_type = QtGui.QComboBox()
        for item in XChemMain.detectorType(): self.radiation_detector_type.addItem(item)
        grid.addWidget(self.radiation_detector_type, 12,1)

        grid.addWidget(QtGui.QLabel('Date'), 13,0)
        self.data_collection_date = QtGui.QLineEdit()
        self.data_collection_date.setText('')
        self.data_collection_date.setFixedWidth(300)
        grid.addWidget(self.data_collection_date, 13,1)
        grid.addWidget(QtGui.QLabel('(e.g. 2004-01-07)'), 13,2)

        grid.addWidget(QtGui.QLabel('Temperature'), 14,0)
        self.data_collection_temperature = QtGui.QLineEdit()
        self.data_collection_temperature.setText('')
        self.data_collection_temperature.setFixedWidth(300)
        grid.addWidget(self.data_collection_temperature, 14,1)
        grid.addWidget(QtGui.QLabel('(e.g. 100) (in Kelvin)'), 14,2)

        grid.addWidget(QtGui.QLabel('Protocol'), 15,0)
        self.data_collection_protocol = QtGui.QLineEdit()
        self.data_collection_protocol.setText('SINGLE WAVELENGTH')
        self.data_collection_protocol.setFixedWidth(300)
        grid.addWidget(self.data_collection_protocol, 15,1)
        grid.addWidget(QtGui.QLabel('(e.g. SINGLE WAVELENGTH, MAD, ...)'), 15,2)

        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict['Methods'][1].addLayout(vb)


        #
        # Software
        #

        vb=QtGui.QVBoxLayout()

        frame=QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel('PDB starting model'), 1,0)
        self.pdbx_starting_model = QtGui.QLineEdit()
        self.pdbx_starting_model.setText('')
        self.pdbx_starting_model.setFixedWidth(300)
        grid.addWidget(self.pdbx_starting_model, 1,1)
        grid.addWidget(QtGui.QLabel('(e.g. 7.5 ...)'), 1,2)

        grid.addWidget(QtGui.QLabel('Data reduction'), 2,0)
        self.data_integration_software = QtGui.QComboBox()
        for item in XChemMain.data_integration_software(): self.data_integration_software.addItem(item)
        grid.addWidget(self.data_integration_software, 2,1)

        grid.addWidget(QtGui.QLabel('Phasing'), 3,0)
        self.phasing_software = QtGui.QComboBox()
        for item in XChemMain.phasing_software(): self.phasing_software.addItem(item)
        grid.addWidget(self.phasing_software, 3,1)

        frame.setLayout(grid)
        vb.addWidget(frame)
        vb.addStretch(1)

        deposit_tab_dict['Software'][1].addLayout(vb)

        vbox.addWidget(deposit_tab_widget)

        hbox=QtGui.QHBoxLayout()
        button=QtGui.QPushButton('Load\nFile')
        button.clicked.connect(self.load_deposit_config_file)
        hbox.addWidget(button)
        button=QtGui.QPushButton('Save\nFile')
        button.clicked.connect(self.save_deposit_config_file)
        hbox.addWidget(button)
        button=QtGui.QPushButton('Load from\nDatabase')
        button.clicked.connect(self.load_deposit_from_database)
        button.setEnabled(False)
        hbox.addWidget(button)
        button=QtGui.QPushButton('Save to\nDatabase')
        button.clicked.connect(self.save_deposit_to_database)
        hbox.addWidget(button)

        vbox.addLayout(hbox)
        depositDataLayout.addLayout(vbox,0,0)

        depositData.exec_()

    def save_deposit_config_file(self):
        self.update_deposit_dict()
        file_name = str(QtGui.QFileDialog.getSaveFileName(self.window,'Save file', self.current_directory))
        #make sure that the file always has .deposit extension
        if str(file_name).rfind('.') != -1:
            file_name=file_name[:file_name.rfind('.')]+'.deposit'
        else:
            file_name=file_name+'.deposit'
        pickle.dump(self.deposit_dict,open(file_name,'wb'))

    def update_database_with_pdb_codes(self):
        self.work_thread=XChemDeposit.import_PDB_IDs(   str(self.pdb_code_entry.toPlainText()),
                                                        os.path.join(self.database_directory,self.data_source_file),
                                                        self.xce_logfile   )
        self.explorer_active=1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()



    def load_deposit_config_file(self):
        file_name_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(self.window,'Open file', self.current_directory,'*.deposit')
        file_name=tuple(file_name_temp)[0]
        self.deposit_dict = pickle.load(open(file_name,"rb"))
        self.update_deposit_input()

    def load_deposit_from_database(self):
        print 'hallo'

    def save_deposit_to_database(self):
        self.update_deposit_dict()
        msgBox = QtGui.QMessageBox()
        msgBox.setText("*** WARNING ***\nAre you sure you want to update the database?\nThis will overwrite previous entries!")
        msgBox.addButton(QtGui.QPushButton('Yes'), QtGui.QMessageBox.YesRole)
        msgBox.addButton(QtGui.QPushButton('No'), QtGui.QMessageBox.RejectRole)
        reply = msgBox.exec_();
        if reply == 0:
            self.work_thread=XChemDeposit.update_depositTable(self.deposit_dict,
                                                              os.path.join(self.database_directory,self.data_source_file),
                                                              self.xce_logfile    )
            self.explorer_active=1
            self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.work_thread.start()

    def update_deposit_input(self):
        try:
            self.contact_author_PI_salutation.setText(self.deposit_dict['contact_author_PI_salutation'])
            self.contact_author_PI_first_name.setText(self.deposit_dict['contact_author_PI_first_name'])
            self.contact_author_PI_last_name.setText(self.deposit_dict['contact_author_PI_last_name'])
            self.contact_author_PI_middle_name.setText(self.deposit_dict['contact_author_PI_middle_name'])
            index = self.contact_author_PI_role.findText(self.deposit_dict['contact_author_PI_role'], QtCore.Qt.MatchFixedString)
            self.contact_author_PI_role.setCurrentIndex(index)
            index = self.contact_author_PI_organization_type.findText(self.deposit_dict['contact_author_PI_organization_type'], QtCore.Qt.MatchFixedString)
            self.contact_author_PI_organization_type.setCurrentIndex(index)
            self.contact_author_PI_organization_name.setText(self.deposit_dict['contact_author_PI_organization_name'])
            self.contact_author_PI_email.setText(self.deposit_dict['contact_author_PI_email'])
            self.contact_author_PI_address.setText(self.deposit_dict['contact_author_PI_address'])
            self.contact_author_PI_city.setText(self.deposit_dict['contact_author_PI_city'])
            self.contact_author_PI_State_or_Province.setText(self.deposit_dict['contact_author_PI_State_or_Province'])
            self.contact_author_PI_Zip_Code.setText(self.deposit_dict['contact_author_PI_Zip_Code'])
            self.contact_author_PI_Country.setText(self.deposit_dict['contact_author_PI_Country'])
            self.contact_author_PI_phone_number.setText(self.deposit_dict['contact_author_PI_phone_number'])

            self.contact_author_salutation.setText(self.deposit_dict['contact_author_salutation'])
            self.contact_author_first_name.setText(self.deposit_dict['contact_author_first_name'])
            self.contact_author_last_name.setText(self.deposit_dict['contact_author_last_name'])
            self.contact_author_middle_name.setText(self.deposit_dict['contact_author_middle_name'])
            index = self.contact_author_role.findText(self.deposit_dict['contact_author_role'], QtCore.Qt.MatchFixedString)
            self.contact_author_role.setCurrentIndex(index)
            index = self.contact_author_organization_type.findText(self.deposit_dict['contact_author_organization_type'], QtCore.Qt.MatchFixedString)
            self.contact_author_organization_type.setCurrentIndex(index)
            self.contact_author_organization_name.setText(self.deposit_dict['contact_author_organization_name'])
            self.contact_author_email.setText(self.deposit_dict['contact_author_email'])
            self.contact_author_address.setText(self.deposit_dict['contact_author_address'])
            self.contact_author_city.setText(self.deposit_dict['contact_author_city'])
            self.contact_author_State_or_Province.setText(self.deposit_dict['contact_author_State_or_Province'])
            self.contact_author_Zip_Code.setText(self.deposit_dict['contact_author_Zip_Code'])
            self.contact_author_Country.setText(self.deposit_dict['contact_author_Country'])
            self.contact_author_phone_number.setText(self.deposit_dict['contact_author_phone_number'])
            index = self.Release_status_for_coordinates.findText(self.deposit_dict['Release_status_for_coordinates'], QtCore.Qt.MatchFixedString)
            self.Release_status_for_coordinates.setCurrentIndex(index)
            index = self.Release_status_for_sequence.findText(self.deposit_dict['Release_status_for_sequence'], QtCore.Qt.MatchFixedString)
            self.Release_status_for_sequence.setCurrentIndex(index)

            self.group_deposition_title.setText(self.deposit_dict['group_deposition_title'])
            self.group_description.setText(self.deposit_dict['group_description'])

            self.structure_title.setText(self.deposit_dict['structure_title'])
            self.structure_title_apo.setText(self.deposit_dict['structure_title_apo'])

            for n,name in enumerate(self.deposit_dict['structure_author_name'].split(';')):
                self.structure_author_name_List[n].setText(name)

            self.primary_citation_id.setText(self.deposit_dict['primary_citation_id'])
            self.primary_citation_journal_abbrev.setText(self.deposit_dict['primary_citation_journal_abbrev'])
            self.primary_citation_title.setText(self.deposit_dict['primary_citation_title'])
            self.primary_citation_year.setText(self.deposit_dict['primary_citation_year'])
            self.primary_citation_journal_volume.setText(self.deposit_dict['primary_citation_journal_volume'])
            self.primary_citation_page_first.setText(self.deposit_dict['primary_citation_page_first'])
            self.primary_citation_page_last.setText(self.deposit_dict['primary_citation_page_last'])

            for n,name in enumerate(self.deposit_dict['primary_citation_author_name'].split(';')):
                self.primary_citation_author_name_List[n].setText(name)

            self.molecule_name.setText(self.deposit_dict['molecule_name'])
            self.fragment_name_one_specific_mutation.setText(self.deposit_dict['fragment_name_one_specific_mutation'])
            index = self.Source_organism_scientific_name.findText(self.deposit_dict['Source_organism_scientific_name'], QtCore.Qt.MatchFixedString)
            self.Source_organism_scientific_name.setCurrentIndex(index)

            self.Source_organism_gene.setText(self.deposit_dict['Source_organism_gene'])
            self.Source_organism_strain.setText(self.deposit_dict['Source_organism_strain'])
            index = self.Expression_system_scientific_name.findText(self.deposit_dict['Expression_system_scientific_name'], QtCore.Qt.MatchFixedString)
            self.Expression_system_scientific_name.setCurrentIndex(index)


            self.Expression_system_strain.setText(self.deposit_dict['Expression_system_strain'])
            self.Expression_system_vector_type.setText(self.deposit_dict['Expression_system_vector_type'])
            self.Expression_system_plasmid_name.setText(self.deposit_dict['Expression_system_plasmid_name'])
            self.Manipulated_source_details.setText(self.deposit_dict['Manipulated_source_details'])

            self.structure_keywords.setText(self.deposit_dict['structure_keywords'])
            self.biological_assembly_chain_number.setText(self.deposit_dict['biological_assembly_chain_number'])
            self.molecule_one_letter_sequence_uniprot_id.setText(self.deposit_dict['molecule_one_letter_sequence_uniprot_id'])
            self.molecule_one_letter_sequence.setText(self.deposit_dict['molecule_one_letter_sequence'])
            self.SG_project_name.setText(self.deposit_dict['SG_project_name'])
            self.full_name_of_SG_center.setText(self.deposit_dict['full_name_of_SG_center'])

            index = self.crystallization_method.findText(self.deposit_dict['crystallization_method'], QtCore.Qt.MatchFixedString)
            self.crystallization_method.setCurrentIndex(index)

            self.crystallization_pH.setText(self.deposit_dict['crystallization_pH'])
            self.crystallization_temperature.setText(self.deposit_dict['crystallization_temperature'])
            self.crystallization_details.setText(self.deposit_dict['crystallization_details'])
            index = self.radiation_source.findText(self.deposit_dict['radiation_source'], QtCore.Qt.MatchFixedString)
            self.radiation_source.setCurrentIndex(index)

            index = self.radiation_source_type.findText(self.deposit_dict['radiation_source_type'], QtCore.Qt.MatchFixedString)
            self.radiation_source_type.setCurrentIndex(index)

            self.radiation_wavelengths.setText(self.deposit_dict['radiation_wavelengths'])
            index = self.radiation_detector.findText(self.deposit_dict['radiation_detector'], QtCore.Qt.MatchFixedString)
            self.radiation_detector.setCurrentIndex(index)

            index = self.radiation_detector_type.findText(self.deposit_dict['radiation_detector_type'], QtCore.Qt.MatchFixedString)
            self.radiation_detector_type.setCurrentIndex(index)

            self.data_collection_date.setText(self.deposit_dict['data_collection_date'])
            self.data_collection_temperature.setText(self.deposit_dict['data_collection_temperature'])
            self.data_collection_protocol.setText(self.deposit_dict['data_collection_protocol'])

            self.pdbx_starting_model.setText(self.deposit_dict['pdbx_starting_model'])
            index = self.data_integration_software.findText(self.deposit_dict['data_integration_software'], QtCore.Qt.MatchFixedString)
            self.data_integration_software.setCurrentIndex(index)
            index = self.phasing_software.findText(self.deposit_dict['phasing_software'], QtCore.Qt.MatchFixedString)
            self.phasing_software.setCurrentIndex(index)

        except ValueError:
            self.update_status_bar('Sorry, this is not a XChemExplorer deposit file!')
            self.update_log.insert('Sorry, this is not a XChemExplorer deposit file!')



    def update_deposit_dict(self):
        self.deposit_dict = {
            'contact_author_PI_salutation':         str(self.contact_author_PI_salutation.text()),
            'contact_author_PI_first_name':         str(self.contact_author_PI_first_name.text()),
            'contact_author_PI_last_name':          str(self.contact_author_PI_last_name.text()),
            'contact_author_PI_middle_name':        str(self.contact_author_PI_middle_name.text()),
            'contact_author_PI_role':               str(self.contact_author_PI_role.currentText()),
            'contact_author_PI_organization_type':  str(self.contact_author_PI_organization_type.currentText()),
            'contact_author_PI_organization_name':  str(self.contact_author_PI_organization_name.text()),
            'contact_author_PI_email':              str(self.contact_author_PI_email.text()),
            'contact_author_PI_address':            str(self.contact_author_PI_address.text()),
            'contact_author_PI_city':               str(self.contact_author_PI_city.text()),
            'contact_author_PI_State_or_Province':  str(self.contact_author_PI_State_or_Province.text()),
            'contact_author_PI_Zip_Code':           str(self.contact_author_PI_Zip_Code.text()),
            'contact_author_PI_Country':            str(self.contact_author_PI_Country.text()),
            'contact_author_PI_phone_number':       str(self.contact_author_PI_phone_number.text()),

            'contact_author_salutation':            str(self.contact_author_salutation.text()),
            'contact_author_first_name':            str(self.contact_author_first_name.text()),
            'contact_author_last_name':             str(self.contact_author_last_name.text()),
            'contact_author_middle_name':           str(self.contact_author_middle_name.text()),
            'contact_author_role':                  str(self.contact_author_role.currentText()),
            'contact_author_organization_type':     str(self.contact_author_organization_type.currentText()),
            'contact_author_organization_name':     str(self.contact_author_organization_name.text()),
            'contact_author_email':                 str(self.contact_author_email.text()),
            'contact_author_address':               str(self.contact_author_address.text()),
            'contact_author_city':                  str(self.contact_author_city.text()),
            'contact_author_State_or_Province':     str(self.contact_author_State_or_Province.text()),
            'contact_author_Zip_Code':              str(self.contact_author_Zip_Code.text()),
            'contact_author_Country':               str(self.contact_author_Country.text()),
            'contact_author_phone_number':          str(self.contact_author_phone_number.text()),

            'Release_status_for_coordinates':       str(self.Release_status_for_coordinates.currentText()),
            'Release_status_for_sequence':          str(self.Release_status_for_sequence.currentText()),

            'group_deposition_title':               str(self.group_deposition_title.text()),
            'group_description':                    str(self.group_description.text()),

            'structure_title':                      str(self.structure_title.text()),
            'structure_title_apo':                  str(self.structure_title_apo.text()),

            'primary_citation_id':                  str(self.primary_citation_id.text()),
            'primary_citation_journal_abbrev':      str(self.primary_citation_journal_abbrev.text()),
            'primary_citation_title':               str(self.primary_citation_title.text()),
            'primary_citation_year':                str(self.primary_citation_year.text()),
            'primary_citation_journal_volume':      str(self.primary_citation_journal_volume.text()),
            'primary_citation_page_first':          str(self.primary_citation_page_first.text()),
            'primary_citation_page_last':           str(self.primary_citation_page_last.text()),

            'molecule_name':                                str(self.molecule_name.text()),
            'Source_organism_scientific_name':              str(self.Source_organism_scientific_name.currentText()),
            'Source_organism_gene':                         str(self.Source_organism_gene.text()),
            'Source_organism_strain':                       str(self.Source_organism_strain.text()),
            'Expression_system_scientific_name':            str(self.Expression_system_scientific_name.currentText()),
            'Expression_system_strain':                     str(self.Expression_system_strain.text()),
            'Expression_system_plasmid_name':               str(self.Expression_system_plasmid_name.text()),
            'Expression_system_vector_type':                str(self.Expression_system_vector_type.text()),
            'Manipulated_source_details':                   str(self.Manipulated_source_details.text()),
            'fragment_name_one_specific_mutation':          str(self.fragment_name_one_specific_mutation.text()),

            'structure_keywords':                           str(self.structure_keywords.text()),
            'biological_assembly_chain_number':             str(self.biological_assembly_chain_number.text()),
            'molecule_one_letter_sequence_uniprot_id':      str(self.molecule_one_letter_sequence_uniprot_id.text()),
            'SG_project_name':                              str(self.SG_project_name.text()),
            'full_name_of_SG_center':                       str(self.full_name_of_SG_center.text()),
            'molecule_one_letter_sequence':                 str(self.molecule_one_letter_sequence.toPlainText()).replace(' ','').replace('\n','').replace('\r',''),

            'crystallization_method':                       str(self.crystallization_method.currentText()),
            'crystallization_pH':                           str(self.crystallization_pH.text()),
            'crystallization_temperature':                  str(self.crystallization_temperature.text()),
            'crystallization_details':                      str(self.crystallization_details.text()),

            'radiation_source':                             str(self.radiation_source.currentText()),
            'radiation_source_type':                        str(self.radiation_source_type.currentText()),
            'radiation_wavelengths':                        str(self.radiation_wavelengths.text()),
            'radiation_detector':                           str(self.radiation_detector.currentText()),
            'radiation_detector_type':                      str(self.radiation_detector_type.currentText()),
            'data_collection_date':                         str(self.data_collection_date.text()),
            'data_collection_temperature':                  str(self.data_collection_temperature.text()),
            'data_collection_protocol':                     str(self.data_collection_protocol.text()),
            'pdbx_starting_model':                          str(self.pdbx_starting_model.text()),
            'data_integration_software':                    str(self.data_integration_software.currentText()),
            'phasing_software':                             str(self.phasing_software.currentText())
        }

        structure_author_name=''
        for widget in self.structure_author_name_List:
            structure_author_name+=str(widget.text())+';'
        self.deposit_dict['structure_author_name']=structure_author_name[:-1]

        primary_citation_author_name=''
        for widget in self.primary_citation_author_name_List:
            primary_citation_author_name+=str(widget.text())+';'
        self.deposit_dict['primary_citation_author_name']=primary_citation_author_name[:-1]

    def set_primary_citation_as_structure_authors(self,state):
        if state == QtCore.Qt.Checked:
            for n,entry in enumerate(self.structure_author_name_List):
                self.primary_citation_author_name_List[n].setText(str(entry.text()))
        else:
            for n,entry in enumerate(self.primary_citation_author_name_List):
                entry.setText('')

###################################################################################################

    def set_xce_logfile(self):
        file_name = str(QtGui.QFileDialog.getSaveFileName(self.window,'Save file', self.current_directory))
        self.xce_logfile=str(file_name)
        self.xce_logfile_label.setText(str(self.xce_logfile))
        if self.xce_logfile=='' or self.xce_logfile[self.xce_logfile.rfind('/')+1:]=='':
           print '==> XCE: invalid file format'
        else:
            XChemLog.startLog(self.xce_logfile).create_logfile(self.xce_version)
            self.update_log=XChemLog.updateLog(self.xce_logfile)


    def select_datasource_columns_to_display(self):
        columns_to_show = QtGui.QMessageBox()
        columns_to_showLayout = columns_to_show.layout()
        columns_in_data_source=self.db.return_column_list()
        try:
            columns_in_data_source=self.db.return_column_list()
        except AttributeError:
            print '==> XCE: please select a datasource file'
            self.status_bar.showMessage('please select a datasource file')
            return

        column_dict={}
        vbox = QtGui.QVBoxLayout()
        number_of_entries=len(columns_in_data_source)
        columns_shown_in_dialog_column=15
        grid = QtGui.QGridLayout()
        x=0
        y=0
        columns_to_ignore=self.db.columns_not_to_display()
        for entries_added in range(number_of_entries):
            if not columns_in_data_source[entries_added][1] in columns_to_ignore:
                data_source_column = QtGui.QCheckBox(columns_in_data_source[entries_added][1])
                column_dict[entries_added]=data_source_column
                if columns_in_data_source[entries_added][1] in self.data_source_columns_to_display:
                    data_source_column.setChecked(True)
                grid.addWidget(data_source_column, y,x)
                y+=1
            if y==columns_shown_in_dialog_column:
                y=0
                x+=1
        vbox.addLayout(grid)
        columns_to_showLayout.addLayout(vbox,0,0)

        columns_to_show.addButton(QtGui.QPushButton('OK'), QtGui.QMessageBox.YesRole)
        columns_to_show.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)
        reply=columns_to_show.exec_();
        if reply == 0:
            columns_to_show_list=['Sample ID']
            for key in column_dict:
                if column_dict[key].isChecked():
                    columns_to_show_list.append(columns_in_data_source[key][1])
            self.data_source_columns_to_display=columns_to_show_list
            self.populate_and_update_data_source_table()

    def update_header_and_data_from_datasource(self):
        self.update_log.insert('getting information for all samples from data source...')
        self.db=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
        self.update_log.insert('creating missing columns in data source')
        self.db.create_missing_columns()
        self.update_log.insert('load header and data from data source')
        self.header,self.data=self.db.load_samples_from_data_source()
        self.update_log.insert('get all samples in data source')
        all_samples_in_db=self.db.execute_statement("select CrystalName from mainTable where CrystalName is not '';")

        self.xtal_db_dict={}
        sampleID_column=0
        for n,entry in enumerate(self.header):
            if entry=='CrystalName':
                sampleID_column=n
                break
        for line in self.data:
            if str(line[sampleID_column]) != '':
                db_dict={}
                for n,entry in enumerate(line):
                    if n != sampleID_column:
                        db_dict[str(self.header[n])]=str(entry)
                self.xtal_db_dict[str(line[sampleID_column])]=db_dict

        print '==> XCE: found '+str(len(self.xtal_db_dict))+' samples'

    def datasource_menu_reload_samples(self):
        self.update_log.insert('reading samples from data source: '+os.path.join(self.database_directory,self.data_source_file))
        self.update_status_bar('reading samples from data source: '+os.path.join(self.database_directory,self.data_source_file))
        self.update_header_and_data_from_datasource()
        self.update_all_tables()

    def datasource_menu_save_samples(self):
        print 'hallo'

    def datasource_menu_export_csv_file(self):
        file_name = str(QtGui.QFileDialog.getSaveFileName(self.window,'Save file', self.database_directory))
        if file_name.rfind('.') != -1:
            file_name=file_name[:file_name.rfind('.')]+'.csv'
        else:
            file_name=file_name+'.csv'
        self.db.export_to_csv_file(file_name)

    def datasource_menu_import_csv_file(self):
        if self.data_source_set:
            file_name = QtGui.QFileDialog.getOpenFileName(self.window,'Open file', self.database_directory)
            self.db.import_csv_file(file_name)
        else:
            self.update_status_bar('Please load a data source file first')


    def datasource_menu_update_datasource(self):
        self.work_thread=XChemThread.synchronise_db_and_filesystem(self.initial_model_directory,os.path.join(self.database_directory,self.data_source_file),self.panddas_directory,self.xce_logfile,'project_directory')
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.connect(self.work_thread, QtCore.SIGNAL("datasource_menu_reload_samples"),self.datasource_menu_reload_samples)
        self.work_thread.start()

    def export_data_for_WONKA(self):
        self.update_log.insert('exporting CSV file for input into WONKA')
        self.db.export_csv_for_WONKA()


    def on_context_menu(self, point):
        # show context menu
        for key in self.dewar_configuration_dict:
            if self.dewar_configuration_dict[key]==self.sender():
                self.dewar_label_active=key
        self.popMenu.exec_(self.sender().mapToGlobal(point))

    def on_context_menu_initial_model(self, point):
        # show context menu
        self.popMenu_for_initial_model_table.exec_(self.sender().mapToGlobal(point))

    def on_context_menu_reprocess_data(self, point):
        # show context menu
        self.popMenu_for_reprocess_datasets_table.exec_(self.sender().mapToGlobal(point))

    def flag_sample_for_recollection(self):
        self.dewar_configuration_dict[self.dewar_label_active].setStyleSheet("background-color: yellow")

    def undo_flag_sample_for_recollection(self):
        self.dewar_configuration_dict[self.dewar_label_active].setStyleSheet("background-color: gray")

    def show_html_summary_in_firefox(self,xtal):
        html_summary=self.albula_button_dict[xtal][2]
        print 'html_summary',html_summary
        new=2
        webbrowser.open(html_summary,new=new)

    def update_pandda_crystal_from_combobox(self):
        self.pandda_analyse_crystal_from_selection_combobox.clear()
        self.pandda_analyse_crystal_from_selection_combobox.addItem('use all datasets')
        if os.path.isfile(os.path.join(self.database_directory,self.data_source_file)):
            self.load_crystal_form_from_datasource()
            if self.xtalform_dict != {}:
                print self.xtalform_dict
                for key in self.xtalform_dict:
                    self.pandda_analyse_crystal_from_selection_combobox.addItem(key)


    def populate_reference_combobox(self,combobox):
        combobox.clear()
        for reference_file in self.reference_file_list:
            combobox.addItem(reference_file[0])

    def set_new_refresh_reference_file_list(self):
        self.reference_file_list=self.get_reference_file_list(' ')
        self.populate_reference_combobox(self.reference_file_selection_combobox)


    def populate_refinement_outcome_combobox(self,combobox):
        combobox.clear()
        for stage in self.refinement_stage:
            combobox.addItem(stage)

    def change_pandda_spg_label(self):
        combo_text=str(self.pandda_reference_file_selection_combobox.currentText())
        for file in self.reference_file_list:
            if file[0] == combo_text:
                self.pandda_reference_file_spg_label.setText(file[1])
                break

    def populate_target_selection_combobox(self,combobox):
        combobox.clear()
        for target in self.target_list:
            combobox.addItem(target)

    def combo_selected(self, text):
        self.map_url = str(self.panddas_directory+'/analyses/html_summaries/pandda_map_' + text + '.html')
        self.pandda_maps_html.load(QtCore.QUrl(self.map_url))
        self.pandda_maps_html.show()

    def add_map_html(self):
        self.map_list = glob.glob(str(self.panddas_directory + '/analyses/html_summaries/pandda_map_*.html'))
        self.list_options = []
        for i in range(0, len(self.map_list)):
            string = self.map_list[i]
            string = string.replace('/analyses/html_summaries/pandda_map_', '')
            string = string.replace('.html', '')
            string = string.replace(self.panddas_directory, '')
            self.list_options.append(string)
        self.pandda_map_list.clear()
        for i in range(0, len(self.list_options)):
            self.pandda_map_list.addItem(self.list_options[i])
        self.connect(self.pandda_map_list, QtCore.SIGNAL('activated(QString)'), self.combo_selected)

    def open_config_file(self):
        file_name_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(self.window,'Open file', self.current_directory,'*.conf')
        file_name=tuple(file_name_temp)[0]
        try:
            pickled_settings = pickle.load(open(file_name,"rb"))
            if pickled_settings['beamline_directory'] != self.beamline_directory:
                self.beamline_directory=pickled_settings['beamline_directory']
                self.target_list,self.visit_list=XChemMain.get_target_and_visit_list(self.beamline_directory)
                self.settings['beamline_directory']=self.beamline_directory
                self.populate_target_selection_combobox(self.target_selection_combobox)

            self.initial_model_directory=pickled_settings['initial_model_directory']
            self.settings['initial_model_directory']=self.initial_model_directory

            self.panddas_directory=pickled_settings['panddas_directory']
            self.settings['panddas_directory']=self.panddas_directory
            if os.path.exists(str(self.panddas_directory + '/interesting_datasets')):
                print('WARNING: USING RESULTS FROM OLD PANDDA ANALYSE! THIS IS NOT FULLY SUPPORTED IN XCE2')
                print('PLEASE CHANGE YOUR PANDDA DIRECTORY TO A NEW RUN, OR USE THE OLD VERSION OF XCE!')
                self.pandda_initial_html_file = str(self.panddas_directory + '/results_summareis/pandda_initial.html')
                self.pandda_analyse_html_file = str(self.panddas_directory + '/results_summaries/pandda_analyse.html')
            self.pandda_initial_html_file=str(self.panddas_directory+'/analyses/html_summaries/'+'pandda_initial.html')
            self.pandda_analyse_html_file=str(self.panddas_directory+'/analyses/html_summaries/'+'pandda_analyse.html')

            self.pandda_inspect_html_file=str(self.panddas_directory+'/analyses/html_summaries/'+'pandda_inspect.html')
            self.show_pandda_html_summary()

            self.html_export_directory=pickled_settings['html_export_directory']
            self.html_export_directory_label.setText(self.html_export_directory)
            self.settings['html_export_directory']=self.html_export_directory
            self.group_deposit_directory=pickled_settings['group_deposit_directory']
            self.group_deposition_directory_label.setText(self.group_deposit_directory)
            self.settings['group_deposit_directory']=self.group_deposit_directory

            self.database_directory=pickled_settings['database_directory']
            self.settings['database_directory']=self.database_directory

            self.data_collection_summary_file=pickled_settings['data_collection_summary']
            self.data_collection_summary_file_label.setText(self.data_collection_summary_file)

            self.data_source_file=pickled_settings['data_source']
            if self.data_source_file != '':
                self.settings['data_source']=os.path.join(self.database_directory,self.data_source_file)
                # this is probably not necessary
                if os.path.isfile(self.settings['data_source']):
                    write_enabled=self.check_write_permissions_of_data_source()
                    if not write_enabled:
                        self.data_source_file_label.setText('')
                        self.data_source_set=False
                    else:
                        self.data_source_file_label.setText(os.path.join(self.database_directory,self.data_source_file))
                        self.data_source_set=True
                        self.db=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
                        self.datasource_menu_reload_samples()

            self.ccp4_scratch_directory=pickled_settings['ccp4_scratch']
            self.settings['ccp4_scratch']=self.ccp4_scratch_directory

            self.allowed_unitcell_difference_percent=pickled_settings['unitcell_difference']
            self.acceptable_low_resolution_limit_for_data=pickled_settings['too_low_resolution_data']

            reference_directory_temp=pickled_settings['reference_directory']
            if reference_directory_temp != self.reference_directory:
                self.reference_directory=reference_directory_temp
                self.settings['reference_directory']=self.reference_directory
                self.update_reference_files(' ')
                for xtal in self.initial_model_dimple_dict:
                    reference_file_selection_combobox=self.initial_model_dimple_dict[xtal][1]
                    self.populate_reference_combobox(reference_file_selection_combobox)

            self.initial_model_directory_label.setText(self.initial_model_directory)
            self.panddas_directory_label.setText(self.panddas_directory)
            self.pandda_output_data_dir_entry.setText(self.panddas_directory)
            self.reference_directory_label.setText(self.reference_directory)
            self.beamline_directory_label.setText(self.beamline_directory)
            self.ccp4_scratch_directory_label.setText(self.ccp4_scratch_directory)
            self.reference_file_list=self.get_reference_file_list(' ')


        except KeyError:
            self.update_status_bar('Sorry, this is not a XChemExplorer config file!')
            self.update_log.insert('Sorry, this is not a XChemExplorer config file!')

    def save_config_file(self):
        file_name = str(QtGui.QFileDialog.getSaveFileName(self.window,'Save file', self.current_directory))
        #make sure that the file always has .conf extension
        if str(file_name).rfind('.') != -1:
            file_name=file_name[:file_name.rfind('.')]+'.conf'
        else:
            file_name=file_name+'.conf'
        pickle.dump(self.settings,open(file_name,'wb'))


    def update_reference_files(self,reference_root):
        self.reference_file_list=self.get_reference_file_list(reference_root)
        self.populate_reference_combobox(self.reference_file_selection_combobox)
        self.populate_reference_combobox(self.pandda_reference_file_selection_combobox)

    def target_selection_combobox_activated(self,text):
        self.target=str(text)


    def check_status_rerun_dimple_on_all_autoprocessing_files(self):
        print 'hallo'


    def rerun_dimple_on_all_autoprocessing_files(self):
        job_list=[]
        self.update_log.insert('preparing to run DIMPLE on all autoprocessing files')
        for xtal in self.data_collection_dict:
            for entry in self.data_collection_dict[xtal]:
                if entry[0]=='logfile':
                    db_dict=entry[6]
                    try:
                        if os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'],db_dict['DataProcessingMTZfileName'])) or \
                            os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'])):
                            job_list=self.get_job_list_for_dimple_rerun(xtal,job_list,db_dict,entry)
                    except KeyError:
                        try:
                            if os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'])):
                                job_list=self.get_job_list_for_dimple_rerun(xtal,job_list,db_dict,entry)
                        except KeyError:
                            continue
        if job_list != []:
            self.update_log.insert('trying to run DIMPLE on ALL auto-processing files')
            self.check_before_running_dimple(job_list)

    def run_dimple_on_selected_autoprocessing_file(self):
        job_list=[]
        for xtal in sorted(self.initial_model_dimple_dict):
            if self.initial_model_dimple_dict[xtal][0].isChecked():
                db_dict=self.xtal_db_dict[xtal]

                # the if statement below is so convoluted, so that it is compatible with older data source files

                if os.path.isfile(os.path.join(db_dict['ProjectDirectory'],xtal,db_dict['DataProcessingPathToMTZfile'],db_dict['DataProcessingMTZfileName'])) or \
                   os.path.isfile(os.path.join(db_dict['ProjectDirectory'],xtal,db_dict['DataProcessingPathToMTZfile'])) or \
                   os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'],db_dict['DataProcessingMTZfileName'])) or \
                   os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'])):

                    if os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'],db_dict['DataProcessingMTZfileName'])):
                        mtzin=os.path.join(db_dict['DataProcessingPathToMTZfile'],db_dict['DataProcessingMTZfileName'])
                    elif os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'])):
                        mtzin=os.path.join(db_dict['DataProcessingPathToMTZfile'])
                    elif os.path.isfile(os.path.join(db_dict['ProjectDirectory'],xtal,db_dict['DataProcessingPathToMTZfile'],db_dict['DataProcessingMTZfileName'])):
                        mtzin=os.path.join(db_dict['ProjectDirectory'],xtal,db_dict['DataProcessingPathToMTZfile'],db_dict['DataProcessingMTZfileName'])
                    elif os.path.isfile(os.path.join(db_dict['ProjectDirectory'],xtal,db_dict['DataProcessingPathToMTZfile'])):
                        mtzin=os.path.join(db_dict['ProjectDirectory'],xtal,db_dict['DataProcessingPathToMTZfile'])


                    reference_file=str(self.initial_model_dimple_dict[xtal][1].currentText())

                    reference_file_pdb=os.path.join(self.reference_directory,reference_file+'.pdb')

                    if not os.path.isfile(reference_file_pdb):
                        continue

                    if os.path.isfile(os.path.join(self.reference_directory,reference_file+'.mtz')):
                        reference_file_mtz=' -R '+os.path.join(self.reference_directory,reference_file+'.mtz')
                    else:
                        reference_file_mtz=''

                    if os.path.isfile(os.path.join(self.reference_directory,reference_file+'.cif')):
                        reference_file_cif=' --libin '+os.path.join(self.reference_directory,reference_file+'.cif')
                    else:
                        reference_file_cif=''

                    job_list.append([   xtal,
                                        'dimple_rerun_on_selected_file',
                                        mtzin,
                                        reference_file_pdb,
                                        reference_file_mtz,
                                        reference_file_cif  ])

        if job_list != []:
            self.update_log.insert('trying to run DIMPLE on SELECTED auto-processing files')
            self.check_before_running_dimple(job_list)


    def remove_selected_dimple_files(self):
        job_list=[]
        for xtal in sorted(self.initial_model_dimple_dict):
            if self.initial_model_dimple_dict[xtal][0].isChecked():
                job_list.append(xtal)

        if job_list != []:
            msgBox = QtGui.QMessageBox()
            msgBox.setText("Do you really want to delete {0!s} Dimple files?".format(len(job_list)))
            msgBox.addButton(QtGui.QPushButton('Go'), QtGui.QMessageBox.YesRole)
            msgBox.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)
            reply = msgBox.exec_();

            if reply == 0:
                self.status_bar.showMessage('preparing to remove DIMPLE files')
                self.update_log.insert('preparing to remove DIMPLE files')
                self.work_thread=XChemThread.remove_selected_dimple_files(  job_list,
                                                                            self.initial_model_directory,
                                                                            self.xce_logfile,
                                                                            self.database_directory,
                                                                            self.data_source_file    )
                self.explorer_active=1
                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                self.connect(self.work_thread, QtCore.SIGNAL("datasource_menu_reload_samples"),self.datasource_menu_reload_samples)
                self.work_thread.start()



    def run_xia2_on_selected_datasets(self,overwrite):

        # check which programs should be run
        protocol=[]
        if self.xia2_3d_checkbox.isChecked():
            protocol.append('3d')
        if self.xia2_3dii_checkbox.isChecked():
            protocol.append('3dii')
        if self.xia2_dials_checkbox.isChecked():
            protocol.append('dials')

        # space group
        spg = []
        if str(self.reprocess_space_group_comboxbox.currentText()) != 'ignore':
            spg.append(str(self.reprocess_space_group_comboxbox.currentText()))

        # reference file
        ref = []
        if os.path.isfile(self.diffraction_data_reference_mtz):
            ref.append(self.diffraction_data_reference_mtz)

        # resolution limit
        reso_limit = []
        if str(self.reprocess_isigma_combobox.currentText()) != 'default':
            reso_limit.append(str(self.reprocess_isigma_combobox.currentText()))

        # cc 1/2
        cc_half = []
        if str(self.reprocess_cc_half_combobox.currentText()) != 'default':
            cc_half.append(str(self.reprocess_cc_half_combobox.currentText()))

        run_dict={}
        allRows = self.reprocess_datasets_table.rowCount()
        for row in xrange(0,allRows):
            dataset_id=str(self.reprocess_datasets_table.item(row,0).text())
            sample_id=str(self.reprocess_datasets_table.item(row,1).text())
            if self.diffraction_data_table_dict[dataset_id][0].isChecked():
                run_dict[sample_id]=self.diffraction_data_dict[dataset_id]

        if protocol != [] and run_dict !={}:
            self.work_thread=XChemProcess.run_xia2( self.initial_model_directory,
                                                    run_dict,
                                                    protocol,
                                                    spg,
                                                    ref,
                                                    reso_limit,
                                                    cc_half,
                                                    self.xce_logfile,
                                                    self.external_software,
                                                    self.ccp4_scratch_directory,
                                                    self.max_queue_jobs,
                                                    os.path.join(self.database_directory,self.data_source_file),
                                                    overwrite   )
            self.explorer_active=1
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.work_thread.start()
        else:
            self.update_log.insert('please select datasets and/ or data processing protocol')
            self.update_status_bar('please select datasets and/ or data processing protocol')

    def update_reprocessing_table(self):
        allRows = self.reprocess_datasets_table.rowCount()
        for row in xrange(0,allRows):
            sample_id=str(self.reprocess_datasets_table.item(row,1).text())
            if sample_id in self.xtal_db_dict:
                db_dict=self.xtal_db_dict[sample_id]
                cell_text=QtGui.QTableWidgetItem()
                cell_text.setText(db_dict['DataProcessingStatus'])
                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                if db_dict['DataProcessingStatus'] == 'running':
                    cell_text.setBackground(QtGui.QColor(100,230,150))
                elif db_dict['DataProcessingStatus'] == 'pending':
                    cell_text.setBackground(QtGui.QColor(20,100,230))
                elif db_dict['DataProcessingStatus'] == 'started':
                    cell_text.setBackground(QtGui.QColor(230,240,110))
                elif db_dict['DataProcessingStatus'] == 'finished':
                    cell_text.setBackground(QtGui.QColor(255,255,255))
                self.reprocess_datasets_table.setItem(row, 7, cell_text)

    def get_job_list_for_dimple_rerun(self,xtal,job_list,db_dict,entry):
        self.status_bar.showMessage('checking: '+str(os.path.join(db_dict['DataProcessingPathToMTZfile'],db_dict['DataProcessingMTZfileName'])))
        suitable_reference=[]
        for reference in self.reference_file_list:
            # first we need one in the same pointgroup
            if reference[5]==db_dict['DataProcessingPointGroup']:
                try:
                    difference=math.fabs(1-(float(db_dict['DataProcessingUnitCellVolume'])/float(reference[4])))
                    suitable_reference.append([reference[0],difference])
                except ValueError:
                    continue
        if suitable_reference != []:
            reference_file=min(suitable_reference,key=lambda x: x[1])[0]
            visit=entry[1]
            run=entry[2]
            autoproc=entry[4]

            reference_file_pdb=os.path.join(self.reference_directory,reference_file+'.pdb')

            if os.path.isfile(os.path.join(self.reference_directory,reference_file+'.mtz')):
                reference_file_mtz=' -R '+os.path.join(self.reference_directory,reference_file+'.mtz')
            else:
                reference_file_mtz=''

            if os.path.isfile(os.path.join(self.reference_directory,reference_file+'.cif')):
                reference_file_cif=' --libin '+os.path.join(self.reference_directory,reference_file+'.cif')
            else:
                reference_file_cif=''

            if os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'],db_dict['DataProcessingMTZfileName'])):
                mtzin=os.path.join(db_dict['DataProcessingPathToMTZfile'],db_dict['DataProcessingMTZfileName'])
            elif os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'])):
                mtzin=os.path.join(db_dict['DataProcessingPathToMTZfile'])

            self.update_log.insert('adding '+xtal+visit+'-'+run+autoproc+' to list')
            job_list.append([   xtal,
                                visit+'-'+run+autoproc,
                                mtzin,
                                reference_file_pdb,
                                reference_file_mtz,
                                reference_file_cif  ])
        self.status_bar.showMessage('idle')
        return job_list


    def check_before_running_dimple(self,job_list):

        msgBox = QtGui.QMessageBox()
        msgBox.setText("Do you really want to run {0!s} Dimple jobs?\nNote: we will not run more than 100 at once on the cluster!".format(len(job_list)))
        msgBox.addButton(QtGui.QPushButton('Go'), QtGui.QMessageBox.YesRole)
        msgBox.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)
        reply = msgBox.exec_();

        if reply == 0:
            self.status_bar.showMessage('preparing {0!s} DIMPLE jobs'.format(len(job_list)))
            self.update_log.insert('preparing to run {0!s} DIMPLE jobs'.format(len(job_list)))
            if self.external_software['qsub_array']:
                self.update_log.insert('we will be running an ARRAY job on the DLS computer cluster')
                self.update_log.insert('please note that the maximum number of jobs that will be running at once is {0!s}'.format(self.max_queue_jobs))
                self.update_log.insert('you can change this in the PREFERENCES menu, but be warned that to high a number might break the cluster!')
            self.update_log.insert('preparing input files for DIMPLE...')
            self.work_thread=XChemThread.run_dimple_on_all_autoprocessing_files(    job_list,
                                                                                    self.initial_model_directory,
                                                                                    self.external_software,
                                                                                    self.ccp4_scratch_directory,
                                                                                    self.database_directory,
                                                                                    self.data_source_file,
                                                                                    self.max_queue_jobs,
                                                                                    self.xce_logfile    )
            self.explorer_active=1
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.connect(self.work_thread, QtCore.SIGNAL("datasource_menu_reload_samples"),self.datasource_menu_reload_samples)
            self.work_thread.start()


    def select_pandda_input_template(self):
        filepath_temp=QtGui.QFileDialog.getOpenFileNameAndFilter(self.window,'Select Example PDB or MTZ File', self.initial_model_directory,'*.pdb;;*.mtz')
        filepath=str(tuple(filepath_temp)[0])
        pdbin=filepath.split('/')[-1]
        if filepath.endswith('.pdb'):
            pdbin=filepath.split('/')[-1]
            mtzin_temp=pdbin.replace('.pdb','.mtz')
            if os.path.isfile(filepath.replace(pdbin,mtzin_temp)):
                mtzin=mtzin_temp
            else:
                mtzin=''
        if filepath.endswith('.mtz'):
            mtzin=filepath.split('/')[-1]
            pdbin_temp=pdbin.replace('.mtz','.pdb')
            if os.path.isfile(filepath.replace(mtzin,pdbin_temp)):
                pdbin=pdbin_temp
            else:
                pdbin=''
        if len(filepath.split('/'))-len(self.initial_model_directory.split('/'))==2:
            self.pandda_input_data_dir_entry.setText(os.path.join(self.initial_model_directory,'*'))
        elif len(filepath.split('/'))-len(self.initial_model_directory.split('/')) > 2:
            subdir=os.path.join(*filepath.split('/')[len(self.initial_model_directory.split('/'))+1:len(filepath.split('/'))-1])
            self.pandda_input_data_dir_entry.setText(os.path.join(self.initial_model_directory,'*',subdir))
        else:
            pass
        self.pandda_pdb_style_entry.setText(pdbin)
        self.pandda_mtz_style_entry.setText(mtzin)

    def select_diffraction_data_directory(self):
        self.diffraction_data_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
        self.diffraction_data_dir_label.setText(self.diffraction_data_directory)
        self.settings['diffraction_data_directory']=self.diffraction_data_directory
        self.update_log.insert('setting diffraction data directory to '+self.diffraction_data_directory)

    def search_for_datasets(self):
        self.update_log.insert('search diffraction data directory for datasets...')
        self.work_thread=XChemMain.find_diffraction_image_directory_fast(self.diffraction_data_directory)
        self.explorer_active=1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.connect(self.work_thread, QtCore.SIGNAL("update_reprocess_datasets_table"),self.update_reprocess_datasets_table)
        self.work_thread.start()

    def translate_datasetID_to_sampleID(self):
        translate = QtGui.QMessageBox()
        translateLayout = translate.layout()
        self.translate_datasetID_to_sampleID_file='-'
        vbox = QtGui.QVBoxLayout()
        button=QtGui.QPushButton('Open CSV')
        button.clicked.connect(self.open_csv_file_translate_datasetID_to_sampleID)
        vbox.addWidget(button)
        self.translate_datasetID_to_sampleID_csv_label=QtGui.QLabel(self.translate_datasetID_to_sampleID_file)
        vbox.addWidget(self.translate_datasetID_to_sampleID_csv_label)
        translateLayout.addLayout(vbox,0,0)
        translate.addButton(QtGui.QPushButton('OK'), QtGui.QMessageBox.YesRole)
        translate.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)
        reply=translate.exec_();
        if reply == 0:
            if os.path.isfile(self.translate_datasetID_to_sampleID_file):
                trans_dict={}
                for line in open(self.translate_datasetID_to_sampleID_file):
                    if len(line.split(','))==2:
                        dataset=line.split(',')[0]
                        new_sample_id=line.split(',')[1]
                        trans_dict[dataset]=new_sample_id
                if len(trans_dict) >= 1:
                    allRows = self.reprocess_datasets_table.rowCount()
                    for row in xrange(0,allRows):
                        dataset_id=str(self.reprocess_datasets_table.item(row,0).text())
                        sample_id=str(self.reprocess_datasets_table.item(row,1).text())
                        if dataset_id in trans_dict:
                            cell_text=QtGui.QTableWidgetItem()
                            cell_text.setText(trans_dict[dataset_id])
                            cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                            self.reprocess_datasets_table.setItem(row, 1, cell_text)
                            self.update_log.insert('dataset: {0!s} -> changing sampleID to: {1!s}'.format(dataset_id, trans_dict[dataset_id]))


    def open_csv_file_translate_datasetID_to_sampleID(self):
        file_name_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(self.window,'Open file', self.current_directory,'*.csv')
        file_name=tuple(file_name_temp)[0]
        self.translate_datasetID_to_sampleID_csv_label.setText(file_name)
        self.translate_datasetID_to_sampleID_file=file_name

    def select_reprocess_reference_mtz(self):
        self.update_log.insert('trying to set new reference mtz file for reprocessing with xia2')
        file_name = str(QtGui.QFileDialog.getOpenFileName(self.window,'Select file', self.database_directory))
        if os.path.isfile(file_name):
            if file_name.endswith('.mtz'):
                self.diffraction_data_reference_mtz=file_name
                self.update_log.insert('new reference file for data processing with xia2: '+self.diffraction_data_reference_mtz)
                self.reprocess_reference_mtz_file_label.setText(self.diffraction_data_reference_mtz)
            else:
                self.update_log.insert('this does not seem to be a mtz file: '+file_name)

    def update_reprocess_datasets_table(self,data_dict):
        self.update_log.insert('updating reprocess datasets table')
        self.diffraction_data_table_dict={}
        self.diffraction_data_dict=data_dict

        self.diffraction_data_search_info='found '+str(len(self.diffraction_data_dict))+' datasets'
        self.diffraction_data_search_label.setText(self.diffraction_data_search_info)
        self.update_log.insert(self.diffraction_data_search_info)
        self.datasource_menu_reload_samples()
        # update table
        column_name=self.db.translate_xce_column_list_to_sqlite(self.reprocess_datasets_column_list)
        # set rows to 0
        self.reprocess_datasets_table.setRowCount(0)
        for entry in sorted(self.diffraction_data_dict):
            self.update_log.insert(str(self.diffraction_data_dict[entry]))
            if entry in self.xtal_db_dict:
                db_dict=self.xtal_db_dict[entry]
            else:
                db_dict={}
            row=self.reprocess_datasets_table.rowCount()
            self.reprocess_datasets_table.insertRow(row)
            for column,header in enumerate(column_name):
                if header[0]=='Dataset ID' or header[0]=='Sample ID':
                    cell_text=QtGui.QTableWidgetItem()
                    cell_text.setText(str(entry))
                    cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                    self.reprocess_datasets_table.setItem(row, column, cell_text)
                elif header[0]=='Run\nxia2':
                    run_xia2 = QtGui.QCheckBox()
                    run_xia2.toggle()
                    self.reprocess_datasets_table.setCellWidget(row, column, run_xia2)
                    run_xia2.setChecked(False)
                    self.diffraction_data_table_dict[entry]=[run_xia2]
                else:
                    cell_text=QtGui.QTableWidgetItem()
                    if db_dict != {}:
                        if header[0]=='DataProcessing\nStatus':
                            if str( db_dict[ header[1] ]  ) == 'running':
                                cell_text.setBackground(QtGui.QColor(100,230,150))
                            elif str( db_dict[ header[1] ]  ) == 'pending':
                                cell_text.setBackground(QtGui.QColor(20,100,230))
                            elif str( db_dict[ header[1] ]  ) == 'started':
                                cell_text.setBackground(QtGui.QColor(230,240,110))
                            elif str( db_dict[ header[1] ]  ) == 'finished':
                                cell_text.setBackground(QtGui.QColor(255,255,255))
                        cell_text.setText(str( db_dict[ header[1] ]  ))
                    else:
                        cell_text.setText('')
                    cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                    self.reprocess_datasets_table.setItem(row, column, cell_text)

    def update_all_tables(self):
        self.update_log.insert('checking for new reference files')
        self.update_status_bar('checking for new reference files')
        self.reference_file_list=self.get_reference_file_list(' ')
        self.update_log.insert('updating Overview table')
        self.update_status_bar('updating Overview table')
        self.populate_and_update_data_source_table()
        self.update_log.insert('updating Maps table')
        self.update_status_bar('updating Maps table')
        self.create_initial_model_table()
        self.update_log.insert('updating PANDDA table')
        self.update_status_bar('updating PANDDA table')
        self.populate_pandda_analyse_input_table()
        self.update_log.insert('updating REFINEMENT table')
        self.update_status_bar('updating REFINEMENT table')
        self.populate_and_update_refinement_table()
        self.update_log.insert('updating REPROCESSING table')
        self.update_status_bar('updating REPROCESSING table')
        self.update_reprocessing_table()
        self.update_status_bar('idle')
        self.update_summary_plot()

    def settings_button_clicked(self):
        if self.sender().text()=='Select Project Directory':
            self.initial_model_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.initial_model_directory_label.setText(self.initial_model_directory)
            self.settings['initial_model_directory']=self.initial_model_directory
        if self.sender().text()=='Select Reference Structure Directory':
            reference_directory_temp = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            if reference_directory_temp != self.reference_directory:
                self.reference_directory=reference_directory_temp
                self.update_reference_files(' ')
            self.reference_directory_label.setText(self.reference_directory)
            self.settings['reference_directory']=self.reference_directory
        if self.sender().text()=='Select Data Source File':
            filepath_temp=QtGui.QFileDialog.getOpenFileNameAndFilter(self.window,'Select File', self.database_directory,'*.sqlite')
            filepath=str(tuple(filepath_temp)[0])
            self.data_source_file =   filepath.split('/')[-1]
            self.database_directory = filepath[:filepath.rfind('/')]
            self.settings['database_directory']=self.database_directory
            self.settings['data_source']=os.path.join(self.database_directory,self.data_source_file)
            write_enabled=self.check_write_permissions_of_data_source()
            if not write_enabled:
                self.data_source_set=False
            else:
                self.data_source_set=True
                self.data_source_file_label.setText(os.path.join(self.database_directory,self.data_source_file))
                self.db=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
                self.db.create_missing_columns()
                self.datasource_menu_reload_samples()
        if self.sender().text()=='Select Data Collection Directory':
            dir_name = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            if dir_name != self.beamline_directory:
                self.beamline_directory=dir_name
                self.target_list,self.visit_list=XChemMain.get_target_and_visit_list(self.beamline_directory)
                self.populate_target_selection_combobox(self.target_selection_combobox)
            self.beamline_directory_label.setText(self.beamline_directory)
            self.settings['beamline_directory']=self.beamline_directory

        if self.sender().text()=='Select Existing\nCollection Summary File':
            if self.data_collection_summary_file != '':
                filepath_temp=QtGui.QFileDialog.getOpenFileNameAndFilter(self.window,'Select File', self.data_collection_summary_file[:self.data_collection_summary_file.rfind('/')],'*.pkl')
            else:
                filepath_temp=QtGui.QFileDialog.getOpenFileNameAndFilter(self.window,'Select File', os.getcwd(),'*.pkl')
            filepath=str(tuple(filepath_temp)[0])
            self.data_collection_summary_file=filepath
            self.data_collection_summary_file_label.setText(self.data_collection_summary_file)
            self.settings['data_collection_summary']=self.data_collection_summary_file

        if self.sender().text()=='Assign New\nCollection Summary File':
            if self.data_collection_summary_file != '':
                file_name = str(QtGui.QFileDialog.getSaveFileName(self.window,'New file', self.data_collection_summary_file[:self.data_collection_summary_file.rfind('/')]))
            else:
                file_name = str(QtGui.QFileDialog.getSaveFileName(self.window,'New file', self.current_directory))
            #make sure that the file always has .pkl extension
            if str(file_name).rfind('.') != -1:
                file_name=file_name[:file_name.rfind('.')]+'.pkl'
            else:
                file_name=file_name+'.pkl'
            self.data_collection_summary_file=file_name
            self.data_collection_summary_file_label.setText(self.data_collection_summary_file)
            self.settings['data_collection_summary']=self.data_collection_summary_file


        if self.sender().text()=='Select CCP4_SCR Directory':
            self.ccp4_scratch_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.ccp4_scratch_directory_label.setText(self.ccp4_scratch_directory)
            self.settings['ccp4_scratch']=self.ccp4_scratch_directory
        if self.sender().text()=='Select PANNDAs Directory':
            self.panddas_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.panddas_directory_label.setText(self.panddas_directory)
            self.pandda_output_data_dir_entry.setText(self.panddas_directory)
            print 'PANDDA',self.panddas_directory
            self.settings['panddas_directory']=self.panddas_directory
            if os.path.exists(str(self.panddas_directory + '/interesting_datasets')):
                print('WARNING: USING RESULTS FROM OLD PANDDA ANALYSE! THIS IS NOT FULLY SUPPORTED IN XCE2')
                print('PLEASE CHANGE YOUR PANDDA DIRECTORY TO A NEW RUN, OR USE THE OLD VERSION OF XCE!')
                self.pandda_initial_html_file = str(self.panddas_directory + '/results_summareis/pandda_initial.html')
                self.pandda_analyse_html_file = str(self.panddas_directory + '/results_summaries/pandda_analyse.html')
            self.pandda_initial_html_file=str(self.panddas_directory+'/analyses/html_summaries/'+'pandda_initial.html')
            self.pandda_analyse_html_file=str(self.panddas_directory+'/analyses/html_summaries/'+'pandda_analyse.html')
            self.pandda_inspect_html_file=str(self.panddas_directory+'/analyses/html_summaries/'+'pandda_inspect.html')

            # update add lead option for proasis if pandda directory is changed
            if os.path.isfile(os.path.join(self.panddas_directory, 'analyses/pandda_analyse_sites.csv')):
                # hide old menu info
                self.proasis_lead.setVisible(False)
                # enable lead adding if pandda_analyse_sites.csv now exists
                self.proasis_lead = QtGui.QAction(str('Create lead from pandda sites...'), self.window)
                self.proasis_lead.triggered.connect(lambda:self.add_lead())
                self.proasis_menu.addAction(self.proasis_lead)
            else:
                # otherwise, keep same as old menu
                self.proasis_lead = QtGui.QAction(str('Site info not found... please run pandda analyse before adding lead'),
                                         self.window)
                self.proasis_menu.addAction(self.proasis_lead)

        if self.sender().text()=='Select HTML Export Directory':
            self.html_export_directory=str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.html_export_directory_label.setText(self.html_export_directory)
            self.settings['html_export_directory']=self.html_export_directory

        if self.sender().text()=='Select Group deposition Directory':
            self.group_deposit_directory=str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.group_deposition_directory_label.setText(self.group_deposit_directory)
            self.settings['group_deposit_directory']=self.group_deposit_directory


    def change_allowed_unitcell_difference_percent(self,text):
        try:
            self.allowed_unitcell_difference_percent=int(text)
            self.settings['unitcell_difference']=self.allowed_unitcell_difference_percent
            self.update_log.insert('changing max allowed unit cell difference between reference and xtal to {0!s} percent'.format(self.allowed_unitcell_difference_percent))
        except ValueError:
            if str(text).find('.') != -1:
                self.allowed_unitcell_difference_percent=int(str(text)[:str(text).find('.')])
                self.settings['unitcell_difference']=self.allowed_unitcell_difference_percent
                self.update_log.insert('changing max allowed unit cell difference between reference and xtal to {0!s} percent'.format(self.allowed_unitcell_difference_percent))
            else:
                pass

    def change_max_queue_jobs(self,text):
        try:
            self.max_queue_jobs=int(text)
            self.settings['max_queue_jobs']=self.max_queue_jobs
            self.update_log.insert('changing max number of jobs running simultaneously on DLS cluster to {0!s}'.format(self.max_queue_jobs))
        except ValueError:
            if str(text).find('.') != -1:
                self.max_queue_jobs=int(str(text)[:str(text).find('.')])
                self.settings['max_queue_jobs']=self.max_queue_jobs
                self.update_log.insert('changing max number of jobs running simultaneously on DLS cluster to {0!s}'.format(self.max_queue_jobs))
            else:
                pass

    def change_acceptable_low_resolution_limit(self,text):
        try:
            self.acceptable_low_resolution_limit_for_data=float(text)
            self.settings['too_low_resolution_data']=self.acceptable_low_resolution_limit_for_data
        except ValueError:
            pass

    def change_filename_root(self,text):
        self.filename_root=str(text)
        self.settings['filename_root']=self.filename_root

    def create_new_data_source(self):
        file_name = str(QtGui.QFileDialog.getSaveFileName(self.window,'Save file', self.database_directory))
        #make sure that the file always has .sqlite extension
        if file_name.rfind('.') != -1:
            file_name=file_name[:file_name.rfind('.')]+'.sqlite'
        else:
            file_name=file_name+'.sqlite'
        self.db=XChemDB.data_source(file_name)
        print '==> XCE: creating new data source'
        self.db.create_empty_data_source_file()
        self.db.create_missing_columns()
        self.database_directory=file_name[:file_name.rfind('/')]
        self.data_source_file=file_name[file_name.rfind('/')+1:]
        self.data_source_file_label.setText(os.path.join(self.database_directory,self.data_source_file))
        self.settings['database_directory']=self.database_directory
        self.settings['data_source']=self.data_source_file
        self.data_source_set=True
        self.datasource_menu_reload_samples()

    def button_clicked(self):

        if self.data_source_set==False:
            if self.sender().text()=="Create New Data\nSource (SQLite)":
                file_name = str(QtGui.QFileDialog.getSaveFileName(self.window,'Save file', self.database_directory))
                #make sure that the file always has .sqlite extension
                if file_name.rfind('.') != -1:
                    file_name=file_name[:file_name.rfind('.')]+'.sqlite'
                else:
                    file_name=file_name+'.sqlite'
                self.db=XChemDB.data_source(file_name)
                print '==> XCE: creating new data source'
                self.db.create_empty_data_source_file()
                self.db.create_missing_columns()
                if self.data_source_file=='':
                    self.database_directory=file_name[:file_name.rfind('/')]
                    self.data_source_file=file_name[file_name.rfind('/')+1:]
                    self.data_source_file_label.setText(os.path.join(self.database_directory,self.data_source_file))
                    self.settings['database_directory']=self.database_directory
                    self.settings['data_source']=self.data_source_file
                    self.data_source_set=True
            else:
                self.no_data_source_selected()
                pass

        # first find out which of the 'Run' or 'Status' buttons is sending
        for item in self.workflow_widget_dict:
            for widget in self.workflow_widget_dict[item]:
                if widget==self.sender():
                    # get index of item in self.workflow; Note this index should be the same as the index
                    # of the self.main_tab_widget which belongs to this task
                    task_index=self.workflow.index(item)
                    instruction =   str(self.workflow_widget_dict[item][0].currentText())
                    action =        str(self.sender().text())
                    if self.main_tab_widget.currentIndex()==task_index:
                        if self.explorer_active==0 and self.data_source_set==True:
                            if action=='Run':
                                self.prepare_and_run_task(instruction)
                            elif action=='Status':
                                self.get_status_of_workflow_milestone(instruction)
                                if os.path.exists(str(self.panddas_directory + '/pandda.done')):
                                    self.pandda_status = 'Finished!'
                                    self.pandda_status_label.setStyleSheet('color: green')
                                if os.path.exists(str(self.panddas_directory + '/pandda.running')):
                                    self.pandda_status = 'Running...'
                                    self.pandda_status_label.setStyleSheet('color: orange')
                                if os.path.exists(str(self.panddas_directory + '/pandda.errored')):
                                    self.pandda_status = 'Error encountered... please check the log files for pandda!'
                                    self.pandda_status_label.setStyleSheet('color: red')
                                self.pandda_status_label.setText(str('STATUS: ' + self.pandda_status))
                    else:
                        self.need_to_switch_main_tab(task_index)

    def get_status_of_workflow_milestone(self,instruction):
        # first update all tables
        self.datasource_menu_reload_samples()

        cluster_dict=XChemMain.get_jobs_running_on_cluster()

        self.update_log.insert('getting status updates...')

        self.status_bar.showMessage('please check terminal window for further information')

        self.update_log.insert('{0!s} samples are currently in database'.format(str(len(self.xtal_db_dict))))

        if 'DIMPLE' in instruction:
            XChemMain.print_cluster_status_message('dimple',cluster_dict,self.xce_logfile)

        elif 'Create CIF/PDB/PNG file' in instruction:
            XChemMain.print_acedrg_status(self.xce_logfile,self.xtal_db_dict)
            XChemMain.print_cluster_status_message('acedrg',cluster_dict,self.xce_logfile)

        elif instruction.startswith('Run xia2 on selected datasets'):
            XChemMain.print_cluster_status_message('xia2',cluster_dict,self.xce_logfile)

        elif 'pandda' in instruction.lower():
            XChemMain.print_cluster_status_message('pandda',cluster_dict,self.xce_logfile)

        elif 'coot' in instruction.lower():
            XChemMain.print_cluster_status_message('refmac',cluster_dict,self.xce_logfile)


    def prepare_and_run_task(self,instruction):

        if instruction=='Get New Results from Autoprocessing':
            self.check_for_new_autoprocessing_or_rescore(False)

        elif instruction=='Rescore Datasets':
            self.check_for_new_autoprocessing_or_rescore(True)

        elif instruction=="Read PKL file":
            summary = pickle.load( open( self.data_collection_summary_file, "rb") )
            self.create_widgets_for_autoprocessing_results_only(summary)

        elif instruction=='Run xia2 on selected datasets':
            self.run_xia2_on_selected_datasets(False)

        elif instruction=='Run xia2 on selected datasets - overwrite':
            self.run_xia2_on_selected_datasets(True)

        elif instruction=='Run DIMPLE on All Autoprocessing MTZ files':
            self.rerun_dimple_on_all_autoprocessing_files()

        elif instruction=='Run DIMPLE on selected MTZ files':
            self.run_dimple_on_selected_autoprocessing_file()

        elif instruction=='Remove selected DIMPLE PDB/MTZ files':
            self.remove_selected_dimple_files()

        elif instruction=='Create CIF/PDB/PNG file of ALL compounds':
            self.create_cif_pdb_png_files('ALL')

        elif instruction=='Create CIF/PDB/PNG file of NEW compounds':
            self.create_cif_pdb_png_files('NEW')

        elif instruction=='Create CIF/PDB/PNG file of SELECTED compounds':
            self.create_cif_pdb_png_files('SELECTED')

        elif instruction=='pandda.analyse':
            self.run_pandda_analyse('production_run')

        elif instruction=='pre-run for ground state model':
            self.run_pandda_analyse('pre_run')

        elif instruction=='pandda.inspect':
            self.run_pandda_inspect()

        elif instruction=='run pandda.inspect at home':
            self.run_pandda_inspect_at_home()

        elif instruction=='Export NEW PANDDA models':
            update_datasource_only=False
            which_models='new'
            self.run_pandda_export(update_datasource_only,which_models)

        elif instruction=='Export ALL PANDDA models':
            update_datasource_only=False
            which_models='all'
            self.run_pandda_export(update_datasource_only,which_models)

        elif instruction=='cluster datasets':
            self.cluster_datasets_for_pandda()

        elif instruction=='Update datasource with results from pandda.inspect':
            update_datasource_only=True
            which_models='all'
            self.run_pandda_export(update_datasource_only,which_models)

        elif instruction=='Show HTML summary':
            self.show_pandda_html_summary()

        elif instruction=='Event Map -> SF':
            self.convert_event_maps_to_SF()

        elif instruction=='check modelled ligands':
            self.compare_modelled_ligands_and_panddaTable()

        elif instruction.startswith("Open COOT") or instruction=='Build ground state model':
            if not self.coot_running:
                self.update_log.insert('starting coot...')
                if instruction=="Open COOT - new interface":
                    interface='new'
                elif instruction=="Open COOT for old PanDDA":
                    interface='panddaV1'
                elif instruction=='Build ground state model':
                    interface='reference'
                else:
                    interface='old'
                self.work_thread=XChemThread.start_COOT(self.settings,interface)
                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                self.work_thread.start()


        elif instruction=='Update Deposition Table':
            self.update_deposition_table()


    def set_new_reference_if_applicable(self):
        print 'hallo'
        reference_root=str(self.reference_file_selection_combobox.currentText())
        pg_ref=''
        ucVol_ref=0.0
        for reference in self.reference_file_list:
            print reference[0],reference_root
            if reference[0]==reference_root:
                pg_ref=reference[5]
                ucVol_ref=reference[4]
                break
        if ucVol_ref==0.0:
            self.update_log.insert('cannot set reference file since unit cell volume of reference pdb is 0!')
            return

        for xtal in self.initial_model_dimple_dict:
            reference_file_selection_combobox=self.initial_model_dimple_dict[xtal][1]
            db_dict=self.xtal_db_dict[xtal]
            pg_xtal=db_dict['DataProcessingPointGroup']
            ucVol_xtal=db_dict['DataProcessingUnitCellVolume']

            try:
                difference=math.fabs(1-(float(ucVol_xtal)/float(ucVol_ref)))*100
            except ValueError:
                self.update_log.insert(xtal+' -> cannot calculate unit cell volume difference')
                continue

            if pg_xtal==pg_ref and difference < self.allowed_unitcell_difference_percent:
                print xtal,pg_xtal,ucVol_xtal
                index = reference_file_selection_combobox.findText(reference_root, QtCore.Qt.MatchFixedString)
                reference_file_selection_combobox.setCurrentIndex(index)
                self.update_log.insert(xtal+' -> setting '+reference_root+' as input PDB file for DIMPLE')


    def check_status_create_png_of_soaked_compound(self):
        number_of_samples=0
        running=0
        timestamp_list=[]
        cif_file_generated=0
        for folder in glob.glob(os.path.join(self.initial_model_directory,'*','compound')):
            number_of_samples += 1
            if os.path.isfile(os.path.join(folder,'RESTRAINTS_IN_PROGRESS')):
                running += 1
                timestamp=datetime.fromtimestamp(os.path.getmtime(os.path.join(folder,'RESTRAINTS_IN_PROGRESS'))).strftime('%Y-%m-%d %H:%M:%S')
                timestamp_list.append(timestamp)
            for cif_file in glob.glob(os.path.join(folder,'*.cif')):
                if os.path.isfile(cif_file):
                    cif_file_generated += 1
        if timestamp_list != []:
            last_timestamp=max(timestamp_list)
        else:
            last_timestamp='n/a'
        message='Datasets: '+str(number_of_samples)+', jobs running: '+str(running)+', jobs finished: '+str(cif_file_generated)+', last job submmitted: '+str(last_timestamp)
        self.status_bar.showMessage(message)

    def check_for_new_autoprocessing_or_rescore(self,rescore_only):
        self.update_log.insert('checking for new data collection')
        start_thread=False
        if rescore_only:
            # first pop up a warning message as this will overwrite all user selections
            msgBox = QtGui.QMessageBox()
            msgBox.setText("*** WARNING ***\nThis will overwrite all your manual selections!\nDo you want to continue?")
            msgBox.addButton(QtGui.QPushButton('Yes'), QtGui.QMessageBox.YesRole)
            msgBox.addButton(QtGui.QPushButton('No'), QtGui.QMessageBox.RejectRole)
            reply = msgBox.exec_();
            if reply == 0:
                start_thread=True
            else:
                start_thread=False
        else:
            start_thread=True


        if start_thread:
            if self.target=='=== SELECT TARGET ===':
                msgBox = QtGui.QMessageBox()
                warning = ( '*** WARNING ***\n'
                            'You did not select a target!\n'
                            'In this case we will only parse the project directory!\n'
                            'Please note that this option is usually only useful in case you reprocessed your data.\n'
                            'Do you want to continue?'  )
                msgBox.setText(warning)
                msgBox.addButton(QtGui.QPushButton('Yes'), QtGui.QMessageBox.YesRole)
                msgBox.addButton(QtGui.QPushButton('No'), QtGui.QMessageBox.RejectRole)
                reply = msgBox.exec_();
                if reply == 0:
                    start_thread=True
                else:
                    start_thread=False
            else:
                start_thread=True

        if start_thread:
            self.work_thread=XChemThread.read_autoprocessing_results_from_disc(self.visit_list,
                                                                                self.target,
                                                                                self.reference_file_list,
                                                                                self.database_directory,
                                                                                self.data_collection_dict,
                                                                                self.preferences,
                                                                                self.data_collection_summary_file,
                                                                                self.initial_model_directory,
                                                                                rescore_only,
                                                                                self.acceptable_low_resolution_limit_for_data,
                                                                                os.path.join(self.database_directory,self.data_source_file),
                                                                                self.xce_logfile    )
            self.explorer_active=1
            self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.connect(self.work_thread, QtCore.SIGNAL("create_widgets_for_autoprocessing_results_only"),
                                                 self.create_widgets_for_autoprocessing_results_only)
            self.work_thread.start()

    def save_files_to_initial_model_folder(self):
        self.work_thread=XChemThread.save_autoprocessing_results_to_disc(self.dataset_outcome_dict,
                                                                             self.data_collection_table_dict,
                                                                             self.data_collection_column_three_dict,
                                                                             self.data_collection_dict,
                                                                             self.database_directory,self.data_source_file,
                                                                             self.initial_model_directory,
                                                                             self.preferences,
                                                                             self.data_collection_summary_file)
        self.explorer_active=1
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def run_pandda_analyse(self,run):
        pandda_params = {
                'data_dir':                             str(self.pandda_input_data_dir_entry.text()),
                'out_dir':                              str(self.pandda_output_data_dir_entry.text()),
                'submit_mode':                          str(self.pandda_submission_mode_selection_combobox.currentText()),
                'nproc':                                str(self.pandda_nproc_entry.text()),
                'min_build_datasets':                   str(self.pandda_min_build_dataset_entry.text()),
                'pdb_style':                            str(self.pandda_pdb_style_entry.text()),
                'mtz_style':                            str(self.pandda_mtz_style_entry.text()),
                'sort_event':                           str(self.pandda_sort_event_combobox.currentText()),
                'max_new_datasets':                     str(self.pandda_max_new_datasets_entry.text()),
                'grid_spacing':                         str(self.pandda_grid_spacing_entry.text()),
                'pandda_dir_structure':                 str(self.pandda_input_data_dir_entry.text()),
                'perform_diffraction_data_scaling':     str(self.wilson_checkbox.isChecked()),
                'filter_pdb':                           str(self.pandda_reference_file_selection_combobox.currentText()),
                'reference_dir':                        self.reference_directory,
                'appendix':                             '',
                'N_datasets':                           len(glob.glob(os.path.join(self.panddas_directory,'*'))),
                'write_mean_map':                       'interesting'
                        }

        if run=='pre_run':
            msgBox = QtGui.QMessageBox()
            msgBoxLayout = msgBox.layout()
            vbox = QtGui.QVBoxLayout()
            text = (    'The aim of the pre-run is NOT to identify bound ligands,\n'
                        'but to create mean ground state maps.  to pre-run will only comprise 100 datasets. The aim is not to identify\n'
                        'bound ligands, but to create the ground-state maps.\n'
                        'You can run m\n'
                        '- select "Build ground state model" \n'
                        '- calculate new maps with the improved reference structure\n'
                        '- run "pandda.analyse\n'    )
            vbox.addWidget(QtGui.QLabel(text))
            hbox=QtGui.QHBoxLayout()
            hbox.addWidget(QtGui.QLabel('appendix:'))
            appendix = QtGui.QLineEdit()
            appendix.setText('pre')
            appendix.setFixedWidth(200)
            hbox.addWidget(appendix)
            vbox.addLayout(hbox)

            msgBoxLayout.addLayout(vbox,0,0)
            msgBox.addButton(QtGui.QPushButton('Go'), QtGui.QMessageBox.YesRole)
            msgBox.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)
            reply = msgBox.exec_();
            if reply == 0:
                pandda_params['appendix']=str(appendix.text())
                pandda_params['max_new_datasets'] = '100'
                pandda_params['N_datasets'] = 100
                pandda_params['write_mean_map'] = 'all'
            else:
                return None

        self.update_log.insert('preparing pandda.analyse input script')
        self.work_thread=XChemPANDDA.run_pandda_analyse(pandda_params,self.xce_logfile,os.path.join(self.database_directory,self.data_source_file))
        self.work_thread=XChemPANDDA.run_pandda_analyse(pandda_params,self.xce_logfile,os.path.join(self.database_directory,self.data_source_file))
        self.connect(self.work_thread, QtCore.SIGNAL("datasource_menu_reload_samples"),self.datasource_menu_reload_samples)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()


    def cluster_datasets_for_pandda(self):

        pandda_params = {
                'out_dir':              str(self.pandda_output_data_dir_entry.text()),
                'pdb_style':            str(self.pandda_pdb_style_entry.text()),
                'mtz_style':            str(self.pandda_mtz_style_entry.text())
                        }
        self.update_log.insert('starting giant.cluster_mtzs_and_pdbs')
        self.work_thread=XChemPANDDA.giant_cluster_datasets(self.initial_model_directory,pandda_params,self.xce_logfile,os.path.join(self.database_directory,self.data_source_file),run_pandda_analyse)
        self.explorer_active=1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("datasource_menu_reload_samples"),self.datasource_menu_reload_samples)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def run_pandda_inspect(self):
        self.settings['panddas_directory']=str(self.pandda_output_data_dir_entry.text())
        print '==> XCE: starting pandda.inspect'
        self.work_thread=XChemThread.start_pandda_inspect(self.settings,self.xce_logfile)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def run_pandda_inspect_at_home(self):
        self.work_thread=XChemPANDDA.run_pandda_inspect_at_home(self.panddas_directory,self.xce_logfile)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)


    def convert_event_maps_to_SF(self):
        self.update_log.insert('converting all event maps in {0!s} to mtz files'.format(self.initial_model_directory))
        self.work_thread=XChemPANDDA.convert_all_event_maps_in_database(self.initial_model_directory,
                                                                        self.xce_logfile,
                                                                        os.path.join(self.database_directory,self.data_source_file))
        self.explorer_active=1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def compare_modelled_ligands_and_panddaTable(self):
        self.update_log.insert('checking agreement of ligands in refine.pdb and entries in panddaTable')
        self.work_thread=XChemPANDDA.check_number_of_modelled_ligands(self.initial_model_directory,
                                                                        self.xce_logfile,
                                                                        os.path.join(self.database_directory,self.data_source_file))
        self.explorer_active=1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.connect(self.work_thread, QtCore.SIGNAL("show_error_dict"), self.show_error_dict)
        self.work_thread.start()


    def run_pandda_export(self,update_datasource_only,which_models):
        self.settings['panddas_directory']=str(self.pandda_output_data_dir_entry.text())
        if update_datasource_only:
            self.update_log.insert('updating data source with results from pandda.inspect')
        else:
            self.update_log.insert('exporting PANDDA models, updating data source and launching inital refinement for new models')

        start_thread=False
        if which_models=='all':
            self.update_log.insert('exporting ALL models! *** WARNING *** This may overwrite previous refinements!!!')
            msgBox = QtGui.QMessageBox()
            msgBox.setText("*** WARNING ***\nThis will overwrite all your manual selections!\nDo you want to continue?")
            msgBox.addButton(QtGui.QPushButton('Yes'), QtGui.QMessageBox.YesRole)
            msgBox.addButton(QtGui.QPushButton('No'), QtGui.QMessageBox.RejectRole)
            reply = msgBox.exec_();
            if reply == 0:
                if update_datasource_only:
                    self.update_log.insert('will update panddaTable in database only')
                else:
                    self.update_log.insert('will export ALL models!')
                start_thread=True
            else:
                start_thread=False
        else:
            self.update_log.insert('exporting new models only')
            start_thread=True

        if start_thread:
            self.work_thread=XChemPANDDA.run_pandda_export(self.panddas_directory,os.path.join(self.database_directory,self.data_source_file),self.initial_model_directory,self.xce_logfile,update_datasource_only,which_models)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.work_thread.start()


    def show_pandda_html_summary(self):
        self.pandda_initial_html.load(QtCore.QUrl(self.pandda_initial_html_file))
        self.pandda_initial_html.show()
        self.pandda_analyse_html.load(QtCore.QUrl(self.pandda_analyse_html_file))
        self.pandda_analyse_html.show()
        self.add_map_html()
        self.pandda_inspect_html.load(QtCore.QUrl(self.pandda_inspect_html_file))
        self.pandda_inspect_html.show()



    def create_cif_pdb_png_files(self,todo):
        tmp=self.db.execute_statement("select CrystalName,CompoundCode,CompoundSmiles from mainTable where CrystalName is not '' and CompoundSmiles is not '' and CompoundSmiles is not NULL;")
        compound_list=[]
        for item in tmp:
            if str(item[1])=='' or str(item[1])=='NULL':
                compoundID='compound'
            else:
                compoundID=str(item[1])

            if todo == 'ALL':
                compound_list.append([str(item[0]),compoundID,str(item[2])])
            elif todo == 'NEW':
                if not os.path.isfile(os.path.join(self.initial_model_directory,str(item[0]),compoundID+'.cif')):
                    compound_list.append([str(item[0]),compoundID,str(item[2])])
            elif todo == 'SELECTED':
                if str(item[0]) in self.initial_model_dimple_dict:
                    if self.initial_model_dimple_dict[str(item[0])][0].isChecked():
                        compound_list.append([str(item[0]),compoundID,str(item[2])])

        if compound_list != []:
            self.update_log.insert('trying to create cif and pdb files for '+str(len(compound_list))+' compounds using ACEDRG...')
            if self.external_software['qsub']:
                self.update_log.insert('will try sending '+str(len(compound_list))+' jobs to your computer cluster!')
            elif self.external_software['qsub_array']:
                self.update_log.insert('will try sending '+str(len(compound_list))+' jobs as part of an ARRAY job to your computer cluster!')
            else:
                self.update_log.insert('apparently no cluster available, so will run '+str(len(compound_list))+' sequential jobs on one core of your local machine.')
                self.update_log.insert('this could take a while...')
            self.explorer_active=1
            self.work_thread=XChemThread.create_png_and_cif_of_compound(self.external_software,
                                                                        self.initial_model_directory,
                                                                        compound_list,
                                                                        self.database_directory,
                                                                        self.data_source_file,
                                                                        todo,
                                                                        self.ccp4_scratch_directory,
                                                                        self.xce_logfile,
                                                                        self.max_queue_jobs,
                                                                        self.restraints_program )
            self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.connect(self.work_thread, QtCore.SIGNAL("datasource_menu_reload_samples"),self.datasource_menu_reload_samples)
            self.work_thread.start()

    def update_deposition_table(self):
        # check if PanDDA models are ready for deposition

        depositChecks=XChemDeposit.update_deposition_table(os.path.join(self.database_directory,self.data_source_file))

        toDeposit,mismatch=depositChecks.PanDDA_models_to_deposit()

        if mismatch != {}:
            self.update_log.insert('The following samples contain ligand that are not ready for deposition:')
            for entry in mismatch:
                self.update_log.insert(entry[0]+' -> site: '+entry[1]+' @ '+entry[2]+' => '+entry[4])
            self.update_log.insert('You need to change this before you can continue!')
            return None

        for xtal in toDeposit:
            self.db.update_insert_depositTable(xtal,{})

    def show_html_summary_and_diffraction_image(self):
        for key in self.albula_button_dict:
            if self.albula_button_dict[key][0]==self.sender():
                print '==> XCE: showing html summary in firefox'
                self.show_html_summary_in_firefox(key)

    def need_to_switch_main_tab(self,task_index):
        msgBox = QtGui.QMessageBox()
        msgBox.setText("Need to switch main tab before you can launch this job")
        msgBox.addButton(QtGui.QPushButton('Yes'), QtGui.QMessageBox.YesRole)
        msgBox.addButton(QtGui.QPushButton('No'), QtGui.QMessageBox.RejectRole)
        reply = msgBox.exec_();
        if reply == 0:
            self.main_tab_widget.setCurrentIndex(task_index)

    def check_write_permissions_of_data_source(self):
        write_enabled=True
        if not os.access(os.path.join(self.database_directory,self.data_source_file),os.W_OK):
            QtGui.QMessageBox.warning(self.window, "Data Source Problem",
                                      '\nData Source is Read-Only\n',
                        QtGui.QMessageBox.Cancel, QtGui.QMessageBox.NoButton,
                        QtGui.QMessageBox.NoButton)
            write_enabled=False
        return write_enabled

    def no_data_source_selected(self):
        QtGui.QMessageBox.warning(self.window, "Data Source Problem",
                                      ('Please set or create a data source file\n')+
                                      ('Options:\n')+
                                      ('1. Use an existing file:\n')+
                                      ('- Settings -> Select Data Source File\n')+
#                                      ('- start XCE with command line argument (-d)\n')+
                                      ('2. Create a new file\n')+
                                      ('- Data Source -> Create New Data\nSource (SQLite)'),
                        QtGui.QMessageBox.Cancel, QtGui.QMessageBox.NoButton,
                        QtGui.QMessageBox.NoButton)

    def update_progress_bar(self,progress):
        self.progress_bar.setValue(progress)

    def update_status_bar(self,message):
        self.status_bar.showMessage(message)

    def thread_finished(self):
        self.explorer_active=0
        self.update_progress_bar(0)
        self.update_status_bar('idle')

    def show_error_dict(self,errorDict):
        text=''
        for key in errorDict:
            text+='{0!s}:\n'.format(key)
            for entry in errorDict[key]:
                text+='  - '+entry+'\n'
        msgBox = QtGui.QMessageBox()
        msgBox.setText(text)
        msgBox.exec_()

    def create_widgets_for_autoprocessing_results_only(self,data_dict):
        self.status_bar.showMessage('Building details table for data processing results')
        self.data_collection_dict=data_dict

        column_name = [ 'Program',
                        'Resolution\nOverall',
                        'Resolution\n[Mn<I/sig(I)> = 1.5]',
                        'DataProcessing\nSpaceGroup',
                        'Mn<I/sig(I)>\nHigh',
                        'Rmerge\nLow',
                        'Completeness\nOverall',
                        'DataProcessing\nUnitCell',
                        'DataProcessing\nRfree',
                        'DataProcessing\nScore'     ]

        # need to do this because db_dict keys are SQLite column names
        diffraction_data_column_name=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file)).translate_xce_column_list_to_sqlite(column_name)

        for xtal in sorted(self.data_collection_dict):
            if os.path.isfile(os.path.join(self.initial_model_directory,xtal,xtal+'.mtz')):
                mtz_already_in_inital_model_directory=True

            # column 2: data collection date
            # this one should always be there; it may need updating in case another run appears
            # first find latest run
            tmp=[]
            for entry in self.data_collection_dict[xtal]:
                if entry[0]=='image':
                    tmp.append( [entry[3],datetime.strptime(entry[3], '%Y-%m-%d %H:%M:%S')])
            latest_run=max(tmp,key=lambda x: x[1])[0]

            # first check if it does already exist
            if xtal not in self.data_collection_column_three_dict:
                # geneerate all the widgets which can later be appended and add them to the dictionary
                data_collection_table=QtGui.QTableWidget()      # table with data processing results for each pipeline
                selection_changed_by_user=False
                self.data_collection_column_three_dict[xtal]=[data_collection_table,selection_changed_by_user]
                xtal_in_table=True
            else:
                data_collection_table =     self.data_collection_column_three_dict[xtal][0]
                selection_changed_by_user = self.data_collection_column_three_dict[xtal][1]

            data_collection_table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
            data_collection_table.setColumnCount(len(column_name))
            font = QtGui.QFont()
            font.setPointSize(8)
            data_collection_table.setFont(font)
            data_collection_table.setHorizontalHeaderLabels(column_name)
            data_collection_table.horizontalHeader().setFont(font)
            data_collection_table.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)

            #############################################################################
            # crystal images
            # first check there are new images that are not displayed yet; i.e. they are not in the self.data_collection_image_dict
            if xtal not in self.data_collection_image_dict:
                # OK this is the first time
                self.data_collection_image_dict[xtal]=[]

            # sort crystal images by timestamp
            # reminder: ['image',visit,run,timestamp,image_list,diffraction_image,run_number]
            # a) get only image entries from self.data_collection_dict
            tmp=[]
            for entry in self.data_collection_dict[xtal]:
                if entry[0]=='image':
                    tmp.append(entry)

            # b) sort by the previously assigned run number
            #    note: entry[6]==run_number
            for entry in sorted(tmp,key=lambda x: x[6]):
                run_number=entry[6]
                images_already_in_table=False
                for image in self.data_collection_image_dict[xtal]:
                    if run_number==image[0]:
                        images_already_in_table=True
                        break
                if not images_already_in_table:
                # not if there is a run, but images are for whatever reason not present in self.data_collection_dict
                # then use image not available from $XChemExplorer_DIR/image/IMAGE_NOT_AVAILABLE.png
                # not sure how to do this at the moment; it will probably trigger an error that I can catch
                    self.data_collection_image_dict[xtal].append([entry[6],entry[1],entry[2],entry[3],entry[5]])

            #############################################################################
            # initialize dataset_outcome_dict for xtal
            if xtal not in self.dataset_outcome_dict:
                self.dataset_outcome_dict[xtal]=[]
                # dataset outcome buttons

            #############################################################################
            # table for data processing results
            # check if results from particular pipeline are already in table;
            # not really looking at the table here, but compare it to self.data_collection_table_dict
            row_position=data_collection_table.rowCount()
            if not xtal in self.data_collection_table_dict:
                self.data_collection_table_dict[xtal]=[]
            # reminder: ['logfile',visit,run,timestamp,autoproc,file_name,aimless_results,<aimless_index>,False]
            logfile_list=[]
            for entry in self.data_collection_dict[xtal]:
                if entry[0]=='logfile':
                    logfile_list.append(entry)
            for entry in sorted(logfile_list,key=lambda x: x[7]):               # sort by aimless_index and so make sure
                entry_already_in_table=False                                    # that aimless_index == row
                for logfile in self.data_collection_table_dict[xtal]:
                    if entry[1]==logfile[1] and entry[2]==logfile[2] and entry[3]==logfile[3] and entry[4]==logfile[4]:
                        entry_already_in_table=True
                        # might have to update Rfree column
                        for column,header in enumerate(diffraction_data_column_name):
                            if header=='DataProcessing\nRfree':
                                # entry[7]==aimless_index, i.e. row number
                                cell_text=QtGui.QTableWidgetItem()
                                cell_text.setText(str( db_dict[ header[1] ]  ))
                                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                                data_collection_table.setItem(entry[7], column, cell_text)
                                break
                        break
                if not entry_already_in_table:
                    data_collection_table.insertRow(row_position)
                    db_dict=entry[6]
                    for column,header in enumerate(diffraction_data_column_name):
                        cell_text=QtGui.QTableWidgetItem()
                        try:
                            cell_text.setText(str( db_dict[ header[1] ]  ))
                        except KeyError:
                            # this may happen if not score exists
                            cell_text.setText('0')
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        data_collection_table.setItem(row_position, column, cell_text)
                    data_collection_table.setRowHeight(row_position,20)
                    row_position+=1

                self.data_collection_table_dict[xtal].append(['logfile',entry[1],entry[2],entry[3],entry[4]])   # 'logfile' is just added to have
                                                                                                                # same index numbers between lists
            data_collection_table.cellClicked.connect(self.user_update_selected_autoproc_data_collection_summary_table)

            # select best resolution file + set data collection outcome
            # the assumption is that index in data_collection_dict and row number are identical
            # the assumption for data collection outcome is that as long as a logfile is found, it's a success
            logfile_found=False
            for entry in self.data_collection_dict[xtal]:
                if entry[0]=='logfile':
                    index=entry[7]
                    best_file=entry[8]
                    logfile_found=True
                    if best_file:
#                        # we change the selection only if the user did not touch it, assuming that he/she knows best
#                        if not selection_changed_by_user:
                        data_collection_table.selectRow(index)

        self.populate_data_collection_summary_table()

    def find_suitable_reference_file(self,db_dict):
        reference_file=[]
        dummy=['...', '', '', '', 0, '0']
        reference_file.append([dummy,999])
#        self.status_bar.showMessage('checking: '+str(os.path.join(db_dict['DataProcessingPathToMTZfile'],db_dict['DataProcessingMTZfileName'])))
        suitable_reference=[]
        for reference in self.reference_file_list:
            # first we need one in the same pointgroup
            if reference[5]==db_dict['DataProcessingPointGroup']:
                try:
                    difference=math.fabs(1-(float(db_dict['DataProcessingUnitCellVolume'])/float(reference[4])))*100
                    reference_file.append([reference,difference])
                except ValueError:
                    continue
        return reference_file


    def create_initial_model_table(self):
        column_name=self.db.translate_xce_column_list_to_sqlite(self.inital_model_column_list)

        for xtal in sorted(self.xtal_db_dict):
            new_xtal=False
            db_dict=self.xtal_db_dict[xtal]
            if str(db_dict['DataCollectionOutcome']).lower().startswith('success'):
                reference_file=self.find_suitable_reference_file(db_dict)
                smallest_uc_difference=min(reference_file,key=lambda x: x[1])
                row=self.initial_model_table.rowCount()
                if xtal not in self.initial_model_dimple_dict:
                    self.initial_model_table.insertRow(row)
                    current_row=row
                    new_xtal=True
                else:
                    for table_row in range(row):
                        if self.initial_model_table.item(table_row,0).text() == xtal:
                            current_row=table_row
                            break
                for column,header in enumerate(column_name):
                    if header[0]=='Sample ID':
                        cell_text=QtGui.QTableWidgetItem()
                        cell_text.setText(str(xtal))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.initial_model_table.setItem(current_row, column, cell_text)
                    elif header[0]=='Select':
                        if new_xtal:
                            run_dimple = QtGui.QCheckBox()
                            run_dimple.toggle()
                            self.initial_model_table.setCellWidget(current_row, column, run_dimple)
                            run_dimple.setChecked(False)
                    elif header[0]=='Reference\nSpaceGroup':
                        cell_text=QtGui.QTableWidgetItem()
                        cell_text.setText(str( smallest_uc_difference[0][1]  ))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.initial_model_table.setItem(current_row, column, cell_text)
                    elif header[0]=='Difference\nUC Volume (%)':
                        cell_text=QtGui.QTableWidgetItem()
                        smallest_uc_difference=min(reference_file,key=lambda x: x[1])
                        cell_text.setText(str( round(float(smallest_uc_difference[1]),1)  ))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.initial_model_table.setItem(current_row, column, cell_text)
                    elif header[0]=='Reference File':
                        if new_xtal:
                            reference_file_selection_combobox = QtGui.QComboBox()
                            self.populate_reference_combobox(reference_file_selection_combobox)
                            if float(smallest_uc_difference[1]) < self.allowed_unitcell_difference_percent:
                                index = reference_file_selection_combobox.findText(str(smallest_uc_difference[0][0]), QtCore.Qt.MatchFixedString)
                                reference_file_selection_combobox.setCurrentIndex(index)
                            else:
                                reference_file_selection_combobox.setCurrentIndex(0)
                            self.initial_model_table.setCellWidget(current_row, column, reference_file_selection_combobox)
                        else:
                            reference_file_selection_combobox=self.initial_model_dimple_dict[xtal][1]
                            self.populate_reference_combobox(reference_file_selection_combobox)
                            if float(smallest_uc_difference[1]) < self.allowed_unitcell_difference_percent:
                                index = reference_file_selection_combobox.findText(str(smallest_uc_difference[0][0]), QtCore.Qt.MatchFixedString)
                                reference_file_selection_combobox.setCurrentIndex(index)
                            else:
                                reference_file_selection_combobox.setCurrentIndex(0)
                    else:
                        cell_text=QtGui.QTableWidgetItem()
                        cell_text.setText(str( db_dict[ header[1] ]  ))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        if header[0]=='Dimple\nStatus':
                            if str( db_dict[ header[1] ]  ) == 'running':
                                cell_text.setBackground(QtGui.QColor(100,230,150))
                            elif str( db_dict[ header[1] ]  ) == 'pending':
                                cell_text.setBackground(QtGui.QColor(20,100,230))
                            elif str( db_dict[ header[1] ]  ) == 'started':
                                cell_text.setBackground(QtGui.QColor(230,240,110))
                            elif str( db_dict[ header[1] ]  ) == 'finished':
                                cell_text.setBackground(QtGui.QColor(255,255,255))
                        if header[0]=='Compound\nStatus':
                            if str( db_dict[ header[1] ]  ) == 'running':
                                cell_text.setBackground(QtGui.QColor(100,230,150))
                            elif str( db_dict[ header[1] ]  ) == 'pending':
                                cell_text.setBackground(QtGui.QColor(20,100,230))
                            elif str( db_dict[ header[1] ]  ) == 'started':
                                cell_text.setBackground(QtGui.QColor(230,240,110))
                            elif str( db_dict[ header[1] ]  ) == 'restraints generated':
                                cell_text.setBackground(QtGui.QColor(255,255,255))
                            elif str( db_dict[ header[1] ]  ) == 'restraints failed':
                                cell_text.setBackground(QtGui.QColor(255,0,0))
                            elif str( db_dict[ header[1] ]  ) == 'missing smiles':
                                cell_text.setBackground(QtGui.QColor(240,150,20))
                        self.initial_model_table.setItem(current_row, column, cell_text)
            if new_xtal:
                self.initial_model_dimple_dict[xtal]=[run_dimple,reference_file_selection_combobox]

    def preferences_data_to_copy_combobox_changed(self,i):
        text = str(self.preferences_data_to_copy_combobox.currentText())
        for item in self.preferences_data_to_copy:
            if item[0] == text:
                self.preferences['processed_data_to_copy']=item[1]
                break

    def preferences_selection_mechanism_combobox_changed(self,i):
        text = str(self.preferences_selection_mechanism_combobox.currentText())
        self.preferences['dataset_selection_mechanism']=text

    def preferences_restraints_generation_combobox_changed(self):
        text = str(self.preferences_restraints_generation_combobox.currentText())
        self.restraints_program=text
        self.update_log.insert('will use {0!s} for generation of ligand coordinates and restraints'.format(text))

    def refinement_outcome_combobox_changed(self):
        for xtal in self.summary_table_dict:
            if self.sender() == self.summary_table_dict[xtal]:
                db_dict={}
                db_dict['RefinementOutcome']=str(self.sender().currentText())
                self.db.create_or_remove_missing_records_in_depositTable(self.xce_logfile,xtal,'ligand_bound',db_dict)


    def get_reference_file_list(self,reference_root):
        # check available reference files
        reference_file_list=[]
        dummy=['...', '', '', '', 0, '0']
        reference_file_list.append(dummy)
        if os.path.isfile(os.path.join(self.reference_directory,reference_root+'.pdb')):
            pdb_reference=parse().PDBheader(os.path.join(self.reference_directory,reference_root+'.pdb'))
            spg_reference=pdb_reference['SpaceGroup']
            unitcell_reference=pdb_reference['UnitCell']
            lattice_reference=pdb_reference['Lattice']
            unitcell_volume_reference=pdb_reference['UnitCellVolume']
            pointgroup_reference=pdb_reference['PointGroup']
            reference_file_list.append([reference_root,
                                        spg_reference,
                                        unitcell_reference,
                                        lattice_reference,
                                        unitcell_volume_reference,
                                        pointgroup_reference])
        else:
            for files in glob.glob(self.reference_directory+'/*'):
                if files.endswith('.pdb'):
                    reference_root=files[files.rfind('/')+1:files.rfind('.')]
                    if os.path.isfile(os.path.join(self.reference_directory,reference_root+'.pdb')):
                        pdb_reference=parse().PDBheader(os.path.join(self.reference_directory,reference_root+'.pdb'))
                        spg_reference=pdb_reference['SpaceGroup']
                        unitcell_reference=pdb_reference['UnitCell']
                        lattice_reference=pdb_reference['Lattice']
                        unitcell_volume_reference=pdb_reference['UnitCellVolume']
                        pointgroup_reference=pdb_reference['PointGroup']
                        reference_file_list.append([reference_root,
                                                    spg_reference,
                                                    unitcell_reference,
                                                    lattice_reference,
                                                    unitcell_volume_reference,
                                                    pointgroup_reference])
        for n,file in enumerate(reference_file_list):
            self.update_log.insert('reference file {0!s}: {1!s}'.format(n, file))
        return reference_file_list

    def dataset_outcome_combobox_change_outcome(self,text):
        outcome=str(text)
        xtal=''
        for key in self.dataset_outcome_combobox_dict:
            if self.dataset_outcome_combobox_dict[key]==self.sender():
                xtal=key
                self.update_log.insert('user changed data collection outcome of {0!s} to {1!s}'.format(xtal, outcome))
                break
        self.dataset_outcome_dict[xtal]=outcome
        if xtal != '':
            # need to also update if not yet done
            user_already_changed_selection=False
            for n,entry in enumerate(self.data_collection_dict[xtal]):
                if entry[0]=='user_changed_selection':
                    user_already_changed_selection=True
                if entry[0]=='logfile':
                    db_dict=entry[6]
                    db_dict['DataCollectionOutcome']=outcome
                    entry[6]=db_dict
                    self.data_collection_dict[xtal][n]=entry
            if not user_already_changed_selection:
                self.data_collection_dict[xtal].append(['user_changed_selection'])
            # finally need to update outcome field in data source accordingly
            self.update_log.insert('updating dataset outcome in datasource for {0!s}'.format(xtal))
            update_dict={}
            update_dict['DataCollectionOutcome']=outcome
            self.db.update_insert_data_source(xtal,update_dict)

    def set_run_dimple_flag(self,state):
        if state == QtCore.Qt.Checked:
            for key in self.initial_model_dimple_dict:
                self.initial_model_dimple_dict[key][0].setChecked(True)
        else:
            for key in self.initial_model_dimple_dict:
                self.initial_model_dimple_dict[key][0].setChecked(False)

    def show_data_collection_details(self,state):
        # first remove currently displayed widget
        if self.data_collection_details_currently_on_display is not None:
            self.data_collection_details_currently_on_display.hide()
            self.data_collection_details_currently_on_display=None

        tmp=[]
        allRows = self.data_collection_summary_table.rowCount()
        for table_row in range(allRows):
            tmp.append([self.data_collection_summary_table.item(table_row,0).text(),table_row])

        for key in self.data_collection_summary_dict:
            if self.data_collection_summary_dict[key][3]==self.sender():
                if self.sender().isChecked():
                    for item in tmp:
                        if item[0]==key:
                            self.data_collection_summary_table.selectRow(item[1])
                    self.data_collection_details_currently_on_display=self.data_collection_column_three_dict[key][0]
                    self.data_collection_summarys_vbox_for_details.addWidget(self.data_collection_details_currently_on_display)
                    self.data_collection_details_currently_on_display.show()
            else:
                # un-check all other ones
                self.data_collection_summary_dict[key][3].setChecked(False)

    def continously_check_for_new_data_collection(self,state):
        self.timer_to_check_for_new_data_collection.timeout.connect(lambda: self.check_for_new_autoprocessing_or_rescore(False))
        if state == QtCore.Qt.Checked:
            print '==> XCE: checking automatically every 120s for new data collection'
            self.timer_to_check_for_new_data_collection.start(120000)
        else:
            print '==> XCE: stopped checking for new data collections'
            self.timer_to_check_for_new_data_collection.stop()

    def populate_data_collection_summary_table(self):
        self.status_bar.showMessage('Building summary table for data processing results; be patient this may take a while')
        row = self.data_collection_summary_table.rowCount()
        column_name=self.db.translate_xce_column_list_to_sqlite(self.data_collection_summary_column_name)

        pinList = self.db.execute_statement("Select CrystalName,PinBarcode,DataCollectionPinBarcode from mainTable where CrystalName is not ''")
        pinDict={}
        for item in pinList:
            pinDict[str(item[0])]=[str(item[1]),str(item[2])]

        for xtal in sorted(self.data_collection_dict):
            new_xtal=False
            if xtal not in self.data_collection_summary_dict:
                row = self.data_collection_summary_table.rowCount()
                self.data_collection_summary_table.insertRow(row)
                self.data_collection_summary_dict[xtal]=[]
                # self.data_collection_summary_dict[xtal]=[outcome,db_dict,image,diffraction_image]
                new_xtal=True

            # check for dataset outcome
            outcome=''
            logfile_found=False
            too_low_resolution=True
            db_dict={}
            for entry in self.data_collection_dict[xtal]:
                if entry[0]=='logfile':
                    logfile_found=True
                    if entry[8]:    # if this was auto-selected best resolution file
                        db_dict=entry[6]
                        try:
                            if float(db_dict['DataProcessingResolutionHigh']) <= float(self.acceptable_low_resolution_limit_for_data):
                                too_low_resolution=False
                        except ValueError:
                            pass

            try:
                outcome=str(self.db.get_value_from_field(xtal,'DataCollectionOutcome')[0])
            except TypeError:
                outcome='Failed - unknown'
                self.update_log.insert('cannot find DataCollectionOutcome for {0!s}'.format(xtal))
            self.dataset_outcome_dict[xtal]=outcome

            # find latest run for crystal and diffraction images
            tmp=[]
            for entry in self.data_collection_dict[xtal]:
                if entry[0]=='image':
                    tmp.append( [entry,datetime.strptime(entry[3], '%Y-%m-%d %H:%M:%S')])
            latest_run=max(tmp,key=lambda x: x[1])[0]


            new_run_for_exisiting_crystal_or_new_sample=True
            if new_xtal:
                self.data_collection_summary_dict[xtal]=[outcome,db_dict,latest_run]
            else:
                # check if newer run appeared
                old_run_timestamp=self.data_collection_summary_dict[xtal][2][3]
                new_run_timestamp=latest_run[3]
                if old_run_timestamp == new_run_timestamp:
                    new_run_for_exisiting_crystal_or_new_sample=False
                else:
                    checkbox_for_details=self.data_collection_summary_dict[xtal][3]
                    self.data_collection_summary_dict[xtal]=[outcome,db_dict,latest_run,checkbox_for_details]

            if new_xtal:
                current_row=row
            else:
                allRows = self.data_collection_summary_table.rowCount()
                for table_row in range(allRows):
                    if self.data_collection_summary_table.item(table_row,0).text() == xtal:
                        current_row=table_row
                        break

            image_number=0
            for column,header in enumerate(column_name):
                if header[0]=='Sample ID':
                    cell_text=QtGui.QTableWidgetItem()
                    cell_text.setText(str(xtal))
                    cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                    self.data_collection_summary_table.setItem(current_row, column, cell_text)
                elif header[0]=='DataCollection\nOutcome':
                    if new_xtal:
                        dataset_outcome_combobox = QtGui.QComboBox()
                        for outcomeItem in self.dataset_outcome:
                            dataset_outcome_combobox.addItem(outcomeItem)
                        self.data_collection_summary_table.setCellWidget(current_row, column, dataset_outcome_combobox)
                        dataset_outcome_combobox.activated[str].connect(self.dataset_outcome_combobox_change_outcome)
                        self.dataset_outcome_combobox_dict[xtal]=dataset_outcome_combobox
                    index = self.dataset_outcome_combobox_dict[xtal].findText(str(outcome), QtCore.Qt.MatchFixedString)
                    self.dataset_outcome_combobox_dict[xtal].setCurrentIndex(index)
                    continue

                elif header[0].startswith('img'):
                    if new_run_for_exisiting_crystal_or_new_sample:
                        img=latest_run[4]
                        pixmap = QtGui.QPixmap()
                        # can do this (img[image_number][1]) because made sure in the threading module
                        # that there are always exactly 5 images in there
                        pixmap.loadFromData(base64.b64decode(img[image_number][1]))
                        image = QtGui.QLabel()
                        image.resize(128,80)
                        image.setPixmap(pixmap.scaled(image.size(), QtCore.Qt.KeepAspectRatio))
                        self.data_collection_summary_table.setCellWidget(current_row, column, image)
                        image_number+=1

                elif header[0].startswith('Show Diffraction\nImage'):
                    if new_run_for_exisiting_crystal_or_new_sample:
                        diffraction_image=latest_run[5]
                        diffraction_image_name=diffraction_image[diffraction_image.rfind('/')+1:]
                        try:    # need to try because older pkl file may not have this item in list
                            html_summary=latest_run[7]
                        except IndexError:
                            html_summary=''
                        if new_xtal:
                            start_albula_button=QtGui.QPushButton('Show: \n'+diffraction_image_name)
                            start_albula_button.clicked.connect(self.show_html_summary_and_diffraction_image)
                            self.albula_button_dict[xtal]=[start_albula_button,diffraction_image,html_summary]
                            self.data_collection_summary_table.setCellWidget(current_row,column,start_albula_button)
                        else:
                            self.albula_button_dict[xtal][1]=diffraction_image
                elif header[0].startswith('Show\nDetails'):
                    if new_xtal:
                        show_data_collection_details_checkbox=QtGui.QCheckBox()
                        show_data_collection_details_checkbox.toggle()
                        show_data_collection_details_checkbox.setChecked(False)
                        show_data_collection_details_checkbox.stateChanged.connect(self.show_data_collection_details)
                        self.data_collection_summary_table.setCellWidget(current_row,column,show_data_collection_details_checkbox)
                        self.data_collection_summary_dict[xtal].append(show_data_collection_details_checkbox)
                elif header[0].startswith('SoakDB\nBarcode'):
                    if new_xtal:
                        cell_text=QtGui.QTableWidgetItem()
                        if xtal in pinDict:
                            cell_text.setText(str(pinDict[xtal][0]))
                            if pinDict[xtal][0] == 'NULL' or pinDict[xtal][1] == 'NULL':
                                cell_text.setBackground(QtGui.QColor(255,215,0))
                            elif pinDict[xtal][0] != pinDict[xtal][1]:
                                cell_text.setBackground(QtGui.QColor(255,0,0))
                        else:
                            cell_text.setText('')
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.data_collection_summary_table.setItem(current_row, column, cell_text)

                elif header[0].startswith('GDA\nBarcode'):
                    if new_xtal:
                        cell_text=QtGui.QTableWidgetItem()
                        if xtal in pinDict:
                            cell_text.setText(str(pinDict[xtal][1]))
                            if pinDict[xtal][0] == 'NULL' or pinDict[xtal][1] == 'NULL':
                                cell_text.setBackground(QtGui.QColor(255,215,0))
                            elif pinDict[xtal][0] != pinDict[xtal][1]:
                                cell_text.setBackground(QtGui.QColor(255,0,0))
                        else:
                            cell_text.setText('')
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.data_collection_summary_table.setItem(current_row, column, cell_text)



                else:
                    cell_text=QtGui.QTableWidgetItem()
                    # in case data collection failed for whatever reason
                    if logfile_found:
                        try:
                            cell_text.setText(str( db_dict[ header[1] ]  ))
                        except KeyError:                # older pkl files may not have all the columns
                            cell_text.setText('n/a')
                    else:
                        if header[0].startswith('Resolution\n[Mn<I/sig(I)> = 1.5]'):
                            cell_text.setText('999')
                        elif header[0].startswith('DataProcessing\nRfree'):
                            cell_text.setText('999')
                        elif header[0].startswith('Rmerge\nLow'):
                            cell_text.setText('999')
                        else:
                            cell_text.setText('')
                    cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                    self.data_collection_summary_table.setItem(current_row, column, cell_text)

            row += 1

        self.data_collection_summary_table.resizeRowsToContents()
        self.data_collection_summary_table.resizeColumnsToContents()

        self.status_bar.showMessage('updating Overview table')

        self.status_bar.showMessage('idle')

        self.save_files_to_initial_model_folder()

        self.update_dewar_configuration_tab()

    def update_dewar_configuration_tab(self):

        # sample status and color code:
        # 1 -   green:    collected and 'good' data
        # 2 -   orange:   collected and some data
        # 3 -   red:      collected, but no logfile
        # 4 -   blue:     sample in dewar, but not yet collected
        # 5 -   grey:     no sample in respective dewar position
        # 6 -   yellow:   flagged for re-collection

        # first find out what is currently in the dewar

        occupied_positions=[]
        for puck_position in self.dewar_sample_configuration_dict:
            sample=self.dewar_sample_configuration_dict[puck_position]
            if sample==[]:
                self.dewar_configuration_dict[puck_position].setStyleSheet("background-color: gray")
                continue
            elif sample not in self.data_collection_dict:
                self.dewar_configuration_dict[puck_position].setStyleSheet("background-color: blue")
                # color and name respective button
            else:
                logfile_found=False
                for entry in self.data_collection_dict[sample]:
                    if entry[0]=='logfile':
                        logfile_found=True
                        if entry[8]:    # if this was auto-selected best resolution file
                            db_dict=entry[6]
                            resolution_high=db_dict['DataProcessingResolutionHigh']
                if not logfile_found:
                    resolution_high='no logfile'
                self.dewar_configuration_dict[puck_position].setText(sample+'\n'+resolution_high+'A')
                outcome=str(self.dataset_outcome_combobox_dict[sample].currentText())
                if outcome=="success":
                    self.dewar_configuration_dict[puck_position].setStyleSheet("background-color: green")
                elif outcome=="Failed - low resolution":
                    self.dewar_configuration_dict[puck_position].setStyleSheet("background-color: orange")
                else:
                    self.dewar_configuration_dict[puck_position].setStyleSheet("background-color: red")
                self.dewar_configuration_dict[puck_position].setStyleSheet("font-size:7px;border-width: 0px")




    def update_outcome_data_collection_summary_table(self,sample,outcome):
        rows_in_table=self.data_collection_summary_table.rowCount()
        for row in range(rows_in_table):
            if self.data_collection_summary_table.item(row,0).text()==sample:
                cell_text=QtGui.QTableWidgetItem()
                cell_text.setText(outcome)
                self.data_collection_summary_table.setItem(row, 3, cell_text)

    def user_update_selected_autoproc_data_collection_summary_table(self):
        for key in self.data_collection_column_three_dict:
            if self.data_collection_column_three_dict[key][0]==self.sender():
                dbTmp=self.xtal_db_dict[key]
                stage=dbTmp['RefinementOutcome'].split()[0]
                print '===>',key,stage
                if int(stage) > 2:
                    msgBox = QtGui.QMessageBox()
                    msgBox.setText("*** WARNING ***\n%s is currently %s\nIt will disappear from the Refinement table,\nwhen you refresh it next time.\nDo you want to continue?" %(key,dbTmp['RefinementOutcome']))
                    msgBox.addButton(QtGui.QPushButton('No'), QtGui.QMessageBox.YesRole)
                    msgBox.addButton(QtGui.QPushButton('Yes'), QtGui.QMessageBox.RejectRole)
                    reply = msgBox.exec_();
                    if reply == 0:
                        self.update_log.insert('will not change data processing selection')
                        # restore previous selection
                        for n,entry in enumerate(self.data_collection_dict[key]):
                            print '==>',n
                            if entry[0]=='logfile':
                                if entry[8]==True:
                                    print '===> found:',n
                                    self.data_collection_column_three_dict[key][0].selectRow(n)
                        break

                indexes=self.sender().selectionModel().selectedRows()
                selected_processing_result=1000000
                for index in sorted(indexes):
                    selected_processing_result=index.row()
                # the user changed the selection, i.e. no automated selection will update it
                self.update_log.insert('user changed selection')
                self.data_collection_column_three_dict[key][1]=True
                # need to also update if not yet done
                user_already_changed_selection=False
                for n,entry in enumerate(self.data_collection_dict[key]):
                    if entry[0]=='user_changed_selection':
                        user_already_changed_selection=True
                    if entry[0]=='logfile':
                        db_dict=entry[6]
                        db_dict['DataProcessingAutoAssigned']='False'
                        if entry[7]==selected_processing_result:
                            db_dict_current=entry[6]
                            program=db_dict['DataProcessingProgram']
                            visit=db_dict['DataCollectionVisit']
                            run=db_dict['DataCollectionRun']
                            self.update_log.insert('user changed data processing files for {0!s} to visit={1!s}, run={2!s}, program={3!s}'.format(key, visit, run, program))
                            # update datasource
                            self.update_log.insert('updating datasource...')
                            self.update_data_source(key,db_dict)
                            entry[8]=True
                        else:
                            entry[8]=False

                        entry[6]=db_dict
                        self.data_collection_dict[key][n]=entry
                if not user_already_changed_selection:
                    self.data_collection_dict[key].append(['user_changed_selection'])
                XChemMain.change_links_to_selected_data_collection_outcome(key,self.data_collection_dict,
                                                                           self.data_collection_column_three_dict,
                                                                           self.dataset_outcome_dict,
                                                                           self.initial_model_directory,
                                                                           os.path.join(self.database_directory,self.data_source_file),
                                                                           self.xce_logfile)

                # update 'Datasets' table
                column_name=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file)).translate_xce_column_list_to_sqlite(self.data_collection_summary_column_name)
                rows_in_table=self.data_collection_summary_table.rowCount()
                for row in range(rows_in_table):
                    if self.data_collection_summary_table.item(row,0).text()==key:
                        for column,header in enumerate(column_name):
                            if header[0]=='Sample ID':
                                continue
                            elif header[0]=='DataCollection\nOutcome':
                                continue
                            elif header[0].startswith('img'):
                                continue
                            elif header[0].startswith('Show'):
                                continue
                            else:
                                cell_text=QtGui.QTableWidgetItem()
                                try:
                                    cell_text.setText(str( db_dict_current[ header[1] ]  ))
                                    cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                                    self.data_collection_summary_table.setItem(row, column, cell_text)
                                except KeyError:
                                    pass



    def update_selected_autoproc_data_collection_summary_table(self):
        for key in self.data_collection_column_three_dict:
            if self.data_collection_column_three_dict[key][0]==self.sender():
                sample=key
                break
        indexes=self.sender().selectionModel().selectedRows()
        for index in sorted(indexes):
            selected_processing_result=index.row()

        for n,entry in enumerate(self.data_collection_dict[sample]):
            if entry[0]=='logfile':
                if entry[7]==selected_processing_result:
                    db_dict=entry[6]
                    program=db_dict['DataProcessingProgram']
                    visit=db_dict['DataCollectionVisit']
                    run=db_dict['DataCollectionRun']
                    self.update_log.insert('user changed data processing files for {0!s} to visit={1!s}, run={2!s}, program={3!s}'.format(sample, visit, run, program))
                    # update datasource
                    self.update_log.insert('updating datasource...')
                    self.update_data_source(sample,db_dict)
                    entry[8]=True
                else:
                    entry[8]=False
                self.data_collection_dict[sample][n]=entry

        # update 'Datasets' table
        column_name=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file)).translate_xce_column_list_to_sqlite(self.data_collection_summary_column_name)
        rows_in_table=self.data_collection_summary_table.rowCount()
        for row in range(rows_in_table):
            if self.data_collection_summary_table.item(row,0).text()==sample:
                for column,header in enumerate(column_name):
                    if header[0]=='Sample ID':
                        continue
                    elif header[0]=='DataCollection\nOutcome':
                        continue
                    elif header[0].startswith('img'):
                        continue
                    elif header[0].startswith('Show'):
                        continue
                    else:
                        cell_text=QtGui.QTableWidgetItem()
                        cell_text.setText(str( db_dict[ header[1] ]  ))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.data_collection_summary_table.setItem(row, column, cell_text)



    def populate_data_source_table(self,header,data):
        self.mounted_crystal_table.setColumnCount(0)
        self.mounted_crystal_table.setColumnCount(len(self.data_source_columns_to_display))
        self.mounted_crystal_table.setRowCount(0)

        columns_to_show=self.get_columns_to_show(self.data_source_columns_to_display,header)
        n_rows=self.get_rows_with_sample_id_not_null(header,data)
        self.mounted_crystal_table.setRowCount(n_rows)
        sample_id_column=self.get_columns_to_show(['Sample ID'],header)

        x=0
        for row in data:
            if str(row[sample_id_column[0]]).lower() == 'none' or str(row[sample_id_column[0]]).replace(' ','') == '':
                continue        # do not show rows where sampleID is null
            for y,item in enumerate(columns_to_show):
                cell_text=QtGui.QTableWidgetItem()
                if row[item] is None:
                    cell_text.setText('')
                else:
                    cell_text.setText(str(row[item]))
                if self.data_source_columns_to_display[y]=='Sample ID':     # assumption is that column 0 is always sampleID
                    cell_text.setFlags(QtCore.Qt.ItemIsEnabled)             # and this field cannot be changed
                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                self.mounted_crystal_table.setItem(x, y, cell_text)
            x+=1
        self.mounted_crystal_table.setHorizontalHeaderLabels(self.data_source_columns_to_display)

    def populate_and_update_data_source_table(self):
        self.mounted_crystal_table.setColumnCount(len(self.data_source_columns_to_display))

        # first get a list of all the samples that are already in the table and which will be updated
        samples_in_table=[]
        current_row = self.mounted_crystal_table.rowCount()
        for row in range(current_row):
            sampleID=str(self.mounted_crystal_table.item(row,0).text())      # this must be the case
            samples_in_table.append(sampleID)

        columns_to_show=self.get_columns_to_show(self.data_source_columns_to_display)
        n_rows=self.get_rows_with_sample_id_not_null_from_datasource()
        sample_id_column=self.get_columns_to_show(['Sample ID'])


        for row in self.data:
            if str(row[sample_id_column[0]]).lower() == 'none' or str(row[sample_id_column[0]]).replace(' ','') == '':
                # do not show rows where sampleID is null
                continue
            else:
                if not str(row[sample_id_column[0]]) in samples_in_table:
                    # insert row, this is a new sample
                    x = self.mounted_crystal_table.rowCount()
                    self.mounted_crystal_table.insertRow(x)
                else:
                    # find row of this sample in data_source_table
                    for present_rows in range(self.mounted_crystal_table.rowCount()):
                        if str(row[sample_id_column[0]])==str(self.mounted_crystal_table.item(present_rows,0).text()):
                            x = present_rows
                            break
            for y,item in enumerate(columns_to_show):
                cell_text=QtGui.QTableWidgetItem()
                if row[item] is None:
                    cell_text.setText('')
                else:
                    cell_text.setText(str(row[item]))
                if self.data_source_columns_to_display[y]=='Sample ID':     # assumption is that column 0 is always sampleID
                    cell_text.setFlags(QtCore.Qt.ItemIsEnabled)             # and this field cannot be changed
                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                self.mounted_crystal_table.setItem(x, y, cell_text)
        self.mounted_crystal_table.setHorizontalHeaderLabels(self.data_source_columns_to_display)

    def populate_pandda_analyse_input_table(self):

        column_name=self.db.translate_xce_column_list_to_sqlite(self.pandda_column_name)
        for xtal in sorted(self.xtal_db_dict):
            new_xtal=False
            db_dict=self.xtal_db_dict[xtal]
            if os.path.isfile(db_dict['DimplePathToPDB']):
                row=self.pandda_analyse_data_table.rowCount()
                if xtal not in self.pandda_analyse_input_table_dict:
                    self.pandda_analyse_data_table.insertRow(row)
                    current_row=row
                    new_xtal=True
                else:
                    for table_row in range(row):
                        if self.pandda_analyse_data_table.item(table_row,0).text() == xtal:
                            current_row=table_row
                            break
                for column,header in enumerate(column_name):
                    if header[0]=='Sample ID':
                        cell_text=QtGui.QTableWidgetItem()
                        cell_text.setText(str(xtal))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.pandda_analyse_data_table.setItem(current_row, column, cell_text)
                    else:
                        cell_text=QtGui.QTableWidgetItem()
                        cell_text.setText(str( db_dict[ header[1] ]  ))
                        if header[0]=='PanDDA\nStatus':
                            if str( db_dict[ header[1] ]  ) == 'running':
                                cell_text.setBackground(QtGui.QColor(100,230,150))
                            elif str( db_dict[ header[1] ]  ) == 'pending':
                                cell_text.setBackground(QtGui.QColor(20,100,230))
                            elif str( db_dict[ header[1] ]  ) == 'started':
                                cell_text.setBackground(QtGui.QColor(230,240,110))
                            elif str( db_dict[ header[1] ]  ) == 'finished':
                                cell_text.setBackground(QtGui.QColor(255,255,255))
                            elif 'problem' in str( db_dict[ header[1] ]  ):
                                cell_text.setBackground(QtGui.QColor(255,0,0))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.pandda_analyse_data_table.setItem(current_row, column, cell_text)
            if new_xtal:
                self.pandda_analyse_input_table_dict[xtal]=[]

    def populate_and_update_refinement_table(self):

        panddaList = self.db.execute_statement("select CrystalName,PANDDA_site_index,PANDDA_site_name,RefinementOutcome from panddaTable where CrystalName is not '' and PANDDA_site_ligand_placed is 'True';")
        panddaDict={}
        for item in panddaList:
            if str(item[0]) not in panddaDict:
                panddaDict[str(item[0])]=[]
            panddaDict[str(item[0])].append([str(item[1]),str(item[2]),str(item[3])])

        column_name=self.db.translate_xce_column_list_to_sqlite(self.summary_column_name)
        for xtal in sorted(self.xtal_db_dict):
            new_xtal=False
            db_dict=self.xtal_db_dict[xtal]
            try:
                stage = int(str(db_dict['RefinementOutcome']).split()[0])
                refinementStage=db_dict['RefinementOutcome']
            except ValueError:
                stage = 0
            except IndexError:
                stage = 0

            if stage >= 3:
                row=self.summary_table.rowCount()
                if xtal not in self.summary_table_dict:
                    self.summary_table.insertRow(row)
                    current_row=row
                    new_xtal=True
                else:
                    for table_row in range(row):
                        if self.summary_table.item(table_row,0).text() == xtal:
                            current_row=table_row
                            break
                for column,header in enumerate(column_name):
                    if header[0]=='Sample ID':
                        cell_text=QtGui.QTableWidgetItem()
                        cell_text.setText(str(xtal))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.summary_table.setItem(current_row, column, cell_text)

                    elif header[0]=='Refinement\nOutcome':
                        if new_xtal:
                            refinement_outcome_combobox = QtGui.QComboBox()
                            self.populate_refinement_outcome_combobox(refinement_outcome_combobox)
                            self.summary_table.setCellWidget(current_row, column, refinement_outcome_combobox)
                        else:
                            refinement_outcome_combobox=self.summary_table_dict[xtal]
                        index = refinement_outcome_combobox.findText(refinementStage, QtCore.Qt.MatchFixedString)
                        refinement_outcome_combobox.setCurrentIndex(index)
                        refinement_outcome_combobox.currentIndexChanged.connect(self.refinement_outcome_combobox_changed)



                    elif header[0]=='PanDDA site details':
                        try:
                            panddaDict[xtal].insert(0,['Index','Name','Status'])
                            outerFrame=QtGui.QFrame()
                            outerFrame.setFrameShape(QtGui.QFrame.Box)
                            grid = QtGui.QGridLayout()
                            for y,entry in enumerate(panddaDict[xtal]):
                                for x,info in enumerate(entry):
                                    frame=QtGui.QFrame()
                                    frame.setFrameShape(QtGui.QFrame.Box)
                                    vbox=QtGui.QVBoxLayout()
                                    vbox.addWidget(QtGui.QLabel(str(entry[x])))
                                    frame.setLayout(vbox)
                                    grid.addWidget(frame,y,x)
                            outerFrame.setLayout(grid)
                            self.summary_table.setCellWidget(current_row, column, outerFrame)
                        except KeyError:
                            cell_text=QtGui.QTableWidgetItem()
                            cell_text.setText('*** N/A ***')
                            cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                            self.summary_table.setItem(current_row, column, cell_text)
                    else:
                        cell_text=QtGui.QTableWidgetItem()
                        cell_text.setText(str( db_dict[ header[1] ]  ))
                        if header[0]=='Refinement\nStatus':
                            if str( db_dict[ header[1] ]  ) == 'running':
                                cell_text.setBackground(QtGui.QColor(100,230,150))
                            elif str( db_dict[ header[1] ]  ) == 'pending':
                                cell_text.setBackground(QtGui.QColor(20,100,230))
                            elif str( db_dict[ header[1] ]  ) == 'started':
                                cell_text.setBackground(QtGui.QColor(230,240,110))
                            elif str( db_dict[ header[1] ]  ) == 'finished':
                                cell_text.setBackground(QtGui.QColor(255,255,255))
                            elif 'problem' in str( db_dict[ header[1] ]  ):
                                cell_text.setBackground(QtGui.QColor(255,0,0))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.summary_table.setItem(current_row, column, cell_text)
            if new_xtal:
                self.summary_table_dict[xtal]=refinement_outcome_combobox

        self.summary_table.resizeColumnsToContents()
        self.summary_table.resizeRowsToContents()



    def get_columns_to_show(self,column_list):
        # maybe I coded some garbage before, but I need to find out which column name in the
        # data source corresponds to the actually displayed column name in the table
        # reason being that the unique column ID for DB may not be nice to look at
        columns_to_show=[]
        for column in column_list:
            # first find out what the column name in the header is:
            column_name=''
            for name in self.all_columns_in_data_source:
                if column==name[1]:
                    column_name=name[0]
            for n,all_column in enumerate(self.header):
                if column_name==all_column:
                    columns_to_show.append(n)
                    break
        return columns_to_show


    def get_rows_with_sample_id_not_null_from_datasource(self):
        sample_id_column=self.get_columns_to_show(['Sample ID'])
        n_rows=0
        for row in self.data:
            if not str(row[sample_id_column[0]]).lower() != 'none' or not str(row[sample_id_column[0]]).replace(' ','') == '':
                n_rows+=1
        return n_rows

    def update_data_source(self,sample,db_dict):
        data_source=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))

    def quit_xce(self):
        # save pkl file
        if self.data_collection_dict != {}:
            if os.path.isfile(self.data_collection_summary_file):
                self.update_log.insert('saving results to PKL file')
                pickle.dump(self.data_collection_dict,open(self.data_collection_summary_file,'wb'))
        self.update_log.insert('quitting XCE... bye,bye!')
        QtGui.qApp.quit()

    def add_lead(self):
        # copy pandda_analyse_sites.csv to proasis directory for lead build
        os.system(str('cp ' + str(os.path.join(self.panddas_directory, 'analyses/pandda_analyse_sites.csv')) + ' ' +
                      str(os.path.join(self.proasis_directory, 'LabXChem', self.proasis_name, 'reference'))))
        # copy reference pdb (from pandda directory to make sure same as in sites file)
        os.system(str('cp ' + str(os.path.join(self.panddas_directory, 'reference/reference.pdb')) + ' ' +
                      str(os.path.join(self.proasis_directory, 'LabXChem', self.proasis_name, 'reference'))))
        # open a temporary job file to write to for proasis scheduling
        temp_job = open(
            os.path.join(self.proasis_directory, 'Scripts/scheduled_jobs/temp_jobs', str(self.proasis_name + '.sh')),
            'w')
        # change file permissions of job
        perm_string = str(
            'chmod 770 ' + os.path.join(self.proasis_directory, 'Scripts/scheduled_jobs/temp_jobs',
                                        str(self.proasis_name + '.sh')))
        os.system(perm_string)
        # string to add leads in temp job file
        job_string = str('/dls/science/groups/proasis/Scripts/generate_leads.py -n ' + self.proasis_name
                         + ' -r '
                         + str(
            os.path.join(self.proasis_directory, 'LabXChem', self.proasis_name, 'reference/reference.pdb'))
                         + ' -p '
                         + str(os.path.join(self.proasis_directory, 'LabXChem', self.proasis_name,
                                            'reference/pandda_analyse_sites.csv'))
                         + ' -d '
                         + str(self.current_directory))

        temp_job.write(str(job_string))
        temp_job.close()
        # remove option from menu so lead can't be added multiple times
        self.proasis_lead.setVisible(False)

    def add_hits(self):
        # open the list of pernament jobs to append
        perm_job = open(os.path.join(self.proasis_directory, 'Scripts/scheduled_jobs/test.sh'), 'a')
        # string for job to add and update hits in proasis
        job_string = (str(os.path.join(self.proasis_directory, 'Scripts/populate_hits.py') + ' -d ' +
                          self.current_directory + ' > ' + os.path.join(self.proasis_directory,
                                                                       'Scripts/scheduled_jobs_logs',
                                                                       str(self.proasis_name + '_proasis.out'))))
        perm_job.write(job_string)
        perm_job.close()
        # remove option from menu so job is not added multiple times
        self.proasis_hits.setVisible(False)


if __name__ == "__main__":
	app=XChemExplorer(sys.argv[1:])


