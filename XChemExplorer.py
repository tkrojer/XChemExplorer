import os, sys, glob
from datetime import datetime

from PyQt4 import QtGui, QtCore, QtWebKit

import pickle
import base64
import math
import multiprocessing
import webbrowser

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))
from XChemUtils import parse
from XChemUtils import external_software
from XChemUtils import helpers
import XChemThread
import XChemDB
import XChemDialogs
import XChemPANDDA
import XChemToolTips
import XChemMain
import XChemPlots
import XChemLog

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import numpy as np

class XChemExplorer(QtGui.QApplication):
    def __init__(self,args):
        QtGui.QApplication.__init__(self,args)

        # general settings
        self.allowed_unitcell_difference_percent=5
        self.acceptable_low_resolution_limit_for_data=3
        self.filename_root='${samplename}'
        self.data_source_set=False

        #
        # directories
        #

        self.current_directory=os.getcwd()
        self.xce_logfile=os.path.join(self.current_directory,'xce.log')
        XChemLog.startLog(self.xce_logfile).create_logfile()
        self.update_log=XChemLog.updateLog(self.xce_logfile)
        self.update_log.insert('new session started')

        if 'labxchem' in self.current_directory:
            self.labxchem_directory='/'+os.path.join(*self.current_directory.split('/')[1:6])    # need splat operator: *
            self.beamline_directory=os.path.join(self.labxchem_directory,'processing','beamline')
            self.initial_model_directory=os.path.join(self.labxchem_directory,'processing','analysis','initial_model')
            self.reference_directory=os.path.join(self.labxchem_directory,'processing','reference')
            self.database_directory=os.path.join(self.labxchem_directory,'processing','database')
            self.panddas_directory=os.path.join(self.labxchem_directory,'processing','analysis','panddas')
            self.data_collection_summary_file=os.path.join(self.database_directory,str(os.getcwd().split('/')[5])+'_summary.pkl')
            self.data_source_file=''
            if os.path.isfile(os.path.join(self.labxchem_directory,'processing','database','soakDBDataFile.sqlite')):
                self.data_source_file='soakDBDataFile.sqlite'
                self.database_directory=os.path.join(self.labxchem_directory,'processing','database')
                self.data_source_set=True
                self.db=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
                self.db.create_missing_columns()
#                self.update_header_and_data_from_datasource()
#                self.header,self.data=self.db.load_samples_from_data_source()
#                XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file)).create_missing_columns()
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

        else:
            self.beamline_directory=self.current_directory
            self.initial_model_directory=self.current_directory
            self.reference_directory=self.current_directory
            self.database_directory=self.current_directory
            self.data_source_file=''
            self.ccp4_scratch_directory=os.getenv('CCP4_SCR')
            self.panddas_directory=self.current_directory
            self.data_collection_summary_file=''


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

        self.settings =     {'current_directory':       self.current_directory,
                             'beamline_directory':      self.beamline_directory,
                             'data_collection_summary': self.data_collection_summary_file,
                             'initial_model_directory': self.initial_model_directory,
                             'panddas_directory':       self.panddas_directory,
                             'reference_directory':     self.reference_directory,
                             'database_directory':      self.database_directory,
                             'data_source':             os.path.join(self.database_directory,self.data_source_file),
                             'ccp4_scratch':            self.ccp4_scratch_directory,
                             'unitcell_difference':     self.allowed_unitcell_difference_percent,
                             'too_low_resolution_data': self.acceptable_low_resolution_limit_for_data,
                             'filename_root':           self.filename_root,
                             'preferences':             self.preferences,
                             'xce_logfile':             self.xce_logfile        }


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
        self.main_data_collection_table_exists=False
        self.timer_to_check_for_new_data_collection = QtCore.QTimer()
#        self.timer_to_check_for_new_data_collection.timeout.connect(self.check_for_new_autoprocessing_or_rescore(False))

        self.target_list,self.visit_list=XChemMain.get_target_and_visit_list(self.beamline_directory)

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

        #
        # checking for external software packages
        #

        self.external_software=external_software().check()

        # start GUI

        self.start_GUI()
        self.exec_()





    def start_GUI(self):


        # GUI setup
        self.window=QtGui.QWidget()
        self.window.setGeometry(0,0, 800,600)
        self.window.setWindowTitle("XChemExplorer")
        self.center_main_window()

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
        quit.triggered.connect(QtGui.qApp.quit)
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

        datasource_menu.addAction(reload_samples_from_datasource)
        datasource_menu.addAction(save_samples_to_datasource)
        datasource_menu.addAction(import_csv_file_into_datasource)
        datasource_menu.addAction(export_csv_file_into_datasource)
        datasource_menu.addAction(update_datasource)
        datasource_menu.addAction(select_columns_to_show)
        datasource_menu.addAction(create_new_data_source)

        preferences_menu = menu_bar.addMenu("&Preferences")
        show_preferences=QtGui.QAction('Edit Preferences',self.window)
        show_preferences.triggered.connect(self.show_preferences)
        preferences_menu.addAction(show_preferences)

        help = menu_bar.addMenu("&Help")


        ######################################################################################
        #
        # Workflow @ Task Containers
        #

        self.workflow =     [   'Overview',         # 0
                                'Datasets',         # 1
                                'Maps',             # 2
                                'PANDDAs',          # 3
                                'Refinement',       # 4
                                'Settings'   ]      # 5

        self.workflow_dict = {  self.workflow[0]:       'Overview',
                                self.workflow[1]:       'Datasets',
                                self.workflow[2]:       'Maps',
                                self.workflow[3]:       'PANDDAs',
                                self.workflow[4]:       'Refinement',
                                self.workflow[5]:       'Settings'   }

        self.workflow_widget_dict = {}

        #
        # @ Update from datasource button ###################################################
        #

        update_from_datasource_button=QtGui.QPushButton("Update From\nDatasource")
        update_from_datasource_button.setToolTip(XChemToolTips.update_from_datasource_button_tip())
        update_from_datasource_button.clicked.connect(self.datasource_menu_reload_samples)

        #
        # @ Datasets ########################################################################
        #

        self.dataset_tasks = [  'Get New Results from Autoprocessing',
                                'Save Files from Autoprocessing to Project Folder',
                                'Rescore Datasets',
                                'Read PKL file'         ]

        frame_dataset_task=QtGui.QFrame()
        frame_dataset_task.setFrameShape(QtGui.QFrame.StyledPanel)
        vboxTask=QtGui.QVBoxLayout()
        label=QtGui.QLabel(self.workflow_dict['Datasets'])
        label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        vboxTask.addWidget(label)
        hboxAction=QtGui.QHBoxLayout()
        self.dataset_tasks_combobox = QtGui.QComboBox()
        for task in self.dataset_tasks:
            self.dataset_tasks_combobox.addItem(task)
        self.dataset_tasks_combobox.setToolTip(XChemToolTips.dataset_task_tip())
        hboxAction.addWidget(self.dataset_tasks_combobox)
        vboxButton=QtGui.QVBoxLayout()
        dataset_task_run_button=QtGui.QPushButton("Run")
        dataset_task_run_button.setToolTip(XChemToolTips.dataset_task_run_button_tip())
        dataset_task_run_button.clicked.connect(self.button_clicked)
        vboxButton.addWidget(dataset_task_run_button)
        dataset_task_status_button=QtGui.QPushButton("Status")
        dataset_task_status_button.setToolTip(XChemToolTips.dataset_task_status_button_tip())
        dataset_task_status_button.clicked.connect(self.button_clicked)
        vboxButton.addWidget(dataset_task_status_button)
        hboxAction.addLayout(vboxButton)
        vboxTask.addLayout(hboxAction)
        frame_dataset_task.setLayout(vboxTask)

        self.workflow_widget_dict['Datasets']=[self.dataset_tasks_combobox,dataset_task_run_button,dataset_task_status_button]


        #
        # @ MAP & CIF files #######################################################################
        #

        self.map_cif_file_tasks = [ 'Run DIMPLE on All Autoprocessing MTZ files',
                                    'Run DIMPLE on selected MTZ files',
                                    'Create CIF/PDB/PNG file of ALL soaked compound',
                                    'Create CIF/PDB/PNG file of NEW soaked compounds'    ]

        frame_map_cif_file_task=QtGui.QFrame()
        frame_map_cif_file_task.setFrameShape(QtGui.QFrame.StyledPanel)
        vboxTask=QtGui.QVBoxLayout()
        label=QtGui.QLabel(self.workflow_dict['Maps'])
        label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        vboxTask.addWidget(label)
        hboxAction=QtGui.QHBoxLayout()
        self.map_cif_file_tasks_combobox = QtGui.QComboBox()
        for task in self.map_cif_file_tasks:
            self.map_cif_file_tasks_combobox.addItem(task)
        self.map_cif_file_tasks_combobox.setToolTip(XChemToolTips.map_cif_file_task_tip())
        hboxAction.addWidget(self.map_cif_file_tasks_combobox)
        vboxButton=QtGui.QVBoxLayout()
        map_cif_file_task_run_button=QtGui.QPushButton("Run")
        map_cif_file_task_run_button.setToolTip(XChemToolTips.map_cif_file_task_run_button_tip())
        map_cif_file_task_run_button.clicked.connect(self.button_clicked)
        vboxButton.addWidget(map_cif_file_task_run_button)
        map_cif_file_task_status_button=QtGui.QPushButton("Status")
        map_cif_file_task_status_button.setToolTip(XChemToolTips.map_cif_file_task_status_button_tip())
        map_cif_file_task_status_button.clicked.connect(self.button_clicked)
        vboxButton.addWidget(map_cif_file_task_status_button)
        hboxAction.addLayout(vboxButton)
        vboxTask.addLayout(hboxAction)
        frame_map_cif_file_task.setLayout(vboxTask)

        self.workflow_widget_dict['Maps']=[self.map_cif_file_tasks_combobox,map_cif_file_task_run_button,map_cif_file_task_status_button]

        #####################################################################################

        #
        # @ PANDDAs #########################################################################
        #

        self.panddas_file_tasks = [ 'pandda.analyse',
                                    'pandda.inspect',
                                    'Export PANDDA models',
                                    'Show HTML summary'     ]

        frame_panddas_file_task=QtGui.QFrame()
        frame_panddas_file_task.setFrameShape(QtGui.QFrame.StyledPanel)
        vboxTask=QtGui.QVBoxLayout()
        label=QtGui.QLabel(self.workflow_dict['PANDDAs'])
        label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        vboxTask.addWidget(label)
        hboxAction=QtGui.QHBoxLayout()
        self.panddas_file_tasks_combobox = QtGui.QComboBox()
        for task in self.panddas_file_tasks:
            self.panddas_file_tasks_combobox.addItem(task)
        self.panddas_file_tasks_combobox.setToolTip(XChemToolTips.panddas_file_task_tip())
        hboxAction.addWidget(self.panddas_file_tasks_combobox)
        vboxButton=QtGui.QVBoxLayout()
        panddas_file_task_run_button=QtGui.QPushButton("Run")
        panddas_file_task_run_button.setToolTip(XChemToolTips.panddas_file_task_run_button_tip())
        panddas_file_task_run_button.clicked.connect(self.button_clicked)
        vboxButton.addWidget(panddas_file_task_run_button)
        panddas_file_task_status_button=QtGui.QPushButton("Status")
        panddas_file_task_status_button.setToolTip(XChemToolTips.panddas_file_task_status_button_tip())
        panddas_file_task_status_button.clicked.connect(self.button_clicked)
        vboxButton.addWidget(panddas_file_task_status_button)
        hboxAction.addLayout(vboxButton)
        vboxTask.addLayout(hboxAction)
        frame_panddas_file_task.setLayout(vboxTask)

        self.workflow_widget_dict['PANDDAs']=[self.panddas_file_tasks_combobox,panddas_file_task_run_button,panddas_file_task_status_button]

        #####################################################################################

        #
        # @ Refine ##########################################################################
        #

        self.refine_file_tasks = [ 'Open COOT'     ]

        frame_refine_file_task=QtGui.QFrame()
        frame_refine_file_task.setFrameShape(QtGui.QFrame.StyledPanel)
        vboxTask=QtGui.QVBoxLayout()
        label=QtGui.QLabel(self.workflow_dict['Refinement'])
        label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        vboxTask.addWidget(label)
        hboxAction=QtGui.QHBoxLayout()
        self.refine_file_tasks_combobox = QtGui.QComboBox()
        for task in self.refine_file_tasks:
            self.refine_file_tasks_combobox.addItem(task)
        self.refine_file_tasks_combobox.setToolTip(XChemToolTips.refine_file_task_tip())
        hboxAction.addWidget(self.refine_file_tasks_combobox)
        vboxButton=QtGui.QVBoxLayout()
        refine_file_task_run_button=QtGui.QPushButton("Run")
        refine_file_task_run_button.setToolTip(XChemToolTips.refine_file_task_run_button_tip())
        refine_file_task_run_button.clicked.connect(self.button_clicked)
        vboxButton.addWidget(refine_file_task_run_button)
        refine_file_task_status_button=QtGui.QPushButton("Status")
        refine_file_task_status_button.setToolTip(XChemToolTips.refine_file_task_status_button_tip())
        refine_file_task_status_button.clicked.connect(self.button_clicked)
        vboxButton.addWidget(refine_file_task_status_button)
        hboxAction.addLayout(vboxButton)
        vboxTask.addLayout(hboxAction)
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

#        self.overview_figure, self.overview_axes = plt.subplots(nrows=1, ncols=1)
        self.overview_figure, self.overview_axes = plt.subplots(nrows=2, ncols=2)
        self.overview_canvas = FigureCanvas(self.overview_figure)
        self.update_summary_plot()
        self.overview_tab_dict['Summary'][1].addWidget(self.overview_canvas)



        ######################################################################################
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
                         'Dewar'    ]
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
                                                        'Rmerge\nLow',
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
#        self.dls_data_collection_vbox.addLayout(data_collection_button_hbox)


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
#                    label=QtGui.QLabel('x')
#                    vbox_for_frame.addWidget(label)
#                    frame.setLayout(vbox_for_frame)
                self.dewar_configuration_layout.addWidget(frame, position, puck)

#        self.data_collection_summarys_vbox_for_details=QtGui.QVBoxLayout()
#        self.data_collection_details_currently_on_display=None
        self.dls_tab_dict['Dewar'][1].addLayout(self.dewar_configuration_layout)




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


        self.reference_file_list=self.get_reference_file_list(' ')
        self.reference_file_selection_combobox = QtGui.QComboBox()
        self.populate_reference_combobox(self.reference_file_selection_combobox)
        initial_model_checkbutton_hbox.addWidget(self.reference_file_selection_combobox)

        self.tab_dict[self.workflow_dict['Maps']][1].addLayout(initial_model_checkbutton_hbox)
        self.initial_model_vbox_for_table=QtGui.QVBoxLayout()
        self.inital_model_column_list=[   'Sample ID',
                        'Run\nDimple',
                        'Resolution\n[Mn<I/sig(I)> = 1.5]',
                        'Dimple\nRcryst',
                        'Dimple\nRfree',
                        'DataProcessing\nSpaceGroup',
                        'Reference\nSpaceGroup',
                        'Difference\nUC Volume (%)',
                        'Reference File',
                        'DataProcessing\nUnitCell'      ]

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
#                                    'Compound\nImage',
                                    'DataCollection\nOutcome',
                                    'DataProcessing\nSpaceGroup',
                                    'DataProcessing\nResolutionHigh',
                                    'Refinement\nOutcome',
                                    'Refinement\nRcryst',
                                    'Refinement\nRfree' ]
        self.summary_table=QtGui.QTableWidget()
        self.summary_table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.summary_table.setSortingEnabled(True)
        self.summary_table.setHorizontalHeaderLabels(self.summary_column_name)
        self.summary_vbox_for_table.addWidget(self.summary_table)
        self.tab_dict[self.workflow_dict['Refinement']][1].addLayout(self.summary_vbox_for_table)

        ######################################################################################


        #
        # @ PANDDAs Tab ######################################################################
        #

        self.panddas_results_vbox=QtGui.QVBoxLayout()
        self.tab_dict[self.workflow_dict['PANDDAs']][1].addLayout(self.panddas_results_vbox)

        pandda_tab_widget = QtGui.QTabWidget()
        pandda_tab_list = [ 'pandda.analyse',
                            'Dataset Summary',
                            'Results Summary',
                            'Inspect Summary'  ]

        self.pandda_tab_dict={}
        for page in pandda_tab_list:
            tab=QtGui.QWidget()
            vbox=QtGui.QVBoxLayout(tab)
            pandda_tab_widget.addTab(tab,page)
            self.pandda_tab_dict[page]=[tab,vbox]

        self.pandda_analyse_hbox=QtGui.QHBoxLayout()
        self.pandda_tab_dict['pandda.analyse'][1].addLayout(self.pandda_analyse_hbox)

        # left hand side: table with information about available datasets
        self.pandda_column_name = [ 'Sample ID',
                                    'Refinement\nSpace Group',
                                    'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                    'Dimple\nRcryst',
                                    'Dimple\nRfree' ]
        self.pandda_analyse_data_table=QtGui.QTableWidget()
        self.pandda_analyse_data_table.setSortingEnabled(True)
        self.pandda_analyse_data_table.resizeColumnsToContents()
        self.pandda_analyse_data_table.setColumnCount(len(self.pandda_column_name))
        self.pandda_analyse_data_table.setHorizontalHeaderLabels(self.pandda_column_name)
        self.pandda_analyse_hbox.addWidget(self.pandda_analyse_data_table)

        # right hand side: input parameters for PANDDAs run

        frame=QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        self.pandda_analyse_input_params_vbox=QtGui.QVBoxLayout()

        pandda_input_dir_hbox=QtGui.QHBoxLayout()
        self.pandda_analyse_input_params_vbox.addWidget(QtGui.QLabel('data directory'))
        self.pandda_input_data_dir_entry = QtGui.QLineEdit()
        self.pandda_input_data_dir_entry.setText(os.path.join(self.initial_model_directory,'*'))
        self.pandda_input_data_dir_entry.setFixedWidth(400)
        pandda_input_dir_hbox.addWidget(self.pandda_input_data_dir_entry)
        self.select_pandda_input_dir_button=QtGui.QPushButton("Select Input Template")
        self.select_pandda_input_dir_button.clicked.connect(self.select_pandda_input_template)
        pandda_input_dir_hbox.addWidget(self.select_pandda_input_dir_button)
        self.pandda_analyse_input_params_vbox.addLayout(pandda_input_dir_hbox)

        pandda_pdb_style_hbox=QtGui.QHBoxLayout()
        pandda_pdb_style_hbox.addWidget(QtGui.QLabel('pdb style'))
        self.pandda_pdb_style_entry=QtGui.QLineEdit()
        self.pandda_pdb_style_entry.setText('dimple.pdb')
        self.pandda_pdb_style_entry.setFixedWidth(200)
        pandda_pdb_style_hbox.addWidget(self.pandda_pdb_style_entry)
        self.pandda_analyse_input_params_vbox.addLayout(pandda_pdb_style_hbox)

        pandda_mtz_style_hbox=QtGui.QHBoxLayout()
        pandda_mtz_style_hbox.addWidget(QtGui.QLabel('mtz style'))
        self.pandda_mtz_style_entry=QtGui.QLineEdit()
        self.pandda_mtz_style_entry.setText('dimple.mtz')
        self.pandda_mtz_style_entry.setFixedWidth(200)
        pandda_mtz_style_hbox.addWidget(self.pandda_mtz_style_entry)
        self.pandda_analyse_input_params_vbox.addLayout(pandda_mtz_style_hbox)

        pandda_output_dir_hbox=QtGui.QHBoxLayout()
        self.pandda_analyse_input_params_vbox.addWidget(QtGui.QLabel('output directory'))
        self.pandda_output_data_dir_entry = QtGui.QLineEdit()
        self.pandda_output_data_dir_entry.setText(self.panddas_directory)
        self.pandda_output_data_dir_entry.setFixedWidth(400)
        pandda_output_dir_hbox.addWidget(self.pandda_output_data_dir_entry)
        self.select_pandda_output_dir_button=QtGui.QPushButton("Select PANNDAs Directory")
        self.select_pandda_output_dir_button.clicked.connect(self.settings_button_clicked)
        pandda_output_dir_hbox.addWidget(self.select_pandda_output_dir_button)
        self.pandda_analyse_input_params_vbox.addLayout(pandda_output_dir_hbox)

        # qstat or local machine
        self.pandda_analyse_input_params_vbox.addWidget(QtGui.QLabel('submit'))
        self.pandda_submission_mode_selection_combobox = QtGui.QComboBox()
        self.pandda_submission_mode_selection_combobox.addItem('local machine')
        if self.external_software['qsub']:
            self.pandda_submission_mode_selection_combobox.addItem('qsub')
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_submission_mode_selection_combobox)

        self.pandda_analyse_input_params_vbox.addWidget(QtGui.QLabel('number of processors'))
        self.pandda_nproc=multiprocessing.cpu_count()-1
        self.pandda_nproc_entry = QtGui.QLineEdit()
        self.pandda_nproc_entry.setText(str(self.pandda_nproc).replace(' ',''))
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_nproc_entry)

        self.pandda_analyse_input_params_vbox.addWidget(QtGui.QLabel('order events by:'))
        self.pandda_sort_event_combobox = QtGui.QComboBox()
        self.pandda_sort_event_combobox.addItem('z_peak')
        self.pandda_sort_event_combobox.addItem('cluster_size')
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_sort_event_combobox)

        # run pandda on specific crystalform only
#        self.pandda_analyse_input_params_vbox.addWidget(QtGui.QLabel('Use Specific Crystal Form Only'))
#        self.pandda_analyse_crystal_from_selection_combobox = QtGui.QComboBox()
#        self.pandda_analyse_crystal_from_selection_combobox.currentIndexChanged.connect(self.pandda_analyse_crystal_from_selection_combobox_changed)
##        self.update_pandda_crystal_from_combobox()
#        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_analyse_crystal_from_selection_combobox)

        self.pandda_analyse_input_params_vbox.addWidget(QtGui.QLabel('\n\n\nExpert Parameters (need rarely changing):\n'))

        # minimum number of datasets
        self.pandda_analyse_input_params_vbox.addWidget(QtGui.QLabel('min_build_datasets'))
        self.pandda_min_build_dataset_entry = QtGui.QLineEdit()
        self.pandda_min_build_dataset_entry.setText('40')
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_min_build_dataset_entry)

        self.pandda_analyse_input_params_vbox.addWidget(QtGui.QLabel('max_new_datasets'))
        self.pandda_max_new_datasets_entry = QtGui.QLineEdit()
        self.pandda_max_new_datasets_entry.setText('300')
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_max_new_datasets_entry)

        self.pandda_analyse_input_params_vbox.addWidget(QtGui.QLabel('grid_spacing (default=0.6; higher values speed up calculations, but maps might be less pretty)'))
        self.pandda_grid_spacing_entry = QtGui.QLineEdit()
        self.pandda_grid_spacing_entry.setText('0.6')
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_grid_spacing_entry)

        self.pandda_analyse_input_params_vbox.addStretch(1)

        frame.setLayout(self.pandda_analyse_input_params_vbox)

#        # green 'Run Pandda' button (which is red when pandda run in progress
#        self.run_panddas_button=QtGui.QPushButton("Run PANDDAs")
#        self.run_panddas_button.clicked.connect(self.button_clicked)
#        self.run_panddas_button.setFixedWidth(200)
#        self.run_panddas_button.setFixedHeight(100)
#        self.color_run_panddas_button()
#        self.pandda_analyse_input_params_vbox.addWidget(self.run_panddas_button)

#        self.pandda_analyse_hbox.addLayout(self.pandda_analyse_input_params_vbox)
        self.pandda_analyse_hbox.addWidget(frame)

        #######################################################
        # next three blocks display html documents created by pandda.analyse
        self.pandda_initial_html_file=os.path.join(self.panddas_directory,'results_summaries','pandda_initial.html')
        self.pandda_analyse_html_file=os.path.join(self.panddas_directory,'results_summaries','pandda_analyse.html')
#        self.pandda_inspect_html_file=os.path.join(self.panddas_directory,'results_summaries','pandda_inspect.html')
        self.pandda_inspect_html_file=os.path.join(self.panddas_directory,'results_summaries','pandda_inspect.html')

        self.pandda_initial_html = QtWebKit.QWebView()
        self.pandda_tab_dict['Dataset Summary'][1].addWidget(self.pandda_initial_html)
        self.pandda_initial_html.load(QtCore.QUrl(self.pandda_initial_html_file))
        self.pandda_initial_html.show()

        self.pandda_analyse_html = QtWebKit.QWebView()
        self.pandda_tab_dict['Results Summary'][1].addWidget(self.pandda_analyse_html)
        self.pandda_analyse_html.load(QtCore.QUrl(self.pandda_analyse_html_file))
        self.pandda_analyse_html.show()

        self.pandda_inspect_html = QtWebKit.QWebView()
        self.pandda_tab_dict['Inspect Summary'][1].addWidget(self.pandda_inspect_html)
        self.pandda_inspect_html.load(QtCore.QUrl(self.pandda_inspect_html_file))
        self.pandda_inspect_html.show()


#        self.pandda_analyse_html = QtWebKit.QWebView()
#        self.pandda_inspect_html = QtWebKit.QWebView()

        self.panddas_results_vbox.addWidget(pandda_tab_widget)


#        panddas_button_hbox=QtGui.QHBoxLayout()
#        show_panddas_results=QtGui.QPushButton("Show PANDDAs Results")
#        show_panddas_results.clicked.connect(self.button_clicked)
#        panddas_button_hbox.addWidget(show_panddas_results)
##        reload_panddas_results=QtGui.QPushButton("Reload PANDDAs Results")
##        reload_panddas_results.clicked.connect(self.button_clicked)
##        panddas_button_hbox.addWidget(reload_panddas_results)
#        launch_panddas_inspect=QtGui.QPushButton("Launch pandda.inspect")
#        launch_panddas_inspect.clicked.connect(self.button_clicked)
#        panddas_button_hbox.addWidget(launch_panddas_inspect)
#        export_panddas_inspect=QtGui.QPushButton("Export PANDDA Models")
#        export_panddas_inspect.clicked.connect(self.button_clicked)
#        panddas_button_hbox.addWidget(export_panddas_inspect)
#        self.tab_dict[self.workflow_dict['PANDDAs']][1].addLayout(panddas_button_hbox)

        ######################################################################################




        #
        # @ Settings Tab #####################################################################
        #

        ######################################################################################
        # Settings Tab
        self.data_collection_vbox_for_settings=QtGui.QVBoxLayout()
        self.tab_dict[self.workflow_dict['Settings']][1].addLayout(self.data_collection_vbox_for_settings)

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nProject Directory:'))
        settings_hbox_initial_model_directory=QtGui.QHBoxLayout()
        self.initial_model_directory_label=QtGui.QLabel(self.initial_model_directory)
        settings_hbox_initial_model_directory.addWidget(self.initial_model_directory_label)
        settings_buttoon_initial_model_directory=QtGui.QPushButton('Select Project Directory')
        settings_buttoon_initial_model_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_initial_model_directory.addWidget(settings_buttoon_initial_model_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_initial_model_directory)


        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nReference Structure Directory:'))
        settings_hbox_reference_directory=QtGui.QHBoxLayout()
        self.reference_directory_label=QtGui.QLabel(self.reference_directory)
        settings_hbox_reference_directory.addWidget(self.reference_directory_label)
        settings_buttoon_reference_directory=QtGui.QPushButton('Select Reference Structure Directory')
        settings_buttoon_reference_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_reference_directory.addWidget(settings_buttoon_reference_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_reference_directory)

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nData Source:'))
#        settings_hbox_database_directory=QtGui.QHBoxLayout()
#        self.database_directory_label=QtGui.QLabel(self.database_directory)
#        settings_hbox_database_directory.addWidget(self.database_directory_label)
#        settings_buttoon_database_directory=QtGui.QPushButton('Select Data Source Directory')
#        settings_buttoon_database_directory.clicked.connect(self.settings_button_clicked)
#        settings_hbox_database_directory.addWidget(settings_buttoon_database_directory)
#        self.data_collection_vbox_for_settings.addLayout(settings_hbox_database_directory)
        settings_hbox_data_source_file=QtGui.QHBoxLayout()
        if self.data_source_file != '':
            self.data_source_file_label=QtGui.QLabel(os.path.join(self.database_directory,self.data_source_file))
        else:
            self.data_source_file_label=QtGui.QLabel('')
        settings_hbox_data_source_file.addWidget(self.data_source_file_label)
        settings_buttoon_data_source_file=QtGui.QPushButton('Select Data Source File')
        settings_buttoon_data_source_file.clicked.connect(self.settings_button_clicked)
        settings_hbox_data_source_file.addWidget(settings_buttoon_data_source_file)
#        create_new_data_source_button=QtGui.QPushButton("Create New Data\nSource (SQLite)")
#        create_new_data_source_button.clicked.connect(self.button_clicked)
#        settings_hbox_data_source_file.addWidget(create_new_data_source_button)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_data_source_file)

        #################
        # Data Collection
        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nData Collection Directory'))

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

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nCCP4_SCR Directory:'))
        settings_hbox_ccp4_scratch_directory=QtGui.QHBoxLayout()
        self.ccp4_scratch_directory_label=QtGui.QLabel(self.ccp4_scratch_directory)
        settings_hbox_ccp4_scratch_directory.addWidget(self.ccp4_scratch_directory_label)
        settings_buttoon_ccp4_scratch_directory=QtGui.QPushButton('Select CCP4_SCR Directory')
        settings_buttoon_ccp4_scratch_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_ccp4_scratch_directory.addWidget(settings_buttoon_ccp4_scratch_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_ccp4_scratch_directory)

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nPANDDAs directory:'))
        settings_hbox_panddas_directory=QtGui.QHBoxLayout()
        self.panddas_directory_label=QtGui.QLabel(self.panddas_directory)
        settings_hbox_panddas_directory.addWidget(self.panddas_directory_label)
        settings_button_panddas_directory=QtGui.QPushButton('Select PANNDAs Directory')
        settings_button_panddas_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_panddas_directory.addWidget(settings_button_panddas_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_panddas_directory)

#        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nSites of Interest:'))
#        self.sites_of_interest_input = QtGui.QTextEdit()
#        self.sites_of_interest_input.setFixedWidth(400)
#        self.data_collection_vbox_for_settings.addWidget(self.sites_of_interest_input)


#        self.data_collection_vbox_for_settings.addStretch(1)
        ######################################################################################

        ######################################################################################
        # Preferences
#        self.vbox_for_preferences=QtGui.QVBoxLayout()
#        self.tab_dict[self.workflow_dict['Preferences']][1].addLayout(self.vbox_for_preferences)
#
#        self.vbox_for_preferences.addWidget(QtGui.QLabel('Select amount of processed data you wish to copy to initial_model directory:'))
#        self.preferences_data_to_copy_combobox = QtGui.QComboBox()
#        for item in self.preferences_data_to_copy:
#            self.preferences_data_to_copy_combobox.addItem(item[0])
#        self.preferences_data_to_copy_combobox.currentIndexChanged.connect(self.preferences_data_to_copy_combobox_changed)
#        self.vbox_for_preferences.addWidget(self.preferences_data_to_copy_combobox)
#
#        self.vbox_for_preferences.addWidget(QtGui.QLabel('Dataset Selection Mechanism:'))
#        self.preferences_selection_mechanism_combobox = QtGui.QComboBox()
#        for item in self.preferences_selection_mechanism:
#            self.preferences_selection_mechanism_combobox.addItem(item)
#        self.preferences_selection_mechanism_combobox.currentIndexChanged.connect(self.preferences_selection_mechanism_combobox_changed)
#        self.vbox_for_preferences.addWidget(self.preferences_selection_mechanism_combobox)
#
#        self.vbox_for_preferences.addStretch(1)
#
        ######################################################################################



        self.status_bar=QtGui.QStatusBar()
        self.progress_bar=QtGui.QProgressBar()
        self.progress_bar.setMaximum(100)
        hbox_status=QtGui.QHBoxLayout()
        hbox_status.addWidget(self.status_bar)
        hbox_status.addWidget(self.progress_bar)

        vbox_main = QtGui.QVBoxLayout()
        vbox_main.addWidget(menu_bar)
        vbox_main.addWidget(self.main_tab_widget)

        hboxTaskFrames=QtGui.QHBoxLayout()
        hboxTaskFrames.addWidget(update_from_datasource_button)
        hboxTaskFrames.addWidget(frame_dataset_task)
        hboxTaskFrames.addWidget(frame_map_cif_file_task)
        hboxTaskFrames.addWidget(frame_panddas_file_task)
        hboxTaskFrames.addWidget(frame_refine_file_task)
        vbox_main.addLayout(hboxTaskFrames)

        vbox_main.addLayout(hbox_status)

        self.window.setLayout(vbox_main)

        # this can be excrutiatingly slow...
#        if self.data_source_set:
#            self.datasource_menu_reload_samples()

        self.status_bar.showMessage('Ready')
#        self.timer = QtCore.QBasicTimer()
#        self.window.showMaximized()
        self.window.show()

        if self.data_source_file != '':
            write_enabled=self.check_write_permissions_of_data_source()
            if not write_enabled:
                self.data_source_set=False


    def select_sample_for_dimple(self):
        indexes = self.initial_model_table.selectionModel().selectedRows()
        for index in sorted(indexes):
            xtal=str(self.initial_model_table.item(index.row(), 0).text())
            self.update_log.insert('%s is marked for DIMPLE' %index.row())
            self.initial_model_dimple_dict[xtal][0].setChecked(True)


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

        hbox=QtGui.QHBoxLayout()
        hbox.addWidget(QtGui.QLabel('XCE logfile:'))
        self.xce_logfile_label=QtGui.QLabel(self.xce_logfile)
        hbox.addWidget(self.xce_logfile_label)
        button=QtGui.QPushButton("Change")
        button.clicked.connect(self.set_xce_logfile)
        hbox.addWidget(button)
        vbox.addLayout(hbox)

        preferencesLayout.addLayout(vbox,0,0)

        preferences.exec_();


    def set_xce_logfile(self):
        file_name = str(QtGui.QFileDialog.getSaveFileName(self.window,'Save file', self.current_directory))
        self.xce_logfile=str(file_name)
        self.xce_logfile_label.setText(str(self.xce_logfile))
        if self.xce_logfile=='' or self.xce_logfile[self.xce_logfile.rfind('/')+1:]=='':
           print '==> XCE: invalid file format'
        else:
            XChemLog.startLog(self.xce_logfile).create_logfile()
            self.update_log=XChemLog.updateLog(self.xce_logfile)


    def select_datasource_columns_to_display(self):
#        self.data_source_columns_to_display, ok = XChemDialogs.select_columns_to_show(
#            os.path.join(self.database_directory,self.data_source_file)).return_selected_columns()
#        self.populate_and_update_data_source_table()
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
#        columns_to_ignore=['Sample ID','ID']
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
        print '==> reading samples from data source: ',os.path.join(self.database_directory,self.data_source_file)
        self.update_status_bar('reading samples from data source: '+os.path.join(self.database_directory,self.data_source_file))
        self.update_header_and_data_from_datasource()
        self.update_all_tables()
#        self.populate_and_update_data_source_table()
#        self.create_initial_model_table()

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
        self.work_thread=XChemThread.update_datasource_from_file_system(self.initial_model_directory,os.path.join(self.database_directory,self.data_source_file))
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.connect(self.work_thread, QtCore.SIGNAL("create_initial_model_table"),self.create_initial_model_table)
        self.work_thread.start()





    def on_context_menu(self, point):
        # show context menu
        for key in self.dewar_configuration_dict:
            if self.dewar_configuration_dict[key]==self.sender():
                self.dewar_label_active=key
        self.popMenu.exec_(self.sender().mapToGlobal(point))

    def on_context_menu_initial_model(self, point):
        # show context menu
        self.popMenu_for_initial_model_table.exec_(self.sender().mapToGlobal(point))


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
        combobox.addItem('...')
        for reference_file in self.reference_file_list:
            combobox.addItem(reference_file[0])

    def populate_target_selection_combobox(self,combobox):
        combobox.clear()
        for target in self.target_list:
            combobox.addItem(target)


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
            self.pandda_initial_html_file=os.path.join(self.panddas_directory,'results_summaries','pandda_initial.html')
            self.pandda_analyse_html_file=os.path.join(self.panddas_directory,'results_summaries','pandda_analyse.html')
            self.pandda_inspect_html_file=os.path.join(self.panddas_directory,'results_summaries','pandda_inspect.html')


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
#                        self.update_header_and_data_from_datasource()
#                        self.populate_and_update_data_source_table()
#                        self.create_initial_model_table()
#                else:
#                    XChemDB.data_source(self.settings['data_source']).create_empty_data_source_file()
#                self.data_source_set=True

            self.ccp4_scratch_directory=pickled_settings['ccp4_scratch']
            self.settings['ccp4_scratch']=self.ccp4_scratch_directory

            self.allowed_unitcell_difference_percent=pickled_settings['unitcell_difference']
            self.acceptable_low_resolution_limit_for_data=pickled_settings['too_low_resolution_data']

            reference_directory_temp=pickled_settings['reference_directory']
            if reference_directory_temp != self.reference_directory:
                self.reference_directory=reference_directory_temp
                self.settings['reference_directory']=self.reference_directory
                self.update_reference_files(' ')

            self.initial_model_directory_label.setText(self.initial_model_directory)
            self.panddas_directory_label.setText(self.panddas_directory)
            self.pandda_output_data_dir_entry.setText(self.panddas_directory)
            self.reference_directory_label.setText(self.reference_directory)
#            self.database_directory_label.setText(self.database_directory)
            self.beamline_directory_label.setText(self.beamline_directory)
            self.ccp4_scratch_directory_label.setText(self.ccp4_scratch_directory)
#            self.adjust_allowed_unit_cell_difference.setText(str(self.allowed_unitcell_difference_percent))
#            self.adjust_acceptable_low_resolution_limit.setText(str(self.acceptable_low_resolution_limit_for_data))
            self.reference_file_list=self.get_reference_file_list(' ')

#            if self.data_source_set:
#                self.update_all_tables()


        except KeyError:
            self.update_status_bar('Sorry, this is not a XChemExplorer config file!')

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


    def target_selection_combobox_activated(self,text):
        self.target=str(text)


    def check_status_rerun_dimple_on_all_autoprocessing_files(self):
        print 'hallo'


    def rerun_dimple_on_all_autoprocessing_files(self):
        print '==> XCE: running DIMPLE on ALL auto-processing files'
        job_list=[]
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
            self.check_before_running_dimple(job_list)

    def run_dimple_on_selected_autoprocessing_file(self):
        job_list=[]
        for xtal in sorted(self.initial_model_dimple_dict):
            if self.initial_model_dimple_dict[xtal][0].isChecked():
                db_dict=self.xtal_db_dict[xtal]
                if os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'],db_dict['DataProcessingMTZfileName'])) or \
                    os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'])):

                    if os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'],db_dict['DataProcessingMTZfileName'])):
                        mtzin=os.path.join(db_dict['DataProcessingPathToMTZfile'],db_dict['DataProcessingMTZfileName'])
                    elif os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'])):
                        mtzin=os.path.join(db_dict['DataProcessingPathToMTZfile'])


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
            self.check_before_running_dimple(job_list)


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
        msgBox.setText("Do you really want to run %s Dimple jobs?\nNote: we will not run more than 100 at once on the cluster!" %len(job_list))
        msgBox.addButton(QtGui.QPushButton('Go'), QtGui.QMessageBox.YesRole)
        msgBox.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)
        reply = msgBox.exec_();

        if reply == 0:
            self.status_bar.showMessage('preparing %s DIMPLE jobs' %len(job_list))
            self.work_thread=XChemThread.run_dimple_on_all_autoprocessing_files(    job_list,
                                                                                    self.initial_model_directory,
                                                                                    self.external_software,
                                                                                    self.ccp4_scratch_directory,
                                                                                    self.database_directory,
                                                                                    self.data_source_file)
            self.explorer_active=1
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.work_thread.start()


    def center_main_window(self):
        screen = QtGui.QDesktopWidget().screenGeometry()
        size = self.window.geometry()
        self.window.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)

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

    def update_all_tables(self):
        print '==> checking for new reference files'
        self.update_status_bar('checking for new reference files')
        self.reference_file_list=self.get_reference_file_list(' ')
        print '==> updating Overview table'
        self.update_status_bar('updating Overview table')
        self.populate_and_update_data_source_table()
        print '==> updating Maps table'
        self.update_status_bar('updating Maps table')
        self.create_initial_model_table()
        print '==> updating PANDDA table'
        self.update_status_bar('updating PANDDA table')
        self.populate_pandda_analyse_input_table()
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
#                self.update_header_and_data_from_datasource()
#                self.populate_and_update_data_source_table()
#                self.create_initial_model_table()
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
            self.pandda_initial_html_file=os.path.join(self.panddas_directory,'results_summaries','pandda_initial.html')
            self.pandda_analyse_html_file=os.path.join(self.panddas_directory,'results_summaries','pandda_analyse.html')
            self.pandda_inspect_html_file=os.path.join(self.panddas_directory,'results_summaries','pandda_inspect.html')

#        if self.data_source_set:
#            self.update_all_tables()



    def change_allowed_unitcell_difference_percent(self,text):
        try:
            self.allowed_unitcell_difference_percent=int(text)
            self.settings['unitcell_difference']=self.allowed_unitcell_difference_percent
        except ValueError:
            if str(text).find('.') != -1:
                self.allowed_unitcell_difference_percent=int(str(text)[:str(text).find('.')])
                self.settings['unitcell_difference']=self.allowed_unitcell_difference_percent
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
#                    self.database_directory_label.setText(str(self.database_directory))
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
                    else:
                        self.need_to_switch_main_tab(task_index)

    def get_status_of_workflow_milestone(self,instruction):
        if instruction=='Run DIMPLE on All Autoprocessing MTZ files':
            self.check_status_rerun_dimple_on_all_autoprocessing_files()

        elif instruction=='Create CIF/PDB/PNG file of ALL soaked compound' or \
             instruction=='Create CIF/PDB/PNG file of NEW soaked compounds':
            self.check_status_create_cif_pdb_png_files()

    def prepare_and_run_task(self,instruction):

        if instruction=='Get New Results from Autoprocessing':
            self.check_for_new_autoprocessing_or_rescore(False)

        elif instruction=="Save Files from Autoprocessing to Project Folder" :
            self.save_files_to_initial_model_folder()

        elif instruction=='Rescore Datasets':
            self.check_for_new_autoprocessing_or_rescore(True)

        elif instruction=="Read PKL file":
            summary = pickle.load( open( self.data_collection_summary_file, "rb") )
            self.create_widgets_for_autoprocessing_results_only(summary)

        elif instruction=='Run DIMPLE on All Autoprocessing MTZ files':
            self.rerun_dimple_on_all_autoprocessing_files()

        elif instruction=='Run DIMPLE on selected MTZ files':
            self.run_dimple_on_selected_autoprocessing_file()

        elif instruction=='Create CIF/PDB/PNG file of ALL soaked compound':
            self.create_cif_pdb_png_files('ALL')

        elif instruction=='Create CIF/PDB/PNG file of NEW soaked compounds':
            self.create_cif_pdb_png_files('NEW')

        elif instruction=='pandda.analyse':
            self.run_pandda_analyse()

        elif instruction=='pandda.inspect':
            self.run_pandda_inspect()

        elif instruction=='Export PANDDA models':
            self.run_pandda_export()

        elif instruction=='Show HTML summary':
            self.show_pandda_html_summary()

        elif instruction=="Open COOT":
            if not self.coot_running:
                print 'starting coot'
                self.work_thread=XChemThread.start_COOT(self.settings)
                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                self.work_thread.start()




    def set_new_reference_if_applicable(self):
        reference_root=str(self.reference_file_selection_combobox.currentText())
        self.update_reference_files(reference_root)
        for xtal in self.initial_model_dimple_dict:
            db_dict=self.xtal_db_dict[xtal]
            reference_file=self.find_suitable_reference_file(db_dict)
            smallest_uc_difference=min(reference_file,key=lambda x: x[1])
            reference_file_selection_combobox=self.initial_model_dimple_dict[xtal][1]
            if float(smallest_uc_difference[1]) < self.allowed_unitcell_difference_percent:
                index = reference_file_selection_combobox.findText(str(reference_file[0][0]), QtCore.Qt.MatchFixedString)
                reference_file_selection_combobox.setCurrentIndex(index)
            else:
                reference_file_selection_combobox.setCurrentIndex(0)

    def check_status_create_png_of_soaked_compound(self):
        number_of_samples=0
        running=0
        timestamp_list=[]
        cif_file_generated=0
        for folder in glob.glob(os.path.join(self.initial_model_directory,'*','compound')):
            number_of_samples += 1
            if os.path.isfile(os.path.join(folder,'ACEDRG_IN_PROGRESS')):
                running += 1
                timestamp=datetime.fromtimestamp(os.path.getmtime(os.path.join(folder,'ACEDRG_IN_PROGRESS'))).strftime('%Y-%m-%d %H:%M:%S')
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
            self.work_thread=XChemThread.NEW_read_autoprocessing_results_from_disc(self.visit_list,
                                                                                self.target,
                                                                                self.reference_file_list,
                                                                                self.database_directory,
                                                                                self.data_collection_dict,
                                                                                self.preferences,
                                                                                self.data_collection_summary_file,
                                                                                self.initial_model_directory,
                                                                                rescore_only,
                                                                                self.acceptable_low_resolution_limit_for_data,
                                                                                os.path.join(self.database_directory,self.data_source_file))
            self.explorer_active=1
            self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.connect(self.work_thread, QtCore.SIGNAL("create_widgets_for_autoprocessing_results_only"),
                                                 self.create_widgets_for_autoprocessing_results_only)
            self.work_thread.start()

    def save_files_to_initial_model_folder(self):
        self.work_thread=XChemThread.LATEST_save_autoprocessing_results_to_disc(self.dataset_outcome_dict,
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

    def run_pandda_analyse(self):
        counter=self.check_data_for_pandda_analyse()
        if counter < 10:
            self.status_bar.showMessage('pandda.analyse: not enough datasets found')
            print '==> XCE: pandda.analyse: not enough datasets found'
            return
        pandda_params = {
                'data_dir':             str(self.pandda_input_data_dir_entry.text()),
                'out_dir':              str(self.pandda_output_data_dir_entry.text()),
                'submit_mode':          str(self.pandda_submission_mode_selection_combobox.currentText()),
                'nproc':                str(self.pandda_nproc_entry.text()),
                'min_build_datasets':   str(self.pandda_min_build_dataset_entry.text()),
                'pdb_style':            str(self.pandda_pdb_style_entry.text()),
                'mtz_style':            str(self.pandda_mtz_style_entry.text()),
                'sort_event':           str(self.pandda_sort_event_combobox.currentText()),
                'N_datasets':           counter,
                'max_new_datasets':     str(self.pandda_max_new_datasets_entry.text()),
                'grid_spacing':         str(self.pandda_grid_spacing_entry.text())
                        }
        self.work_thread=XChemPANDDA.run_pandda_analyse(pandda_params)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def check_data_for_pandda_analyse(self):
        counter=0
        for file in glob.glob(os.path.join(str(self.pandda_input_data_dir_entry.text()),str(self.pandda_pdb_style_entry.text()))):
            if os.path.isfile(file):
                counter+=1
        self.status_bar.showMessage('pandda.analyse: found %s useable datasets' %counter)
        print '==> XCE: pandda.analyse: found %s useable datasets' %counter
        return counter

    def run_pandda_inspect(self):
        self.settings['panddas_directory']=str(self.pandda_output_data_dir_entry.text())
        print '==> XCE: starting pandda.inspect'
        self.work_thread=XChemThread.start_pandda_inspect(self.settings)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def run_pandda_export(self):
        self.settings['panddas_directory']=str(self.pandda_output_data_dir_entry.text())
        print '==> XCE: exporting PANDDA models'
        self.work_thread=XChemPANDDA.run_pandda_export(self.panddas_directory,os.path.join(self.database_directory,self.data_source_file),self.initial_model_directory)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()


    def show_pandda_html_summary(self):
        self.pandda_initial_html.load(QtCore.QUrl(self.pandda_initial_html_file))
        self.pandda_initial_html.show()
        self.pandda_analyse_html.load(QtCore.QUrl(self.pandda_analyse_html_file))
        self.pandda_analyse_html.show()
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
            compound_list.append([str(item[0]),compoundID,str(item[2])])
        if compound_list != []:
            print '==> XCE: trying to create cif and pdb files for ',len(compound_list),' compounds using ACEDRG...'
            if self.external_software['qsub']:
                print '==> XCE: will try sending ',len(compound_list),' jobs to your computer cluster!'
            else:
                print '==> XCE: apparently no cluster available, so will run',len(compound_list),' sequential jobs on one core of your local machine.'
                print '==> XCE: this could take a while...'
            self.explorer_active=1
            self.work_thread=XChemThread.create_png_and_cif_of_compound(self.external_software,
                                                                        self.initial_model_directory,
                                                                        compound_list,
                                                                        self.database_directory,
                                                                        self.data_source_file,
                                                                        todo,
                                                                        self.ccp4_scratch_directory )
            self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.work_thread.start()


    def check_status_create_cif_pdb_png_files(self):
        self.update_status_bar('Please check terminal window for details!')
#        samples_in_db=self.db.execute_statement("select CrystalName from mainTable where CrystalName is not NULL;")
#        smiles_for_sample=self.db.execute_statement("select CrystalName,compoundSMILES from mainTable where compoundSMILES is not NULL or compoundSMILES is not '';")
#        samples_with_data=self.db.execute_statement("select CrystalName from mainTable where DataCollectionOutcome is 'success';")
#        cif_files=self.db.execute_statement("select CrystalName,RefinementCIF from mainTable where RefinementCIF is not Null or RefinementCIF is not '';")
        print '==> XCE: summary for compounds:'
#        print '    * nr samples in datasource:',len(samples_in_db)
#        print '    * nr SMILES for samples:   ',len(smiles_for_sample)
#        print '    * nr samples with data:    ',len(samples_with_data)
#        print '    * nr CIF files created:    ',len(cif_files)
        print 'here:'
        print XChemMain.get_jobs_running_on_cluster()
        print XChemMain.get_datasource_summary(os.path.join(self.database_directory,self.data_source_file))
#        out_bytes = subprocess.check_output(['qstat'])
#        out_text = out_bytes.decode('utf-8')
#        jobs_on_cluster = subprocess.check_output(['qstat'])
#        jobs_running=0
#        for n,line in enumerate(jobs_on_cluster):
#
#        jobs_running=n
#        print '==> XCE: job info'
#        print '    * nr jobs currently running on cluster:',jobs_running
#        print '    * nr ACEDRG jobs submitted'
#        print '    * nr ACEDRG jobs waiting'
#        print '    * nr ACEDRG jobs finished'
#        print '    * time ACEDRG queue started'
#        print '    * expected time to finish'

    def show_html_summary_and_diffraction_image(self):
        for key in self.albula_button_dict:
            if self.albula_button_dict[key][0]==self.sender():
                print '==> XCE: showing html summary in firefox'
                self.show_html_summary_in_firefox(key)
                # dials image viewer is unavailable at the moment (23/06/2016)
                # because of a clash with the PANDDA ccp4 installation
#                print '==> XCE: starting dials.image_viewer'
#                self.work_thread=XChemThread.start_dials_image_viewer(self.albula_button_dict[key][1])
#                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
#                self.work_thread.start()



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
##            font.setFamily(_fromUtf8("Verdana"))
##            font =  self.horizontalHeader().font()
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
#            layout=self.data_collection_image_dict[xtal][0]
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
#            data_collection_table.setFixedHeight(300)
#            data_collection_table.horizontalHeader().setStretchLastSection(False)
#            data_collection_table.verticalHeader().setStretchLastSection(False)
            data_collection_table.itemSelectionChanged.connect(self.update_selected_autoproc_data_collection_summary_table)
            data_collection_table.cellClicked.connect(self.user_update_selected_autoproc_data_collection_summary_table)
            data_collection_table.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Minimum)

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
                        # we change the selection only if the user did not touch it, assuming that he/she knows best
                        if not selection_changed_by_user:
                            data_collection_table.selectRow(index)

        self.populate_data_collection_summary_table()

        #-----------------------------------------------------------------------------------------------


    def find_suitable_reference_file(self,db_dict):
        reference_file=[]
        dummy=['...', '', '', '', 0, '0']
        reference_file.append([dummy,999])
#        self.status_bar.showMessage('checking: '+str(os.path.join(db_dict['DataProcessingPathToMTZfile'],db_dict['DataProcessingMTZfileName'])))
        suitable_reference=[]
        for reference in self.reference_file_list:
            # first we need one in the same pointgroup
            if reference[5]==db_dict['DataProcessingPointGroup']:
#                try:
#                    difference=math.fabs(1-(float(db_dict['DataProcessingUnitCellVolume'])/float(reference[4])))*100
#                    if difference < self.allowed_unitcell_difference_percent:
#                        suitable_reference.append([reference,difference])
#                except ValueError:
#                    continue
                try:
                    difference=math.fabs(1-(float(db_dict['DataProcessingUnitCellVolume'])/float(reference[4])))*100
                    reference_file.append([reference,difference])
                except ValueError:
                    continue
#        if suitable_reference != []:
#            reference_file=min(suitable_reference,key=lambda x: x[1])
        return reference_file


    def create_initial_model_table(self):
#        self.update_header_and_data_from_datasource()
        column_name=self.db.translate_xce_column_list_to_sqlite(self.inital_model_column_list)

        for xtal in sorted(self.xtal_db_dict):
            new_xtal=False
            db_dict=self.xtal_db_dict[xtal]
            if str(db_dict['DataCollectionOutcome']).lower().startswith('success'):
                reference_file=self.find_suitable_reference_file(db_dict)
                print 'here'
                print reference_file
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
                    elif header[0]=='Run\nDimple':
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
                        self.initial_model_table.setItem(current_row, column, cell_text)
            if new_xtal:
                self.initial_model_dimple_dict[xtal]=[run_dimple,reference_file_selection_combobox]




    def pandda_analyse_crystal_from_selection_combobox_changed(self,i):
        crystal_form = self.pandda_analyse_crystal_from_selection_combobox.currentText()
        self.populate_pandda_analyse_input_table(crystal_form)

    def preferences_data_to_copy_combobox_changed(self,i):
        text = str(self.preferences_data_to_copy_combobox.currentText())
        for item in self.preferences_data_to_copy:
            if item[0] == text:
                self.preferences['processed_data_to_copy']=item[1]
                break

    def preferences_selection_mechanism_combobox_changed(self,i):
        text = str(self.preferences_selection_mechanism_combobox.currentText())
        self.preferences['dataset_selection_mechanism']=text

    def get_reference_file_list(self,reference_root):
        # check available reference files
        reference_file_list=[]
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
        for i in reference_file_list:
            print i
        return reference_file_list




    def dataset_outcome_combobox_change_outcome(self,text):
        outcome=str(text)
        for key in self.dataset_outcome_combobox_dict:
            if self.dataset_outcome_combobox_dict[key]==self.sender():
                dataset=key
        self.dataset_outcome_dict[key]=outcome

#        for button in self.dataset_outcome_dict[dataset]:
#            if str(button.text())==outcome:
#                if outcome.startswith('success'):
#                    button.setStyleSheet("font-size:9px;background-color: rgb(0,255,0)")
#                else:
#                    button.setStyleSheet("font-size:9px;background-color: rgb(255,0,0)")
#            else:
#
#                button.setStyleSheet("font-size:9px;background-color: "+self.dataset_outcome[str(button.text())])

    def set_run_dimple_flag(self,state):
        if state == QtCore.Qt.Checked:
            for key in self.initial_model_dimple_dict:
                self.initial_model_dimple_dict[key][0].setChecked(True)
        else:
            for key in self.initial_model_dimple_dict:
                self.initial_model_dimple_dict[key][0].setChecked(False)

    def show_data_collection_details(self,state):
        # first remove currently displayed widget
        if self.data_collection_details_currently_on_display != None:
            self.data_collection_details_currently_on_display.hide()
#            self.data_collection_summarys_vbox_for_details.removeWidget(self.data_collection_details_currently_on_display)
#            self.data_collection_details_currently_on_display.deleteLater()
#            self.data_collection_details_currently_on_display.setParent(None)
#            sip.delete(self.data_collection_details_currently_on_display)
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
#                            for column in range(self.data_collection_summary_table.columnCount()):
#                                try:
#                                    self.data_collection_summary_table.item(item[1], column).setBackground(QtGui.QColor(255,255,150))
#                                except AttributeError:
#                                    pass
#                    print self.data_collection_column_three_dict[key][0].frameGeometry().height()
                    self.data_collection_details_currently_on_display=self.data_collection_column_three_dict[key][0]
                    self.data_collection_summarys_vbox_for_details.addWidget(self.data_collection_details_currently_on_display)
#                    self.data_collection_summarys_vbox_for_details.setSizeConstraint(QtGui.QLayout.SetMinimumSize)
                    self.data_collection_details_currently_on_display.show()
            else:
                # un-check all other ones
                self.data_collection_summary_dict[key][3].setChecked(False)
#        print self.data_collection_summary_table.columnCount()
#            if xtal_in_table and mtz_already_in_inital_model_directory:
#                self.main_data_collection_table.item(current_row, 0).setBackground(QtGui.QColor(100,100,150))
#                self.main_data_collection_table.item(current_row, 1).setBackground(QtGui.QColor(100,100,150))

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
        column_name=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file)).translate_xce_column_list_to_sqlite(self.data_collection_summary_column_name)
#        new_xtal=False
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
            if not logfile_found:
                db_dict={}
            if logfile_found and not too_low_resolution:
                outcome="success"
            elif logfile_found and too_low_resolution:
                outcome="Failed - low resolution"
            else:
                outcome="Failed - unknown"
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

            # update data source
#            if db_dict != {}:
#                self.update_data_source(xtal,db_dict)

            row += 1

        self.data_collection_summary_table.resizeRowsToContents()
        self.data_collection_summary_table.resizeColumnsToContents()

        self.status_bar.showMessage('updating Overview table')
#        self.update_header_and_data_from_datasource()
#        self.populate_and_update_data_source_table()

        self.status_bar.showMessage('idle')

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
#        self.dewar_sample_configuration_dict=self.get_dewar_configuration()

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
#            self.data_collection_summary_table.resizeRowsToContents()
#            self.data_collection_summary_table.resizeColumnsToContents()

    def user_update_selected_autoproc_data_collection_summary_table(self):
        for key in self.data_collection_column_three_dict:
            if self.data_collection_column_three_dict[key][0]==self.sender():
                # the user changed the selection, i.e. no automated selection will update it
                self.update_log.insert('user changed selection')
                self.data_collection_column_three_dict[key][1]=True
                # need to also update if not yet done
                user_already_changed_selection=False
                for entry in self.data_collection_dict[key]:
                    if entry[0]=='user_changed_selection':
                        user_already_changed_selection=True
                if not user_already_changed_selection:
                    self.data_collection_dict[key].append(['user_changed_selection'])
                XChemMain.change_links_to_selected_data_collection_outcome(key,self.data_collection_dict,
                                                                           self.data_collection_column_three_dict,
                                                                           self.dataset_outcome_dict,
                                                                           self.initial_model_directory,
                                                                           os.path.join(self.database_directory,self.data_source_file),
                                                                           self.xce_logfile)

    def update_selected_autoproc_data_collection_summary_table(self):
        for key in self.data_collection_column_three_dict:
            if self.data_collection_column_three_dict[key][0]==self.sender():
                sample=key
                break
        indexes=self.sender().selectionModel().selectedRows()
        for index in sorted(indexes):
            selected_processing_result=index.row()

        for entry in self.data_collection_dict[sample]:
            if entry[0]=='logfile':
                if entry[7]==selected_processing_result:
                    db_dict=entry[6]
                    # update datasource
#                    print db_dict
                    self.update_data_source(sample,db_dict)

        # update Overview table
#        self.update_header_and_data_from_datasource()
#        self.populate_and_update_data_source_table()

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
                if row[item]==None:
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
#        self.mounted_crystal_table.setColumnCount(0)
        self.mounted_crystal_table.setColumnCount(len(self.data_source_columns_to_display))
#        self.mounted_crystal_table.setRowCount(0)

        # first get a list of all the samples that are already in the table and which will be updated
        samples_in_table=[]
        current_row = self.mounted_crystal_table.rowCount()
        for row in range(current_row):
            sampleID=str(self.mounted_crystal_table.item(row,0).text())      # this must be the case
            samples_in_table.append(sampleID)

        columns_to_show=self.get_columns_to_show(self.data_source_columns_to_display)
        n_rows=self.get_rows_with_sample_id_not_null_from_datasource()
#        self.mounted_crystal_table.setRowCount(n_rows)
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
                if row[item]==None:
                    cell_text.setText('')
                else:
                    cell_text.setText(str(row[item]))
                if self.data_source_columns_to_display[y]=='Sample ID':     # assumption is that column 0 is always sampleID
                    cell_text.setFlags(QtCore.Qt.ItemIsEnabled)             # and this field cannot be changed
                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                self.mounted_crystal_table.setItem(x, y, cell_text)
        self.mounted_crystal_table.setHorizontalHeaderLabels(self.data_source_columns_to_display)


    def populate_summary_table(self,header,data):
        self.summary_table.setColumnCount(len(self.summary_column_name))
        self.summary_table.setRowCount(0)
        self.summary_table.setRowCount(len(data))

        columns_to_show=self.get_columns_to_show(self.summary_column_name,header)
        for x,row in enumerate(data):
            for y,item in enumerate(columns_to_show):
                cell_text=QtGui.QTableWidgetItem()
#                if item=='Image':
#                    cell_text.setText('')
                if row[item]==None:
                    cell_text.setText('')
                else:
                    cell_text.setText(str(row[item]))
                if self.summary_column_name[y]=='Sample ID':     # assumption is that column 0 is always sampleID
                    cell_text.setFlags(QtCore.Qt.ItemIsEnabled)             # and this field cannot be changed
                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                self.summary_table.setItem(x, y, cell_text)
        self.summary_table.setHorizontalHeaderLabels(self.summary_column_name)


    def populate_pandda_analyse_input_table(self):
#        self.update_header_and_data_from_datasource()
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
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.pandda_analyse_data_table.setItem(current_row, column, cell_text)
            if new_xtal:
                self.pandda_analyse_input_table_dict[xtal]=[]



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

#    def get_columns_to_show(self,column_list,header_of_current_datasource):
#        # maybe I coded some garbage before, but I need to find out which column name in the
#        # data source corresponds to the actually displayed column name in the table
#        # reason being that the unique column ID for DB may not be nice to look at
#        columns_to_show=[]
#        for column in column_list:
#            # first find out what the column name in the header is:
#            column_name=''
#            for name in self.all_columns_in_data_source:
#                if column==name[1]:
#                    column_name=name[0]
#            for n,all_column in enumerate(header_of_current_datasource):
#                if column_name==all_column:
#                    columns_to_show.append(n)
#                    break
#        return columns_to_show

#    def get_rows_with_sample_id_not_null(self,header,data):
#        sample_id_column=self.get_columns_to_show(['Sample ID'],header)
#        n_rows=0
#        for row in data:
#            if not str(row[sample_id_column[0]]).lower() != 'none' or not str(row[sample_id_column[0]]).replace(' ','') == '':
#                n_rows+=1
#        return n_rows

    def get_rows_with_sample_id_not_null_from_datasource(self):
        sample_id_column=self.get_columns_to_show(['Sample ID'])
        n_rows=0
        for row in self.data:
            if not str(row[sample_id_column[0]]).lower() != 'none' or not str(row[sample_id_column[0]]).replace(' ','') == '':
                n_rows+=1
        return n_rows

    def update_data_source(self,sample,db_dict):
        data_source=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
#        try:
#            data_source=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
#            data_source.update_insert_data_source(sample,db_dict)
#        except sqlite3.OperationalError,NameError:
#            pass



if __name__ == "__main__":
#    print '\n\n\n'
#    print '     ######################################################################'
#    print '     #                                                                    #'
#    print '     # XCHEMEXPLORER - multi dataset analysis                             #'
#    print '     #                                                                    #'
#    print '     # Version: 0.1                                                       #'
#    print '     #                                                                    #'
#    print '     # Date:                                                              #'
#    print '     #                                                                    #'
#    print '     # Author: Tobias Krojer, Structural Genomics Consortium, Oxford, UK  #'
#    print '     #         tobias.krojer@sgc.ox.ac.uk                                 #'
#    print '     #                                                                    #'
#    print '     ######################################################################'
#    print '\n\n\n'


    app=XChemExplorer(sys.argv[1:])

