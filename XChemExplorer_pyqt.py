import os, sys, glob

sys.path.append('/dls_sw/apps/albula/3.1/dectris/albula/3.1/python')
import dectris.albula
#try:
#    import dectris.albula
#except ImportError:
#    pass

#from datetime import datetime
from PyQt4 import QtGui, QtCore
#from PyQt4.QtCore import QThread, SIGNAL

#import time
import pickle
import base64
#import math
import subprocess

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))
#from XChemUtils import process
from XChemUtils import parse
#from XChemUtils import queue
#from XChemUtils import mtztools
from XChemUtils import external_software
import XChemThread
import XChemDB
import XChemDialogs



class XChemExplorer(QtGui.QApplication):
    def __init__(self,args):
        QtGui.QApplication.__init__(self,args)

        # checking for external software packages
        self.external_software=external_software().check()

        # Settings @ Directories
        self.current_directory=os.getcwd()
        if 'labxchem' in self.current_directory:
            self.project_directory='/'+os.path.join(*self.current_directory.split('/')[1:6])    # need splat operator: *
            self.beamline_directory=os.path.join(self.project_directory,'processing','beamline')
            self.initial_model_directory=os.path.join(self.project_directory,'processing','analysis','initial_model')
            self.refine_model_directory=os.path.join(self.project_directory,'processing','analysis','refine_model')
            self.reference_directory=os.path.join(self.project_directory,'processing','reference')
            self.database_directory=os.path.join(self.project_directory,'processing','database')
            self.data_source_file='XChemExplorer.csv'
            self.ccp4_scratch_directory=os.path.join(self.project_directory,'processing','tmp')

            if not os.path.isdir(self.beamline_directory):
                os.mkdir(self.beamline_directory)
            if not os.path.isdir(os.path.join(self.project_directory,'processing','analysis')):
                os.mkdir(os.path.join(self.project_directory,'processing','analysis'))
            if not os.path.isdir(self.initial_model_directory):
                os.mkdir(self.initial_model_directory)
            if not os.path.isdir(self.reference_directory):
                os.mkdir(self.reference_directory)
            if not os.path.isdir(self.database_directory):
                os.mkdir(self.database_directory)
            if not os.path.isdir(self.ccp4_scratch_directory):
                os.mkdir(self.ccp4_scratch_directory)


        else:
            self.project_directory=self.current_directory
            self.beamline_directory=self.current_directory
            self.initial_model_directory=self.current_directory
            self.refine_model_directory=self.current_directory
            self.reference_directory=self.current_directory
            self.database_directory=self.current_directory
            self.data_source_file=''
            self.ccp4_scratch_directory=os.getenv('CCP4_SCR')


        # general settings
        self.allowed_unitcell_difference_percent=12

        self.settings =     {'current_directory':       self.current_directory,
                             'project_directory':       self.project_directory,
                             'beamline_directory':      self.beamline_directory,
                             'initial_model_directory': self.initial_model_directory,
                             'refine_model_directory':  self.refine_model_directory,
                             'reference_directory':     self.reference_directory,
                             'database_directory':      self.database_directory,
                             'data_source':             os.path.join(self.database_directory,self.data_source_file),
                             'ccp4_scratch':            self.ccp4_scratch_directory,
                             'external_software':       self.external_software,
                             'unitcell_difference':     self.allowed_unitcell_difference_percent    }

#        self.FindHitsDir=self.project_directory+'/processing/analysis/find_hits'
#        self.DatabaseDir=self.project_directory+'/processing/database'
#        self.reference_directory=self.project_directory+'/processing/reference'
#        self.reference_file_root='reference'


        # Settings @ Lists
        self.data_collection_list=[]
        self.visit_list=[]
        self.target=''
        self.dataset_outcome_dict={}            # contains the dataset outcome buttons
        self.data_collection_table_dict={}      # contains the dataset table
        self.data_collection_dict={}
        self.data_collection_statistics_dict={}
        self.initial_model_dimple_dict={}       # contains toggle button if dimple should be run
        self.reference_file_list=[]
        self.all_columns_in_data_source=XChemDB.data_source(os.path.join(self.database_directory,
                                                                         self.data_source_file)).return_column_list()

        self.target_list=[]

        try:
            if args[0]=='-h' or args[0]=='--help':
                print 'help'
            if args[0].startswith('/'):
                self.beamline_directory=args[0]
                self.target_list.append('*')
        except IndexError:
            pass

        for dir in glob.glob(self.beamline_directory+'/*'):
            self.visit_list.append(os.path.realpath(dir))
            for target in glob.glob(os.path.realpath(dir)+'/processed/*'):
                if target[target.rfind('/')+1:] not in ['results','README-log','edna-latest.html']:
                    if target[target.rfind('/')+1:] not in self.target_list:
                        self.target_list.append(target[target.rfind('/')+1:])


        # Settings @ Switches
        self.explorer_active=0
        self.coot_running=0
        self.progress_bar_start=0
        self.progress_bar_step=0
        self.albula = None
        self.show_diffraction_image = None
        # check if qstat is available
        try:
            subprocess.call(['qstat'])
            self.queueing_system_available=True
        except OSError:
            self.queueing_system_available=False


        self.start_GUI()
        self.exec_()


    def start_GUI(self):


        # GUI setup
        self.window=QtGui.QWidget()
        self.window.setGeometry(0,0, 1800,1100)
        self.window.setWindowTitle("XChemExplorer")
        self.center_main_window()
        
        # Menu Widget
        load=QtGui.QAction("Open Config File", self.window)
        load.setShortcut('Ctrl+O')
        load.triggered.connect(self.open_config_file)
        save=QtGui.QAction("Save Config File", self.window)
        save.setShortcut('Ctrl+S')
        save.triggered.connect(self.save_config_file)
        quit=QtGui.QAction("Quit", self.window)
        quit.setShortcut('Ctrl+Q')
        quit.triggered.connect(QtGui.qApp.quit)

        menu_bar = QtGui.QMenuBar()
        file = menu_bar.addMenu("&File")
#        settings = menu_bar.addMenu("&Settings")
        help = menu_bar.addMenu("&Help")
        
        file.addAction(load)
        file.addAction(save)
        file.addAction(quit)
        
        # Tab widget
        tab_widget = QtGui.QTabWidget()
        tab_list = [    'Data Source',
                        'DLS @ Data Collection',
                        'DLS @ Summary',
                        'Initial Model',
                        'PANDDAs',
                        'Summary & Refine',
                        'Queue Control',
                        'Settings'  ]
        self.tab_dict={}
        for page in tab_list:
            tab=QtGui.QWidget()
            vbox=QtGui.QVBoxLayout(tab)
            tab_widget.addTab(tab,page)
            self.tab_dict[page]=[tab,vbox]

        # Data Source Tab
        self.data_source_columns_to_display=[   'Sample ID',
                                                'Compound ID',
                                                'Smiles',
                                                'Compound Name',
                                                'Tag'           ]
        self.mounted_crystal_table=QtGui.QTableWidget()
        self.mounted_crystal_table.setSortingEnabled(True)
#        self.mounted_crystal_table.setColumnWidth(0,250)
        self.mounted_crystal_table.resizeColumnsToContents()
        self.mounted_crystals_vbox_for_table=QtGui.QVBoxLayout()
        self.tab_dict['Data Source'][1].addLayout(self.mounted_crystals_vbox_for_table)
        self.mounted_crystals_vbox_for_table.addWidget(self.mounted_crystal_table)

        mounted_crystals_button_hbox=QtGui.QHBoxLayout()
        get_mounted_crystals_button=QtGui.QPushButton("Load Samples From Datasource")
        get_mounted_crystals_button.clicked.connect(self.button_clicked)
        mounted_crystals_button_hbox.addWidget(get_mounted_crystals_button)
        save_mounted_crystals_button=QtGui.QPushButton("Save Samples To Datasource")
        save_mounted_crystals_button.clicked.connect(self.button_clicked)
        mounted_crystals_button_hbox.addWidget(save_mounted_crystals_button)
        create_png_of_soaked_compound_button=QtGui.QPushButton("Create PDB/CIF/PNG files of Compound")
        create_png_of_soaked_compound_button.clicked.connect(self.button_clicked)
        mounted_crystals_button_hbox.addWidget(create_png_of_soaked_compound_button)

        create_new_data_source_button=QtGui.QPushButton("Create New CSV Data Source")
        create_new_data_source_button.clicked.connect(self.button_clicked)
        mounted_crystals_button_hbox.addWidget(create_new_data_source_button)

        select_data_source_columns_to_display_button=QtGui.QPushButton("Select Columns")
        select_data_source_columns_to_display_button.clicked.connect(self.button_clicked)
        mounted_crystals_button_hbox.addWidget(select_data_source_columns_to_display_button)

        self.tab_dict['Data Source'][1].addLayout(mounted_crystals_button_hbox)


        # DLS @ Data Collection Tab
        self.data_collection_vbox_for_table=QtGui.QVBoxLayout()
        self.tab_dict['DLS @ Data Collection'][1].addLayout(self.data_collection_vbox_for_table)
        data_collection_button_hbox=QtGui.QHBoxLayout()
        get_data_collection_button=QtGui.QPushButton("Get New Results from Autoprocessing")
        get_data_collection_button.clicked.connect(self.button_clicked)
        data_collection_button_hbox.addWidget(get_data_collection_button)
        write_files_button=QtGui.QPushButton("Save Files from Autoprocessing in 'inital_model' Folder")
        write_files_button.clicked.connect(self.button_clicked)
        data_collection_button_hbox.addWidget(write_files_button)
        target_selection_combobox = QtGui.QComboBox()
        for target in self.target_list:
            target_selection_combobox.addItem(target)
        target_selection_combobox.activated[str].connect(self.target_selection_combobox_activated)
        data_collection_button_hbox.addWidget(target_selection_combobox)
        self.tab_dict['DLS @ Data Collection'][1].addLayout(data_collection_button_hbox)
        self.target=str(target_selection_combobox.currentText())
#        print self.target
#        self.target='ATAD2A'

        # DLS @ Summary
        data_collection_summary_list=[]
        self.data_collection_summary_column_name=[      'Sample ID',
                                                        'Puck',
                                                        'Position',
                                                        'Date',
                                                        'img1',
                                                        'img2',
                                                        'img3',
                                                        'img4',
                                                        'Program',
                                                        'Space\nGroup',
                                                        'Dataset\nOutcome',
                                                        'Resolution\nOverall',
                                                        'Rmerge\nInner Shell',
                                                        'Mn(I/sig(I))\nOuter Shell',
                                                        'Completeness\nOverall',
                                                        'Show Diffraction\nImage'
                                                        ]


        data_collection_summary_list.append(['']*len(self.data_collection_summary_column_name))
        self.data_collection_summary_table=QtGui.QTableWidget()
        self.data_collection_summary_table.setRowCount(len(data_collection_summary_list))
#        self.data_collection_summary_table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.data_collection_summary_table.setColumnCount(len(self.data_collection_summary_column_name))
        self.data_collection_summary_table.setSortingEnabled(True)
        for row,line in enumerate(data_collection_summary_list):
            for column,item in enumerate(line):
                cell_text=QtGui.QTableWidgetItem()
                cell_text.setText(str(item))
                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                self.data_collection_summary_table.setItem(row, column, cell_text)
        self.data_collection_summary_table.setHorizontalHeaderLabels(self.data_collection_summary_column_name)
        self.data_collection_summarys_vbox_for_table=QtGui.QVBoxLayout()
        self.tab_dict['DLS @ Summary'][1].addLayout(self.data_collection_summarys_vbox_for_table)
        self.data_collection_summarys_vbox_for_table.addWidget(self.data_collection_summary_table)


        # Initial Model Tab
        initial_model_checkbutton_hbox=QtGui.QHBoxLayout()
        select_sample_for_dimple = QtGui.QCheckBox('(de-)select all samples for DIMPLE')
        select_sample_for_dimple.toggle()
#        select_sample_for_dimple.connect(self.set_run_dimple_flag)
        initial_model_checkbutton_hbox.addWidget(select_sample_for_dimple)
        self.tab_dict['Initial Model'][1].addLayout(initial_model_checkbutton_hbox)
        self.initial_model_vbox_for_table=QtGui.QVBoxLayout()
        self.tab_dict['Initial Model'][1].addLayout(self.initial_model_vbox_for_table)
        initial_model_button_hbox=QtGui.QHBoxLayout()
        get_initial_model_button=QtGui.QPushButton("Check for inital Refinement")
        get_initial_model_button.clicked.connect(self.button_clicked)
        initial_model_button_hbox.addWidget(get_initial_model_button)
        run_dimple_button=QtGui.QPushButton("Run Dimple")
        run_dimple_button.clicked.connect(self.button_clicked)
        initial_model_button_hbox.addWidget(run_dimple_button)
        refresh_inital_model_button=QtGui.QPushButton("Refresh")
        refresh_inital_model_button.clicked.connect(self.button_clicked)
        initial_model_button_hbox.addWidget(refresh_inital_model_button)
        self.reference_file_list=self.get_reference_file_list()
        reference_file_selection_combobox = QtGui.QComboBox()
        for reference_file in self.reference_file_list:
            reference_file_selection_combobox.addItem(reference_file[0])
        initial_model_button_hbox.addWidget(reference_file_selection_combobox)
        set_new_reference_button=QtGui.QPushButton("Set New Reference (if applicable)")
        set_new_reference_button.clicked.connect(self.button_clicked)
        initial_model_button_hbox.addWidget(set_new_reference_button)
        self.tab_dict['Initial Model'][1].addLayout(initial_model_button_hbox)

        # Summary & Refine Tab

        # prgress table
        # Sample ID
        # tag
        # compound
        # dataset outcome
        # space group
        # resolution
        # R/Rfree
        # PANDDAs/ Averaging
        # Refinement outcome
        # ligand CC




        self.summary_vbox_for_table=QtGui.QVBoxLayout()
        self.tab_dict['Summary & Refine'][1].addLayout(self.summary_vbox_for_table)
        summary_button_hbox=QtGui.QHBoxLayout()
        load_all_samples_button=QtGui.QPushButton("Load All Samples")
        load_all_samples_button.clicked.connect(self.button_clicked)
        summary_button_hbox.addWidget(load_all_samples_button)
        refresh_all_samples_button=QtGui.QPushButton("Refresh All Samples")
        refresh_all_samples_button.clicked.connect(self.button_clicked)
        summary_button_hbox.addWidget(refresh_all_samples_button)
        open_cootl_button=QtGui.QPushButton("Open COOT")
        open_cootl_button.clicked.connect(self.button_clicked)
        summary_button_hbox.addWidget(open_cootl_button)
        self.tab_dict['Summary & Refine'][1].addLayout(summary_button_hbox)

        # Settings Tab
        self.data_collection_vbox_for_settings=QtGui.QVBoxLayout()
        self.tab_dict['Settings'][1].addLayout(self.data_collection_vbox_for_settings)
        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('Project Directory:'))
        settings_hbox_project_directory=QtGui.QHBoxLayout()
        self.project_directory_label=QtGui.QLabel(self.project_directory)
        settings_hbox_project_directory.addWidget(self.project_directory_label)
        settings_buttoon_project_directory=QtGui.QPushButton('Select Project Directory')
        settings_buttoon_project_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_project_directory.addWidget(settings_buttoon_project_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_project_directory)
        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nInitial Model Directory:'))
        settings_hbox_initial_model_directory=QtGui.QHBoxLayout()
        self.initial_model_directory_label=QtGui.QLabel(self.initial_model_directory)
        settings_hbox_initial_model_directory.addWidget(self.initial_model_directory_label)
        settings_buttoon_initial_model_directory=QtGui.QPushButton('Select Initial Model Directory')
        settings_buttoon_initial_model_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_initial_model_directory.addWidget(settings_buttoon_initial_model_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_initial_model_directory)

        settings_hbox_adjust_allowed_unit_cell_difference=QtGui.QHBoxLayout()
        self.adjust_allowed_unit_cell_difference_label=QtGui.QLabel('Max. Allowed Unit Cell Difference (%):')
        settings_hbox_adjust_allowed_unit_cell_difference.addWidget(self.adjust_allowed_unit_cell_difference_label)
        self.adjust_allowed_unit_cell_difference = QtGui.QLineEdit()
        self.adjust_allowed_unit_cell_difference.setText(self.allowed_unitcell_difference_percent)
        settings_hbox_adjust_allowed_unit_cell_difference.addWidget(self.adjust_allowed_unit_cell_difference)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_adjust_allowed_unit_cell_difference)

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nRefine Model Directory:'))
        settings_hbox_refine_model_directory=QtGui.QHBoxLayout()
        self.refine_model_directory_label=QtGui.QLabel(self.refine_model_directory)
        settings_hbox_refine_model_directory.addWidget(self.refine_model_directory_label)
        settings_buttoon_refine_model_directory=QtGui.QPushButton('Select Refine Model Directory')
        settings_buttoon_refine_model_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_refine_model_directory.addWidget(settings_buttoon_refine_model_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_refine_model_directory)
        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nReference Structure Directory:'))
        settings_hbox_reference_directory=QtGui.QHBoxLayout()
        self.reference_directory_label=QtGui.QLabel(self.reference_directory)
        settings_hbox_reference_directory.addWidget(self.reference_directory_label)
        settings_buttoon_reference_directory=QtGui.QPushButton('Select Reference Structure Directory')
        settings_buttoon_reference_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_reference_directory.addWidget(settings_buttoon_reference_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_reference_directory)
        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nData Source:'))
        settings_hbox_database_directory=QtGui.QHBoxLayout()

        self.database_directory_label=QtGui.QLabel(self.database_directory)
        settings_hbox_database_directory.addWidget(self.database_directory_label)
        settings_buttoon_database_directory=QtGui.QPushButton('Select Data Source Directory')
        settings_buttoon_database_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_database_directory.addWidget(settings_buttoon_database_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_database_directory)

        settings_hbox_data_source_file=QtGui.QHBoxLayout()
        self.data_source_file_label=QtGui.QLabel(self.data_source_file)
        settings_hbox_data_source_file.addWidget(self.data_source_file_label)
        settings_buttoon_data_source_file=QtGui.QPushButton('Select Data Source File')
        settings_buttoon_data_source_file.clicked.connect(self.settings_button_clicked)
        settings_hbox_data_source_file.addWidget(settings_buttoon_data_source_file)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_data_source_file)
        
        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nBeamline Directory:'))
        settings_hbox_beamline_directory=QtGui.QHBoxLayout()
        self.beamline_directory_label=QtGui.QLabel(self.beamline_directory)
        settings_hbox_beamline_directory.addWidget(self.beamline_directory_label)
        settings_buttoon_beamline_directory=QtGui.QPushButton('Select Beamline Directory')
        settings_buttoon_beamline_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_beamline_directory.addWidget(settings_buttoon_beamline_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_beamline_directory)

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nCCP4_SCR Directory:'))
        settings_hbox_ccp4_scratch_directory=QtGui.QHBoxLayout()
        self.ccp4_scratch_directory_label=QtGui.QLabel(self.ccp4_scratch_directory)
        settings_hbox_ccp4_scratch_directory.addWidget(self.ccp4_scratch_directory_label)
        settings_buttoon_ccp4_scratch_directory=QtGui.QPushButton('Select CCP4_SCR Directory')
        settings_buttoon_ccp4_scratch_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_ccp4_scratch_directory.addWidget(settings_buttoon_ccp4_scratch_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_ccp4_scratch_directory)

        self.data_collection_vbox_for_settings.addStretch(1)
        # ----------------------------------------------------

        self.status_bar=QtGui.QStatusBar()
        self.progress_bar=QtGui.QProgressBar()
        self.progress_bar.setMaximum(100)
        hbox_status=QtGui.QHBoxLayout()
        hbox_status.addWidget(self.status_bar)
        hbox_status.addWidget(self.progress_bar)

        vbox_main = QtGui.QVBoxLayout()
        vbox_main.addWidget(menu_bar)
        vbox_main.addWidget(tab_widget)
        vbox_main.addLayout(hbox_status)

        self.window.setLayout(vbox_main)

        self.status_bar.showMessage('Ready')
#        self.timer = QtCore.QBasicTimer()
        self.window.show()

    def open_config_file(self):
        file_name = QtGui.QFileDialog.getOpenFileName(self.window,'Open file', self.current_directory)
        try:
            pickled_settings = pickle.load(open(file_name,"rb"))
            self.project_directory=pickled_settings['project_directory']
            self.beamline_directory=pickled_settings['beamline_directory']
            self.initial_model_directory=pickled_settings['initial_model_directory']
            self.refine_model_directory=pickled_settings['refine_model_directory']
            self.reference_directory=pickled_settings['reference_directory']
            self.database_directory=pickled_settings['database_directory']
            self.data_source_file=pickled_settings['data_source']
            self.ccp4_scratch_directory=pickled_settings['ccp4_scratch']
            self.project_directory_label.setText(self.project_directory)
            self.initial_model_directory_label.setText(self.initial_model_directory)
            self.refine_model_directory_label.setText(self.refine_model_directory)
            self.reference_directory_label.setText(self.reference_directory)
            self.database_directory_label.setText(self.database_directory)
            self.beamline_directory_label.setText(self.beamline_directory)
            self.data_source_file_label.setText(self.data_source_file)
            self.ccp4_scratch_directory_label.setText(self.ccp4_scratch_directory)
            self.reference_file_list=self.get_reference_file_list()
        except KeyError:
            self.update_status_bar('Sorry, this is not a XChemExplorer config file!')

    def save_config_file(self):
        file_name = QtGui.QFileDialog.getSaveFileName(self.window,'Save file', self.current_directory)
        pickle.dump(self.settings,open(file_name,'wb'))



    def target_selection_combobox_activated(self,text):
        #str(self.initial_model_dimple_dict[sample][1].currentText()
        print str(text)
        self.target=str(text)

    def center_main_window(self):
        screen = QtGui.QDesktopWidget().screenGeometry()
        size = self.window.geometry()
        self.window.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)

    def settings_button_clicked(self):
        if self.sender().text()=='Select Project Directory':
            self.project_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.project_directory_label.setText(self.project_directory)
        if self.sender().text()=='Select Initial Model Directory':
            self.initial_model_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.initial_model_directory_label.setText(self.initial_model_directory)
        if self.sender().text()=='Select Refine Model Directory':
            self.refine_model_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.refine_model_directory_label.setText(self.refine_model_directory)
        if self.sender().text()=='Select Reference Structure Directory':
            self.reference_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.reference_directory_label.setText(self.reference_directory)
        if self.sender().text()=='Select Data Source Directory':
            self.database_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.database_directory_label.setText(self.database_directory)
        if self.sender().text()=='Select Data Source File':
            filepath=str(QtGui.QFileDialog.getOpenFileName(self.window,'Select File', self.database_directory))
            self.data_source_file =   filepath.split('/')[-1]
            self.database_directory = filepath[:filepath.rfind('/')]
            self.data_source_file_label.setText(self.data_source_file)
            self.database_directory_label.setText(str(self.database_directory))
        if self.sender().text()=='Select Beamline Directory':
            self.beamline_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.beamline_directory_label.setText(self.beamline_directory)
        if self.sender().text()=='Select CCP4_SCR Directory':
            self.ccp4_scratch_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.ccp4_scratch_directory_label.setText(self.ccp4_scratch_directory)
        self.settings =     {'current_directory':       self.current_directory,
                             'project_directory':       self.project_directory,
                             'beamline_directory':      self.beamline_directory,
                             'initial_model_directory': self.initial_model_directory,
                             'refine_model_directory':  self.refine_model_directory,
                             'reference_directory':     self.reference_directory,
                             'database_directory':      self.database_directory,
                             'data_source':             self.data_source_file,
                             'ccp4_scratch':            self.ccp4_scratch_directory }


    def button_clicked(self):
#        if self.target != '' and self.explorer_active==0:
        if self.explorer_active==0:
            if self.sender().text()=='Get New Results from Autoprocessing':
                print 'target: ',self.target
                self.work_thread=XChemThread.read_autoprocessing_results_from_disc(self.visit_list,
                                                                                   self.target,
                                                                                   self.reference_file_list,
                                                                                   self.database_directory)
                self.explorer_active=1
                self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                self.connect(self.work_thread, QtCore.SIGNAL("create_widgets_for_autoprocessing_results"),
                                                         self.create_widgets_for_autoprocessing_results)
                self.work_thread.start()

            if self.sender().text()=="Save Files from Autoprocessing in 'inital_model' Folder":
                self.work_thread=XChemThread.save_autoprocessing_results_to_disc(self.dataset_outcome_dict,
                                                                     self.data_collection_table_dict,
                                                                     self.data_collection_statistics_dict,
                                                                     self.database_directory,self.data_source_file,
                                                                     self.initial_model_directory)
                self.explorer_active=1
                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                self.work_thread.start()

            if self.sender().text()=="Load All Samples":
                print 'hallo'

            if self.sender().text()=="Check for inital Refinement":
                # first check if there is already content in the table and if so
                # delete checkbox and combobox widgets
                if self.initial_model_dimple_dict != {}:
                    for key in self.initial_model_dimple_dict:
                        print key, self.initial_model_dimple_dict[key]

                self.explorer_active=1
                self.work_thread=XChemThread.read_intial_refinement_results(self.initial_model_directory,
                                                                            self.reference_file_list,
                                                                            os.path.join(self.database_directory,
                                                                                         self.data_source_file),
                                                                            self.allowed_unitcell_difference_percent)
                self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                self.connect(self.work_thread, QtCore.SIGNAL("create_initial_model_table"),self.create_initial_model_table)
                self.work_thread.start()

            if self.sender().text()=="Refresh":
                for key in self.initial_model_dimple_dict:
                    print key, self.initial_model_dimple_dict[key]

            if self.sender().text()=="Open COOT":
                print 'found'
                if not self.coot_running:
                    print 'starting coot'
                    self.work_thread=XChemThread.start_COOT(self.settings)
                    self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                    self.work_thread.start()

            if self.sender().text()=="Run Dimple":
                self.explorer_active=1
                self.work_thread=XChemThread.run_dimple_on_selected_samples(self.settings,
                                                                self.initial_model_dimple_dict,
                                                                self.external_software,
                                                                self.ccp4_scratch_directory )
                self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                self.work_thread.start()

            if self.sender().text()=="Load Samples From Datasource":
                content=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file)).load_samples_from_data_source()
                header=content[0]
                data=content[1]
                self.populate_data_source_table(header,data)


            if self.sender().text()=="Create New CSV Data Source":
                file_name = str(QtGui.QFileDialog.getSaveFileName(self.window,'Save file', self.database_directory))
                #make sure that the file always has .csv extension
                if file_name.rfind('.') != -1:
                   file_name=file_name[:file_name.rfind('.')]+'.csv'
                else:
                    file_name=file_name+'.csv'
                print file_name
                XChemDB.data_source(os.path.join(file_name)).create_empty_data_source_file()


            if self.sender().text()=="Save Samples To Datasource":
                allRows = self.mounted_crystal_table.rowCount()
                out=''
                for row in range(allRows):
                    for i in range(5):
                        if self.mounted_crystal_table.item(row,i).text()=='':
                            print 'hallo'
                            break
                        else:
                            out+=self.mounted_crystal_table.item(row,i).text()+','
                    out+='\n'
                print out

            if self.sender().text()=="Select Columns":
                print 'hallo'
                # QDialog is the kind of widget that will help here, but for now I park this
#                self.dialogTextBrowser = XChemDialogs.select_data_source_columns()
#                self.dialogTextBrowser.exec_()

            if self.sender().text()=="Create PDB/CIF/PNG files of Compound":
                columns_to_read=['Sample ID',
                                 'Compound ID',
                                 'Smiles']
                # it's a bit pointless now to search for the position of the columns,
                # but later I may want to give the user the option to change the column order
                headercount = self.mounted_crystal_table.columnCount()
                column_positions=[]
                for item in columns_to_read:
                    for x in range(headercount):
                        headertext = self.mounted_crystal_table.horizontalHeaderItem(x).text()
                        if headertext==item:
                            column_positions.append(x)
                            break
                allRows = self.mounted_crystal_table.rowCount()
                compound_list=[]
                for row in range(allRows):
                    if self.mounted_crystal_table.item(row,column_positions[2]).text() != '':
                        compound_list.append([str(self.mounted_crystal_table.item(row,column_positions[0]).text()),
                                              str(self.mounted_crystal_table.item(row,column_positions[1]).text()),
                                              str(self.mounted_crystal_table.item(row,column_positions[2]).text())] )
                if compound_list != []:
                    self.explorer_active=1
                    self.work_thread=XChemThread.create_png_and_cif_of_compound(self.external_software,
                                                                                self.initial_model_directory,
                                                                                compound_list)
                    self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
                    self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
                    self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                    self.work_thread.start()
            if str(self.sender().text()).startswith('Show'):
                print str(self.sender().text())
                visit=str(self.sender().text()).split()[1]
                diffraction_image=str(self.sender().text()).split()[3]
                xtal=diffraction_image[:diffraction_image.find('_')]
                print diffraction_image
                print xtal
                print visit
                print os.path.join(self.beamline_directory,visit,self.target,xtal,diffraction_image)

                try:
                    self.show_diffraction_image.loadFile('')
                except dectris.albula.viewer.DNoObject:
                    self.albula = dectris.albula.openMainFrame()
                    self.show_diffraction_image = self.albula.openSubFrame()
                    self.show_diffraction_image.loadFile(os.path.join(self.beamline_directory,visit,self.target,xtal,diffraction_image))

#                self.albula = dectris.albula.openMainFrame()
#                self.show_diffraction_image = self.albula.openSubFrame()
#                self.show_diffraction_image.loadFile(os.path.join(self.beamline_directory,visit,self.target,xtal,diffraction_image))




#                    print self.mounted_crystal_table.item(row,0).text()
#                    for i in range(5):
#                        print
#                        if self.mounted_crystal_table.item(row,i).text()=='':
#                            print 'hallo'
#                            break
#                        else:
#                            out+=self.mounted_crystal_table.item(row,i).text()+','
#                    out+='\n'
#                   print out



    def update_progress_bar(self,progress):
        self.progress_bar.setValue(progress)

    def update_status_bar(self,message):
        self.status_bar.showMessage(message)

    def thread_finished(self):
        self.explorer_active=0
        self.update_progress_bar(0)
        self.update_status_bar('idle')


    def create_widgets_for_autoprocessing_results(self,dict_list):
        self.data_collection_dict=dict_list[0]
        self.data_collection_statistics_dict=dict_list[1]

        # reset the two dictionaries which contain the buttons and tables for each data collection
        self.dataset_outcome_dict={}
        self.data_collection_table_dict={}


        diffraction_data_column_name = ['Program',
                                        'Run',
                                        'Space\nGroup',
                                        'Unit Cell',
                                        'Resolution\nOverall',
                                        'Resolution\nInner Shell',
                                        'Resolution\nOuter Shell',
                                        'Rmerge\nOverall',
                                        'Rmerge\nInner Shell',
                                        'Rmerge\nOuter Shell',
                                        'Mn(I/sig(I))\nOverall',
                                        'Mn(I/sig(I))\nInner Shell',
                                        'Mn(I/sig(I))\nOuter Shell',
                                        'Completeness\nOverall',
                                        'Completeness\nInner Shell',
                                        'Completeness\nOuter Shell',
                                        'Multiplicity\nOverall',
                                        'Multiplicity\nInner Shell',
                                        'Multiplicity\nOuter Shell',    ]



        self.dataset_outcome = {    "success":                      "rgb(200,200,200)",
                                    "Failed - centring failed":     "rgb(200,200,200)",
                                    "Failed - no diffraction":      "rgb(200,200,200)",
                                    "Failed - processing barfs":    "rgb(200,200,200)",
                                    "Failed - loop empty":          "rgb(200,200,200)",
                                    "Failed - low resolution":      "rgb(200,200,200)",
                                    "Failed - no X-rays":           "rgb(200,200,200)",
                                    "Failed - unknown":             "rgb(200,200,200)"  }


        table=QtGui.QTableWidget()
        table.setSortingEnabled(True)
#        table.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
#        table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
#        table.setSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.MinimumExpanding)
#        table.resizeColumnsToContents()

        table.setRowCount(len(self.data_collection_dict))
        table.setColumnCount(3)

        for row,key in enumerate(sorted(self.data_collection_statistics_dict)):
            self.dataset_outcome_dict[key]=[]
            # this is the main table

            # column 1: sample ID
            sample_ID=QtGui.QTableWidgetItem(key)
            sample_ID.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
            table.setItem(row, 0, sample_ID)

            # column 2: data collection date
            data_collection_date_time=QtGui.QTableWidgetItem(self.data_collection_dict[key][0][0][1])
            data_collection_date_time.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
            table.setItem(row, 1, data_collection_date_time)

            # column 3:
            # ---------------------------------------------------------|
            # |                                                        |
            # | crystal images                                         |
            # |                                                        |
            # |--------------------------------------------------------|
            # |         |                                              |
            # | dataset |                                              |
            # | outcome | table for data processing results            |
            # | buttons |                                              |
            # |         |                                              |
            # |--------------------------------------------------------|

            cell_widget=QtGui.QWidget()
            vbox_cell=QtGui.QVBoxLayout(cell_widget)        # main vbox

            # crystal images
            layout = QtGui.QGridLayout()
            for run_number,run in enumerate(self.data_collection_dict[key][0]):
                #label = QtGui.QLabel(run[0]+' ('+run[1]+' @ '+run[2]+')')
                #layout.addWidget(label,(run_number)*2,0)
                if len(self.data_collection_dict[key][3]) != 0:
                    #label = QtGui.QLabel(run[0]+' ('+run[1]+' @ '+run[2]+')')
                    label = QtGui.QLabel(str(run))
                    layout.addWidget(label,(run_number)*2,0)
                    label.setSizePolicy ( QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
                    column_number=0
                    for column in sorted(self.data_collection_dict[key][3]):
                        if run[0] in column[0]:
                            pixmap = QtGui.QPixmap()
                            pixmap.loadFromData(base64.b64decode(column[1]))
                            label = QtGui.QLabel()
                            label.resize(320,200)
                            label.setPixmap(pixmap.scaled(label.size(), QtCore.Qt.KeepAspectRatio))
                            layout.addWidget(label, (run_number)*2+1, column_number)
                            column_number+=1
                vbox_cell.addLayout(layout)

            hbox_for_button_and_table=QtGui.QHBoxLayout()

            # dataset outcome buttons
            dataset_outcome_groupbox=QtGui.QGroupBox()
            dataset_outcome_vbox=QtGui.QVBoxLayout()
            found_successful_data_collection=1
            for outcome in sorted(self.dataset_outcome):
                button=QtGui.QPushButton(outcome)
                button.setAutoExclusive(True)
                button.setCheckable(True)
                button.setStyleSheet("font-size:9px;background-color: "+self.dataset_outcome[outcome])
#                button.setFixedHeight(14)
                button.clicked.connect(self.dataset_outcome_button_change_color)
                self.dataset_outcome_dict[key].append(button)
                dataset_outcome_vbox.addWidget(button)

            if not str(self.data_collection_statistics_dict[key][0][0]).startswith('#'):
                for button in self.dataset_outcome_dict[key]:
                    if button.text()=='success':
                        button.setChecked(True)
                        button.setStyleSheet("background-color: rgb(0,255,0)")
            else:
                for button in self.dataset_outcome_dict[key]:
                    if button.text()=='Failed - unknown':
                        button.setChecked(True)
                        button.setStyleSheet("background-color: rgb(255,0,0)")

            dataset_outcome_groupbox.setLayout(dataset_outcome_vbox)
            hbox_for_button_and_table.addWidget(dataset_outcome_groupbox)

            # table for data processing results
            data_collection_table=QtGui.QTableWidget()
            data_collection_table.setRowCount(len(self.data_collection_statistics_dict[key]))
            data_collection_table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
            data_collection_table.setColumnCount(len(diffraction_data_column_name))

            font = QtGui.QFont()
#            font.setFamily(_fromUtf8("Verdana"))
#            font =  self.horizontalHeader().font()
            font.setPointSize(8)
#            font.setPointSize(5)
#            self.setFont(font)
            data_collection_table.setFont(font)

            data_processing_success=True
            for n,sample in enumerate(self.data_collection_statistics_dict[key]):
                # failed data processing
                if str(self.data_collection_statistics_dict[key][0][0]).startswith('#'):
                    for column,header in enumerate(diffraction_data_column_name):
                        cell_text=QtGui.QTableWidgetItem()
                        cell_text.setText('#')
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        data_collection_table.setItem(n, column, cell_text)
                    data_processing_success=False
#                if not data_processing_success:
#                    break               # otherwise a table with 60 lines appears
                # successful data processing
                if data_processing_success:
                    for column,header in enumerate(diffraction_data_column_name):
#                    for item in self.data_collection_statistics_dict[key]:
                        for item in sample:
                            if isinstance(item, list):
                                if len(item)==3:
                                    if item[0]==header:
                                        cell_text=QtGui.QTableWidgetItem()
                                        cell_text.setText(str(item[1]))
                                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                                        r=item[2][0]
                                        g=item[2][1]
                                        b=item[2][2]
                                        data_collection_table.setItem(n, column, cell_text)
                                        data_collection_table.item(n,column).setBackground(QtGui.QColor(r,g,b))
                                if len(item)==2:
                                    if item[0]=='best file':
                                        if item[1]==True:
                                            data_collection_table.selectRow(n)
            # some_list[start:stop:step]
            data_collection_table.setHorizontalHeaderLabels(diffraction_data_column_name)
            data_collection_table.horizontalHeader().setFont(font)
#            table.horizontalHeaderItem(1).setTextAlignment(Qt.AlignHCenter)
            data_collection_table.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)

            # this is necessary to render table properly
#            data_collection_table.verticalHeader().setStretchLastSection(False)
            data_collection_table.resizeRowsToContents()
            data_collection_table.resizeColumnsToContents()
            data_collection_table.horizontalHeader().setStretchLastSection(False)
            data_collection_table.verticalHeader().setStretchLastSection(True)
            data_collection_table.itemSelectionChanged.connect(self.update_selected_autoproc_data_collection_summary_table)
            hbox_for_button_and_table.addWidget(data_collection_table)
            vbox_cell.addLayout(hbox_for_button_and_table)
            self.data_collection_table_dict[key]=data_collection_table

            cell_widget.setLayout(vbox_cell)
            table.setCellWidget(row, 2, cell_widget)
#            table.stretchLastSection()
            table.setColumnWidth(2, 1000)

#            table.setItem(row, 2, vbox)
        table.setHorizontalHeaderLabels(['Sample','Date',''])
#        table.resizeColumnsToContents()
        table.resizeRowsToContents()
        table.setLineWidth(10)
        self.data_collection_vbox_for_table.addWidget(table)

        self.populate_data_collection_summary_table()

        #-----------------------------------------------------------------------------------------------

    def create_initial_model_table(self,initial_model_list):

        self.initial_model_dimple_dict={}
        initial_model_column_name = [   'SampleID',
                                        'Run\nDimple',
                                        'Resolution',
                                        'Rcryst',
                                        'Rfree',
                                        'Space Group\nautoprocessing',
                                        'Space Group\nreference',
                                        'Difference\nUnit Cell Volume (%)',
                                        'Reference File',
                                        'Unit Cell\nautoprocessing',
                                        'Unit Cell\nreference'  ]

        initial_model_table=QtGui.QTableWidget()
        initial_model_table.setRowCount(len(initial_model_list))
        initial_model_table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        initial_model_table.setColumnCount(len(initial_model_list[0])-1)
        initial_model_table.setSortingEnabled(True)
        for n,line in enumerate(initial_model_list):
            for column,item in enumerate(line[:-1]):
                if column==1:
                    run_dimple = QtGui.QCheckBox()
                    run_dimple.toggle()
                    initial_model_table.setCellWidget(n, column, run_dimple)
                    run_dimple.setChecked(line[1])
#                    self.initial_model_dimple_dict[line[0]]=run_dimple
                elif column==8:
                    # don't need to connect, because only the displayed text will be read out
                    reference_file_selection_combobox = QtGui.QComboBox()
                    for reference_file in self.reference_file_list:
                        reference_file_selection_combobox.addItem(reference_file[0])
                    initial_model_table.setCellWidget(n, column, reference_file_selection_combobox)
                    index = reference_file_selection_combobox.findText(str(line[8]), QtCore.Qt.MatchFixedString)
                    reference_file_selection_combobox.setCurrentIndex(index)
                else:
                    cell_text=QtGui.QTableWidgetItem()
                    cell_text.setText(str(item))
                    cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                    initial_model_table.setItem(n, column, cell_text)
            self.initial_model_dimple_dict[line[0]]=[run_dimple,reference_file_selection_combobox]
#                r=item
#                g=item
#                b=item
#                initial_model_table.item(n,column).setBackground(QtGui.QColor(r,g,b))
        initial_model_table.setHorizontalHeaderLabels(initial_model_column_name)
        initial_model_table.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        initial_model_table.resizeColumnsToContents()
        self.initial_model_vbox_for_table.addWidget(initial_model_table)

    def get_reference_file_list(self):
        # check available reference files
        reference_file_list=[]
        for files in glob.glob(self.reference_directory+'/*'):
            if files.endswith('.pdb'):
                reference_root=files[files.rfind('/')+1:files.rfind('.')]
#                if os.path.isfile(self.reference_directory+'/'+reference_root+'.mtz'):
#                    mtz_reference=mtztools(self.reference_directory+'/'+reference_root+'.mtz').get_all_values_as_dict()
#                    spg_reference=mtz_reference['spacegroup']
#                    unitcell_reference=mtz_reference['unitcell']
#                    lattice_reference=mtz_reference['bravais_lattice']
#                    unitcell_volume_reference=mtz_reference['unitcell_volume']
#                    reference_file_list.append([reference_root,
#                                                spg_reference,
#                                                unitcell_reference,
#                                                lattice_reference,
#                                                unitcell_volume_reference])

                if os.path.isfile(os.path.join(self.reference_directory,reference_root+'.pdb')):
                    pdb_reference=parse().PDBheader(os.path.join(self.reference_directory,reference_root+'.pdb'))
                    spg_reference=pdb_reference['SpaceGroup']
                    unitcell_reference=pdb_reference['UnitCell']
                    lattice_reference=pdb_reference['Lattice']
                    unitcell_volume_reference=pdb_reference['UnitCellVolume']
                    reference_file_list.append([reference_root,
                                                spg_reference,
                                                unitcell_reference,
                                                lattice_reference,
                                                unitcell_volume_reference])



        return reference_file_list


    def dataset_outcome_button_change_color(self):
#        print self.sender().text()
        outcome=''
        for key in self.dataset_outcome_dict:
            for button in self.dataset_outcome_dict[key]:
                if button==self.sender():
                    dataset=key
        for button in self.dataset_outcome_dict[dataset]:
            if button==self.sender():
                outcome=self.sender().text()
                if str(self.sender().text()).startswith('success'):
                    button.setStyleSheet("font-size:9px;background-color: rgb(0,255,0)")
                else:
                    button.setStyleSheet("font-size:9px;background-color: rgb(255,0,0)")
#                button.setStyleSheet("border-style: inset")
            else:
                print self.dataset_outcome[str(button.text())]
                button.setStyleSheet("font-size:9px;background-color: "+self.dataset_outcome[str(button.text())])
        self.update_outcome_data_collection_summary_table(dataset,outcome)

    def set_run_dimple_flag(self,state):
        if state == QtCore.Qt.Checked:
            print 'checked'
        else:
            print 'not checked'

    def populate_data_collection_summary_table(self):
        # 1. get length of table
        # 2. delete all entries
        self.data_collection_summary_table.setRowCount(0)

        self.data_collection_summary_table.setRowCount(len(self.data_collection_statistics_dict))
        for row,key in enumerate(sorted(self.data_collection_statistics_dict)):
            # find which dataset_outcome_button is checked
            outcome=''
            for button in self.dataset_outcome_dict[key]:
                if button.isChecked():
                    outcome=button.text()

            # find which autoprocessing run was thought to be the best
            selected_processing_result=0
            for n,sample in enumerate(self.data_collection_statistics_dict[key]):
                # check which row was auto-selected
                for item in sample:
                    if isinstance(item, list):
                        if len(item)==2:
                            if item[0]=='best file':
                                if item[1]==True:
                                    selected_processing_result=n

#            print self.data_collection_statistics_dict[key][selected_processing_result]
            # find latest run
            # take all images from latest run and put in table
            tmp=[]
            for runs in self.data_collection_dict[key][0]:
                tmp.append(runs[1])
            latest_run=''
            visit=''
            for n,run in enumerate(self.data_collection_dict[key][0]):
                if run[1]==max(tmp):
                    latest_run=run[0]
                    visit=run[2]
            images_to_show=[]
            if latest_run != '':
                for image in self.data_collection_dict[key][3]:
                    if latest_run in image[0] and image[0].endswith('t.png'):
                        images_to_show.append(image)
#                    for column in sorted(self.data_collection_dict[key][3]):
#                        if run[0] in column[0]:
#                            pixmap = QtGui.QPixmap()
#                            pixmap.loadFromData(base64.b64decode(column[1]))
#                            label = QtGui.QLabel()
#                            label.resize(320,200)
#                            label.setPixmap(pixmap.scaled(label.size(), QtCore.Qt.KeepAspectRatio))
#                            layout.addWidget(label, (run_number)*2+1, column_number)
#                            column_number+=1

            image_number=0
            for column,header in enumerate(self.data_collection_summary_column_name):
                cell_text=QtGui.QTableWidgetItem()
                if header=='Sample ID':
                    cell_text.setText(str(key))
                if header=='Dataset\nOutcome':
                    cell_text.setText(outcome)
                if header.startswith('img'):
                    if len(images_to_show) > image_number:
                        pixmap = QtGui.QPixmap()
                        pixmap.loadFromData(base64.b64decode(images_to_show[image_number][1]))
                        image = QtGui.QLabel()
                        image.resize(80,50)
                        image.setPixmap(pixmap.scaled(image.size(), QtCore.Qt.KeepAspectRatio))
                        self.data_collection_summary_table.setCellWidget(row, column, image)
                        image_number+=1
                        continue
                    else:
                        cell_text.setText('')
                if header=='Puck':
                    puck='n/a'
                    if len(self.data_collection_dict[key])==5:
                        puck=self.data_collection_dict[key][4][0]
                    cell_text.setText(puck)
                if header=='Position':
                    position='n/a'
                    if len(self.data_collection_dict[key])==5:
                        position=self.data_collection_dict[key][4][1]
                    cell_text.setText(position)
                if header.startswith('Show'):
                    start_albula_button=QtGui.QPushButton('Show: '+visit+" @ "+latest_run+'0001.cbf')
                    start_albula_button.clicked.connect(self.button_clicked)
                    self.data_collection_summary_table.setCellWidget(row,column,start_albula_button)
                for item in self.data_collection_statistics_dict[key][selected_processing_result]:
                    if isinstance(item, list):
                        if len(item)==3:
                            if item[0]==header:
                                cell_text.setText(str(item[1]))

                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                self.data_collection_summary_table.setItem(row, column, cell_text)

        self.data_collection_summary_table.resizeRowsToContents()
        self.data_collection_summary_table.resizeColumnsToContents()

    def update_outcome_data_collection_summary_table(self,sample,outcome):
#        print 'hallo update'
#	    allRows=self.data_collection_summary_table.rowCount()
        rows_in_table=self.data_collection_summary_table.rowCount()
        for row in range(rows_in_table):
            if self.data_collection_summary_table.item(row,0).text()==sample:
#                print self.data_collection_summary_table.item(row,0).text()
                cell_text=QtGui.QTableWidgetItem()
                cell_text.setText(outcome)
                self.data_collection_summary_table.setItem(row, 3, cell_text)
#            self.data_collection_summary_table.resizeRowsToContents()
#            self.data_collection_summary_table.resizeColumnsToContents()

    def update_selected_autoproc_data_collection_summary_table(self):
        for key in self.data_collection_table_dict:
            if self.data_collection_table_dict[key]==self.sender():
                sample=key
                break

        indexes=self.sender().selectionModel().selectedRows()
        for index in sorted(indexes):
            selected_processing_result=index.row()


        rows_in_table=self.data_collection_summary_table.rowCount()
        for row in range(rows_in_table):
            if self.data_collection_summary_table.item(row,0).text()==sample:
#                print self.data_collection_summary_table.item(row,0).text()
                for column,header in enumerate(self.data_collection_summary_column_name):
                    cell_text=QtGui.QTableWidgetItem()
                    if header=='Sample ID':
                        continue
                    if header=='Dataset\nOutcome':
                        continue
                    for item in self.data_collection_statistics_dict[sample][selected_processing_result]:
                        if isinstance(item, list):
                            if len(item)==3:
                                if item[0]==header:
                                    cell_text.setText(str(item[1]))
                    cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                    self.data_collection_summary_table.setItem(row, column, cell_text)

    def populate_data_source_table(self,header,data):
        self.mounted_crystal_table.setColumnCount(len(self.data_source_columns_to_display))
        self.mounted_crystal_table.setRowCount(0)
        self.mounted_crystal_table.setRowCount(len(data))

        # maybe I coded some garbage before, but I need to find out which column name in the
        # data source corresponds to the actually displayed column name in the table
        # reason being that the unique column ID for DB may not be nice to look at
        columns_to_show=[]
        for column in self.data_source_columns_to_display:
#            print column
            for n,all_column in enumerate(self.all_columns_in_data_source):
                if column==all_column[1]:
                    columns_to_show.append(n)
                    break

        for x,row in enumerate(data):
#            print row
            y=0
            for item in columns_to_show:
                cell_text=QtGui.QTableWidgetItem()
                if row[item]==None:
                    cell_text.setText('')
                else:
                    cell_text.setText(str(row[item]))
                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                self.mounted_crystal_table.setItem(x, y, cell_text)
                y+=1
        self.mounted_crystal_table.setHorizontalHeaderLabels(self.data_source_columns_to_display)



if __name__ == "__main__":
    app=XChemExplorer(sys.argv[1:])

