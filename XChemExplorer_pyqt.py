import os, sys, glob
from datetime import datetime
from PyQt4 import QtGui, QtCore
#from PyQt4.QtCore import QThread, SIGNAL

import time
import pickle
import base64
import math

sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
from XChemUtils import process
from XChemUtils import parse
from XChemUtils import queue
from XChemUtils import mtztools


class XChemExplorer(QtGui.QApplication):
    def __init__(self,args):
        QtGui.QApplication.__init__(self,args)
        self.start_GUI()
        self.exec_()


    def start_GUI(self):

        # Settings @ Directories
        self.current_directory=os.getcwd()
        self.project_directory='/'+os.path.join(*self.current_directory.split('/')[1:6])    # need splat operator: *
        self.beamline_directory=os.path.join(self.project_directory,'processing','beamline')
        self.initial_model_directory=os.path.join(self.project_directory,'processing','analysis','initial_model')
        self.refine_model_directory=os.path.join(self.project_directory,'processing','analysis','refine_model')
        self.reference_directory=os.path.join(self.project_directory,'processing','reference')
        self.database_directory=os.path.join(self.project_directory,'database')
        self.data_source_file='XChemExplorer.csv'

        self.settings =     {'current_directory':       self.current_directory,
                             'project_directory':       self.project_directory,
                             'beamline_directory':      self.beamline_directory,
                             'initial_model_directory': self.initial_model_directory,
                             'refine_model_directory':  self.refine_model_directory,
                             'reference_directory':     self.reference_directory,
                             'data_source':             os.path.join(self.database_directory,self.data_source_file) }




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
        self.data_collection_statistics_dict={}
        self.initial_model_dimple_dict={}       # contains toggle button if dimple should be run
        self.reference_file_list=[]

        self.target_list=[]
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


        # GUI setup
        self.window=QtGui.QWidget()

        self.window.setGeometry(0,0, 1200,800)
        self.window.setWindowTitle("XChemExplorer")
        self.center_main_window()
        
        # Menu Widget
        load=QtGui.QAction("Load Config File", self.window)
        save=QtGui.QAction("Save Config File", self.window)
        quit=QtGui.QAction("Quit", self.window)
        quit.setShortcut('Ctrl+Q')
        quit.triggered.connect(QtGui.qApp.quit)

        menu_bar = QtGui.QMenuBar()
        file = menu_bar.addMenu("&File")
        settings = menu_bar.addMenu("&Settings")
        help = menu_bar.addMenu("&Help")
        
        file.addAction(load)
        file.addAction(save)
        file.addAction(quit)
        
        # Tab widget
        tab_widget = QtGui.QTabWidget()
        tab_list = [    'Mounted Crystals',
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

        # Mounted Crystals Tab
        self.mounted_crystals_vbox_for_table=QtGui.QVBoxLayout()
        self.tab_dict['Mounted Crystals'][1].addLayout(self.mounted_crystals_vbox_for_table)
        mounted_crystals_button_hbox=QtGui.QHBoxLayout()


        get_mounted_crystals_button=QtGui.QPushButton("Load Samples From Datasource")
        get_mounted_crystals_button.clicked.connect(self.button_clicked)
        mounted_crystals_button_hbox.addWidget(get_mounted_crystals_button)

        create_png_of_soaked_compound_button=QtGui.QPushButton("Create PNG file of Soaked Compound")
        create_png_of_soaked_compound_button.clicked.connect(self.button_clicked)
        mounted_crystals_button_hbox.addWidget(create_png_of_soaked_compound_button)


        self.tab_dict['Mounted Crystals'][1].addLayout(mounted_crystals_button_hbox)


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
        print self.target
        self.target='ATAD2A'


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
        self.timer = QtCore.QBasicTimer()
#        self.window.addWidget(vbox)
#        statusBar=QtGui.QStatusBar()
#        self.window.setStatusBar(statusBar)
        self.window.show()

    def target_selection_combobox_activated(self,text):
        print text
        self.target=text

    def buttonClicked(self):
        print 'hallo'
#        self.statusBar().showMessage('HALLLO')
        self.statusBar().showMessage('HALLLO')
#        sender = self.sender()
#        self.statusBar().showMessage(sender.text() + ' was pressed')    
    
    def center_main_window(self):
        screen = QtGui.QDesktopWidget().screenGeometry()
        size = self.window.geometry()
        self.window.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)

    def button_clicked(self):
        if self.target != '' and self.explorer_active==0:
### --- for offline testing -------------------------------------------
#            if self.sender().text()=='Get New Results from Autoprocessing':
#                dict_list=[]
#                self.create_widgets_for_autoprocessing_results(dict_list)
### -------------------------------------------------------------------
### --- this works but disabled so that stuff can be tested offline ---
            if self.sender().text()=='Get New Results from Autoprocessing':
#                reference_file_list=self.get_reference_file_list()
                self.work_thread=read_autoprocessing_results_from_disc(self.visit_list,self.target,self.reference_file_list)
                self.explorer_active=1
                self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                self.connect(self.work_thread, QtCore.SIGNAL("create_widgets_for_autoprocessing_results"),
                                                         self.create_widgets_for_autoprocessing_results)
                self.work_thread.start()
### -------------------------------------------------------------------
            if self.sender().text()=="Save Files from Autoprocessing in 'inital_model' Folder":
                self.work_thread=save_autoprocessing_results_to_disc(self.dataset_outcome_dict,
                                                                     self.data_collection_table_dict,
                                                                     self.data_collection_statistics_dict)
                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                self.work_thread.start()

            if self.sender().text()=="Load All Samples":
                print 'hallo'

            if self.sender().text()=="Check for inital Refinement":
#                reference_file_list=self.get_reference_file_list()
                print "checking for initial refinement"
                self.explorer_active=1
                self.work_thread=read_intial_refinement_results(self.initial_model_directory,self.reference_file_list)
                self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                self.connect(self.work_thread, QtCore.SIGNAL("create_initial_model_table"),self.create_initial_model_table)
                self.work_thread.start()

            if self.sender().text()=="Open COOT":
                print 'found'
                if not self.coot_running:
                    print 'starting coot'
                    self.work_thread=start_COOT(self.settings)
                    self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                    self.work_thread.start()



#                and self.sender().text()=='Get New Results from Autoprocessing':
#            self.explorer_active=1
#            threading.Thread(target=self.read_autoprocessing_results_from_disc, args=()).start()
#            self.status_bar.showMessage('COOL')
#            self.task='read_autoprocessing_results_from_disc()'
#            self.workThread=read_autoprocessing_results_from_disc(self.visit_list,self.target)
#            self.connect(self.workThread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
#            self.connect(self.workThread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
#            self.workThread.start()
#            self.get_autoprocessing_results=self.read_autoprocessing_results_from_disc()
#            self.read_autoprocessing_results_from_disc()

#        elif self.target != '' and data=='load_mounted_crystals' and self.explorer_active==0:
#            self.explorer_active=1
#            threading.Thread(target=self.load_mounted_crystals, args=()).start()
#        elif self.target != '' and data=='load_samples' and self.explorer_active==0 \
#             and self.load_initial_models_first_time==0:
#            self.explorer_active=1
#            threading.Thread(target=self.InitialRefinement, args=()).start()
#            self.load_initial_models_first_time=1
#        elif self.target != '' and data=='write_files' and self.explorer_active==0:
#            self.explorer_active=1
#            threading.Thread(target=self.WRITEFILES, args=()).start()
#        elif self.target != '' and data=='run_dimple' and self.explorer_active==0:
#            self.explorer_active=1
#            threading.Thread(target=self.RunDimple, args=()).start()
#        elif self.target != '' and data=='refresh_inital_refinement' and self.explorer_active==0:
#            self.initial_model_liststore.clear()
#            self.initial_model_scrolled.remove(self.initial_model_treeview)
#            self.explorer_active=1
#            threading.Thread(target=self.InitialRefinement, args=()).start()
#        elif data=='get_queue_jobs' and self.explorer_active==0:
#            self.explorer_active=1
#            threading.Thread(target=self.show_queue_jobs, args=()).start()
#        elif data=='refresh_queue_jobs' and self.explorer_active==0:
#            self.explorer_active=1
#            self.queue_jobs_liststore.clear()
#            self.queue_control_scrolled.remove(self.queue_jobs_treeview)
#            threading.Thread(target=self.show_queue_jobs, args=()).start()
#        elif data=='remove_queue_jobs' and self.explorer_active==0:
#            self.explorer_active=1
#            threading.Thread(target=self.remove_queue_jobs, args=()).start()
#        elif data=='cancel':
#            self.window.destroy()
#            quit()
#        elif data=='initial_refinement_set_reference' and self.explorer_active==0:
#            self.explorer_active=1
#            try:
#                if self.reference_file_root != '':
#                    self.initial_model_liststore.clear()
#                    self.initial_model_scrolled.remove(self.initial_model_treeview)
#                    self.explorer_active=1
#                    threading.Thread(target=self.InitialRefinement, args=()).start()
#            except AttributeError:
#                print 'no samples loaded'
#        else:
#            buff = '-> please select a target first'
#            self.status_bar.push(self.context_id, buff)


    def update_progress_bar(self,progress):
        self.progress_bar.setValue(progress)

    def update_status_bar(self,message):
        self.status_bar.showMessage(message)

    def thread_finished(self):
        self.explorer_active=0
        self.update_progress_bar(0)
        self.update_status_bar('idle')

#    def create_table_for_dimple(self):
#
#        dimple_table_column_name = [    'SampleID',
#                                        'Run\nDimple',
#                                        'Resolution',
#                                        'Rcryst',
#                                        'Rfree',
#                                        'Space Group\nautoprocessing',
#                                        'Space Group\nreference',
#                                        'Difference\nUnit Cell Volume (%)',
#                                        'Unit Cell\nautoprocessing',
#                                        'Unit Cell\nreference',
#                                        'Reference File'    ]
#
#
#        table=QtGui.QTableWidget()
#        table.setSortingEnabled(True)
#        table.setRowCount(len(self.initial_model_directory+'/*'))
#        table.setColumnCount(11)





    def create_widgets_for_autoprocessing_results(self,dict_list):
        data_collection_dict=dict_list[0]
        self.data_collection_statistics_dict=dict_list[1]

        # reset the two dictionaries which contain the buttons and tables for each data collection
        self.dataset_outcome_dict={}
        self.data_collection_table_dict={}

### --- used temporarily to be able to test stuff offline ---
        pickle.dump(data_collection_dict,open('data_collection_dict.p','wb'))
        pickle.dump(self.data_collection_statistics_dict,open('data_collection_statistics_dict.p','wb'))
#        data_collection_dict = pickle.load( open(os.getenv('XChemExplorer_DIR')+"/tmp/data_collection_dict.p", "rb" ) )
#        self.data_collection_statistics_dict= pickle.load( open(os.getenv('XChemExplorer_DIR')+"/tmp/data_collection_statistics_dict.p", "rb" ) )
### ---------------------------------------------------------


        diffraction_data_column_name = ['Program',                       (255,0,40),
                                        'Run',                           (100,230,40),
                                        'SpaceGroup',                    (100,230,40),
                                        'Unit Cell',                     (100,230,40),
                                        'Resolution\nOverall',           (100,230,40),
                                        'Resolution\nInner Shell',       (100,230,40),
                                        'Resolution\nOuter Shell',       (100,230,40),
                                        'Rmerge\nOverall',               (100,230,40),
                                        'Rmerge\nInner Shell',           (100,230,40),
                                        'Rmerge\nOuter Shell',           (100,230,40),
                                        'Mn(I/sig(I))\nOverall',         (100,230,40),
                                        'Mn(I/sig(I))\nInner Shell',     (100,230,40),
                                        'Mn(I/sig(I))\nOuter Shell',     (100,230,40),
                                        'Completeness\nOverall',         (100,230,40),
                                        'Completeness\nInner Shell',     (100,230,40),
                                        'Completeness\nOuter Shell',     (100,230,40),
                                        'Multiplicity\nOverall',         (100,230,40),
                                        'Multiplicity\nInner Shell',     (100,230,40),
                                        'Multiplicity\nOuter Shell',     (100,230,40)  ]

        diffraction_data_column_name = ['Program',
                                        'Run',
                                        'SpaceGroup',
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

        table.setRowCount(len(data_collection_dict))
        table.setColumnCount(3)

        for row,key in enumerate(sorted(self.data_collection_statistics_dict)):
            self.dataset_outcome_dict[key]=[]
            # this is the main table

            # column 1: sample ID
            table.setItem(row, 0, QtGui.QTableWidgetItem(key))

            # column 2: data collection date
            table.setItem(row, 1, QtGui.QTableWidgetItem(data_collection_dict[key][0][0][1]))

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
            for run_number,run in enumerate(data_collection_dict[key][0]):
                label = QtGui.QLabel(run[0])
                layout.addWidget(label,(run_number)*2,0)
                if len(data_collection_dict[key][3]) != 0:
                    for column_number,column in enumerate(sorted(data_collection_dict[key][3])):
                        if run[0] in column[0]:
                            pixmap = QtGui.QPixmap()
                            pixmap.loadFromData(base64.b64decode(column[1]))
                            label = QtGui.QLabel()
                            label.resize(320,200)
                            label.setPixmap(pixmap.scaled(label.size(), QtCore.Qt.KeepAspectRatio))
                            layout.addWidget(label, (run_number)*2+1, column_number)
                vbox_cell.addLayout(layout)

            hbox_for_button_and_table=QtGui.QHBoxLayout()
            # dataset outcome buttons1
            dataset_outcome_groupbox=QtGui.QGroupBox()
            dataset_outcome_vbox=QtGui.QVBoxLayout()
            found_successful_data_collection=1
            for outcome in sorted(self.dataset_outcome):
                button=QtGui.QPushButton(outcome)
                button.setAutoExclusive(True)
                button.setCheckable(True)
                button.setStyleSheet("background-color: "+self.dataset_outcome[outcome])
                button.clicked.connect(self.dataset_outcome_button_change_color)
                self.dataset_outcome_dict[key].append(button)
                if outcome=='success':
                    button.setChecked(True)
                    button.setStyleSheet("background-color: rgb(0,255,0)")
                    found_successful_data_collection=1
                dataset_outcome_vbox.addWidget(button)
            if not found_successful_data_collection:
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

            for n,sample in enumerate(self.data_collection_statistics_dict[key]):
                for column,item in enumerate(diffraction_data_column_name):
                    for line in self.data_collection_statistics_dict[key]:
#                        print len(line),line
                        if len(line)==3:
                            print line[0]
                            print line[1]
                            print line[2]
                            if line[0]==item:
                                cell_text=QtGui.QTableWidgetItem()
                                cell_text.setText(str(line[1]))
                                data_collection_table.setItem(n, column, cell_text)
#                                r=line[2][0]
#                                g=line[2][1]
#                                b=line[2][2]
#                                data_collection_table.item(n,column).setBackground(QtGui.QColor(r,g,b))
                                break
                        if len(line)==2:
                            if line[0]=='best file':
                                if line[1]==True:
                                    data_collection_table.selectRow(n)

#                for column,item in enumerate(line[2:20]):
#                    cell_text=QtGui.QTableWidgetItem()
#                    cell_text.setText(str(item))
#                    data_collection_table.setItem(n, column, cell_text)
#                    data_collection_table.item(n,column).setBackground(QtGui.QColor(100,100,150))
#                    r=diffraction_data_column_name[2*column+1][0]
#                    g=diffraction_data_column_name[2*column+1][1]
#                    b=diffraction_data_column_name[2*column+1][2]
#                    data_collection_table.item(n,column).setBackground(QtGui.QColor(r,g,b))
#                    print 'hallo'
#                    print diffraction_data_column_name[2*column+1]
#                    data_collection_table.item(n,column).setBackground(QtGui.QColor(diffraction_data_column_name[2*column+1]))
#                    if line[27]==True:
#                        data_collection_table.selectRow(n)
            # some_list[start:stop:step]
            data_collection_table.setHorizontalHeaderLabels(diffraction_data_column_name[0::2])
            data_collection_table.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)

            # this is necessary to render table properly
#            data_collection_table.verticalHeader().setStretchLastSection(False)
            data_collection_table.resizeRowsToContents()
            data_collection_table.resizeColumnsToContents()
            data_collection_table.horizontalHeader().setStretchLastSection(False)
            data_collection_table.verticalHeader().setStretchLastSection(True)

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
                    print line[1]
                    run_dimple.setChecked(line[1])
                    self.initial_model_dimple_dict[line[0]]=run_dimple
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
                if os.path.isfile(self.reference_directory+'/'+reference_root+'.mtz'):
                    mtz_reference=mtztools(self.reference_directory+'/'+reference_root+'.mtz').get_all_values_as_dict()
                    spg_reference=mtz_reference['spacegroup']
                    unitcell_reference=mtz_reference['unitcell']
                    lattice_reference=mtz_reference['bravais_lattice']
                    unitcell_volume_reference=mtz_reference['unitcell_volume']
                    reference_file_list.append([reference_root,
                                                spg_reference,
                                                unitcell_reference,
                                                lattice_reference,
                                                unitcell_volume_reference])
        return reference_file_list


    def dataset_outcome_button_change_color(self):
#        print self.sender().text()
        for key in self.dataset_outcome_dict:
            for button in self.dataset_outcome_dict[key]:
                if button==self.sender():
                    dataset=key
        for button in self.dataset_outcome_dict[dataset]:
            if button==self.sender():
                if str(self.sender().text()).startswith('success'):
                    button.setStyleSheet("background-color: rgb(0,255,0)")
                else:
                    button.setStyleSheet("background-color: rgb(255,0,0)")
#                button.setStyleSheet("border-style: inset")
            else:
                print self.dataset_outcome[str(button.text())]
                button.setStyleSheet("background-color: "+self.dataset_outcome[str(button.text())])
#        print 'hallo'
#        print self.sender()

    def set_run_dimple_flag(self,state):
        if state == QtCore.Qt.Checked:
            print 'checked'
        else:
            print 'not checked'


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
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'parsing initial_models folder ->'+sample)

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

                    if 'xia2' in self.data_collection_statistics_dict[key][index.row()][1]:
#                        print self.data_collection_statistics_dict[key][index.row()][1]
                        print os.path.join(*self.data_collection_statistics_dict[key][index.row()][1].split('/')[:13])
#                print('Row %d is selected' % index.row())
                    if 'fast_dp' in self.data_collection_statistics_dict[key][index.row()][1]:
#                        print self.data_collection_statistics_dict[key][index.row()][1]
                        print os.path.join(*self.data_collection_statistics_dict[key][index.row()][1].split('/')[:12])


            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'writing files from data processing to inital_model folder ->'+key)

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
        print csv_out

        self.emit(QtCore.SIGNAL("finished()"))


class read_autoprocessing_results_from_disc(QtCore.QThread):
    def __init__(self,visit_list,target,reference_file_list):
        QtCore.QThread.__init__(self)
        self.visit_list=visit_list
        self.target=target
        self.reference_file_list=reference_file_list
        self.data_collection_dict={}
        self.data_collection_statistics_dict={}

    def run(self):
        number_of_visits_to_search=len(self.visit_list)
        search_cycle=1
        for visit_directory in sorted(self.visit_list):
            progress_step=100/float(len(glob.glob(os.path.join(visit_directory,'processed',self.target,'*'))))
            progress=0
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
                    run=runs[runs.rfind('/')+1:]
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
                    self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'Step 2 of 3: parsing aimless logfiles ->'+sample)
                    aimless_results=parse().GetAimlessLog(logfile)

#        diffraction_data_column_name = ['Program',                       (255,0,40),
#                                        'Run',                           (100,230,40),
#                                        'SpaceGroup',                    (100,230,40),
#                                        'Unit Cell',                     (100,230,40),
#                                        'Resolution\nOverall',           (100,230,40),
#                                        'Resolution\nInner Shell',       (100,230,40),
#                                        'Resolution\nOuter Shell',       (100,230,40),
#                                        'Rmerge\nOverall',               (100,230,40),
#                                        'Rmerge\nInner Shell',           (100,230,40),
#                                        'Rmerge\nOuter Shell',           (100,230,40),
#                                        'Mn(I/sig(I))\nOverall',         (100,230,40),
#                                        'Mn(I/sig(I))\nInner Shell',     (100,230,40),
#                                        'Mn(I/sig(I))\nOuter Shell',     (100,230,40),
#                                        'Completeness\nOverall',         (100,230,40),
#                                        'Completeness\nInner Shell',     (100,230,40),
#                                        'Completeness\nOuter Shell',     (100,230,40),
#                                        'Multiplicity\nOverall',         (100,230,40),
#                                        'Multiplicity\nInner Shell',     (100,230,40),
#                                        'Multiplicity\nOuter Shell',     (100,230,40)  ]

                    self.data_collection_statistics_dict[sample].append([
                        index,                                                                                      # 0
                        logfile,                                                                                    # 1
                        ['Program',                     aimless_results['AutoProc'],                                                        (255,0,40)],
                        ['Run',                         aimless_results['Run'],                                                             (100,230,40)],
                        ['SpaceGroup',                  aimless_results['SpaceGroup'],                                                      (100,230,40)],
                        ['Unit Cell',                   aimless_results['UnitCell'],                                                        (100,230,40)],
                        ['Resolution\nOverall',         aimless_results['ResolutionLow']+'-'+aimless_results['ResolutionHigh'],             (100,230,40)],
                        ['Resolution\nInner Shell',     aimless_results['ResolutionLow']+'-'+aimless_results['ResolutionLowInnerShell'],    (100,230,40)],
                        ['Resolution\nOuter Shell',     aimless_results['ResolutionHighOuterShell']+'-'+aimless_results['ResolutionHigh'],  (100,230,40)],
                        ['Rmerge\nOverall',             aimless_results['RmergeOverall'],                                                   (100,230,40)],
                        ['Rmerge\nInner Shell',         aimless_results['RmergeLow'],                                                       (100,230,40)],
                        ['Rmerge\nOuter Shell',         aimless_results['RmergeHigh'],                                                      (100,230,40)],
                        ['Mn(I/sig(I))\nOverall',       aimless_results['IsigOverall'],                                                     (100,230,40)],
                        ['Mn(I/sig(I))\nInner Shell',   aimless_results['IsigLow'],                                                         (100,230,40)],
                        ['Mn(I/sig(I))\nOuter Shell',   aimless_results['IsigHigh'],                                                        (100,230,40)],
                        ['Completeness\nOverall',       aimless_results['CompletenessOverall'],                                             (100,230,40)],
                        ['Completeness\nInner Shell',   aimless_results['CompletenessLow'],                                                 (100,230,40)],
                        ['Completeness\nOuter Shell',   aimless_results['CompletenessHigh'],                                                (100,230,40)],
                        ['Multiplicity\nOverall',       aimless_results['MultiplicityOverall'],                                             (100,230,40)],
                        ['Multiplicity\nInner Shell',   aimless_results['MultiplicityLow'],                                                 (100,230,40)],
                        ['Multiplicity\nOuter Shell',   aimless_results['MultiplicityHigh'],                                                (100,230,40)],
                        aimless_results['Lattice'],                                                                 # 21
                        float(aimless_results['UniqueReflectionsOverall']),                                         # 22
                        float(aimless_results['CompletenessOverall']),                                              # 23
                        float(aimless_results['IsigOverall']),                                                      # 24
                        float(aimless_results['UnitCellVolume']),                                                   # 25
                        float(aimless_results['RmergeLow']),                                                        # 26
                        ['best file',False]                                                                                       # 27
                                        ])


#                    self.data_collection_statistics_dict[sample].append([
#                                index,                                                                                      # 0
#                                logfile,                                                                                    # 1
#                                aimless_results['AutoProc'],                                                                # 2
#                                aimless_results['Run'],                                                                     # 3
#                                aimless_results['SpaceGroup'],                                                              # 4
#                                aimless_results['UnitCell'],                                                                # 5
#                                aimless_results['ResolutionLow']+'-'+aimless_results['ResolutionHigh'],                     # 6
#                                aimless_results['ResolutionLow']+'-'+aimless_results['ResolutionLowInnerShell'],            # 7
#                                aimless_results['ResolutionHighOuterShell']+'-'+aimless_results['ResolutionHigh'],          # 8
#                                aimless_results['RmergeOverall'],                                                           # 9
#                                aimless_results['RmergeLow'],                                                               # 10
#                                aimless_results['RmergeHigh'],                                                              # 11
#                                aimless_results['IsigOverall'],                                                             # 12
#                                aimless_results['IsigLow'],                                                                 # 13
#                                aimless_results['IsigHigh'],                                                                # 14
#                                aimless_results['CompletenessOverall'],                                                     # 15
#                                aimless_results['CompletenessLow'],                                                         # 16
#                                aimless_results['CompletenessHigh'],                                                        # 17
#                                aimless_results['MultiplicityOverall'],                                                     # 18
#                                aimless_results['MultiplicityLow'],                                                         # 19
#                                aimless_results['MultiplicityHigh'],                                                        # 20
#                                aimless_results['Lattice'],                                                                 # 21
#                                float(aimless_results['UniqueReflectionsOverall']),                                         # 22
#                                float(aimless_results['CompletenessOverall']),                                              # 23
#                                float(aimless_results['IsigOverall']),                                                      # 24
#                                float(aimless_results['UnitCellVolume']),                                                   # 25
#                                float(aimless_results['RmergeLow']),                                                        # 26
#                                False                                                                                       # 27
#                                        ])
            else:
                self.data_collection_statistics_dict[sample]+='###'*20
 #               print sample
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
                        self.data_collection_statistics_dict[sample][index][27]=True
            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.emit(QtCore.SIGNAL('create_widgets_for_autoprocessing_results'), [self.data_collection_dict,
                                                                               self.data_collection_statistics_dict])


if __name__ == "__main__":
    app=XChemExplorer(sys.argv)

