import os, sys, glob
from datetime import datetime
from PyQt4 import QtGui, QtCore
#from PyQt4.QtCore import QThread, SIGNAL

import time
import pickle
import base64

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

#        self.FindHitsDir=self.project_directory+'/processing/analysis/find_hits'
#        self.DatabaseDir=self.project_directory+'/processing/database'
#        self.reference_directory=self.project_directory+'/processing/reference'
#        self.reference_file_root='reference'


        # Settings @ Lists
        self.data_collection_list=[]
        self.visit_list=[]
        self.target=''
        self.dataset_outcome_dict={}        # contains the dataset outcome buttons
        self.target_list=[]
        for dir in glob.glob(self.beamline_directory+'/*'):
            self.visit_list.append(os.path.realpath(dir))
            for target in glob.glob(os.path.realpath(dir)+'/processed/*'):
                if target[target.rfind('/')+1:] not in ['results','README-log','edna-latest.html']:
                    if target[target.rfind('/')+1:] not in self.target_list:
                        self.target_list.append(target[target.rfind('/')+1:])


        # Settings @ Switches
        self.explorer_active=0
        self.progress_bar_start=0
        self.progress_bar_step=0


        # GUI setup
        self.window=QtGui.QWidget()
#        self.window=QtGui.QMainWindow()
        
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
                        'Initial Model',
                        'PANDDAs',
                        'Summary',
                        'Queue Control',
                        'Settings'  ]
        self.tab_dict={}
        for page in tab_list:
            tab=QtGui.QWidget()
            vbox=QtGui.QVBoxLayout(tab)
            tab_widget.addTab(tab,page)
            self.tab_dict[page]=[tab,vbox]

        # Mounted Crystals Tab
        self.mounted_crystals_widget=QtGui.QWidget()
        self.tab_dict['Mounted Crystals'][1].addWidget(self.mounted_crystals_widget)
        get_mounted_crystals_button=QtGui.QPushButton("\nLoad Samples\n")
        get_mounted_crystals_button.clicked.connect(self.button_clicked)
        self.tab_dict['Mounted Crystals'][1].addWidget(get_mounted_crystals_button)


        # DLS @ Data Collection Tab
        # main vbox in tab: self.tab_dict['DLS @ Data Collection'][1]
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


#        self.data_collection_widget=QtGui.QWidget()
#        self.tab_dict['DLS @ Data Collection'][1].addWidget(self.data_collection_widget)
#        self.data_collection_vbox_for_table=QtGui.QVBoxLayout()
#        self.data_collection_scrolled_area=QtGui.QScrollArea()
#        self.data_collection_vbox_for_table.setSizeConstraint(QtGui.QLayout.SetMinAndMaxSize)

#        self.tab_dict['DLS @ Data Collection'][1].addLayout(self.data_collection_vbox_for_table)
#        self.data_collection_scrolled_area=QtGui.QScrollArea()
#        self.data_collection_scrolled_area.setWidget(self.tab_dict['DLS @ Data Collection'][1])

#        self.tab_dict['DLS @ Data Collection'][1].addLayout(self.data_collection_vbox_for_table)
#        data_collection_vbox.addLayout(self.data_collection_vbox_for_table)

#        data_collection_button_hbox=QtGui.QHBoxLayout()
#        get_data_collection_button=QtGui.QPushButton("Get New Results from Autoprocessing")
#        get_data_collection_button.clicked.connect(self.button_clicked)
#        data_collection_button_hbox.addWidget(get_data_collection_button)
#        write_files_button=QtGui.QPushButton("Save Files from Autoprocessing in 'inital_model' Folder")
#        data_collection_button_hbox.addWidget(write_files_button)
#        data_collection_vbox.addLayout(data_collection_button_hbox)
#        self.tab_dict['DLS @ Data Collection'][1].addLayout(data_collection_button_hbox)


#        # DLS @ Data Collection Tab
#        data_collection_vbox=QtGui.QVBoxLayout()
#
#        self.data_collection_vbox_for_table=QtGui.QVBoxLayout()
##        self.data_collection_scrolled_area=QtGui.QScrollArea()
#        self.data_collection_vbox_for_table.setSizeConstraint(QtGui.QLayout.SetMinAndMaxSize)
#
##        self.tab_dict['DLS @ Data Collection'][1].addLayout(self.data_collection_vbox_for_table)
##        self.data_collection_scrolled_area=QtGui.QScrollArea()
##        self.data_collection_scrolled_area.setWidget(self.tab_dict['DLS @ Data Collection'][1])
#
##        self.tab_dict['DLS @ Data Collection'][1].addLayout(self.data_collection_vbox_for_table)
#        data_collection_vbox.addLayout(self.data_collection_vbox_for_table)
#
#        data_collection_button_hbox=QtGui.QHBoxLayout()
#        get_data_collection_button=QtGui.QPushButton("Get New Results from Autoprocessing")
#        get_data_collection_button.clicked.connect(self.button_clicked)
#        data_collection_button_hbox.addWidget(get_data_collection_button)
#        write_files_button=QtGui.QPushButton("Save Files from Autoprocessing in 'inital_model' Folder")
#        data_collection_button_hbox.addWidget(write_files_button)
#        data_collection_vbox.addLayout(data_collection_button_hbox)
#        self.tab_dict['DLS @ Data Collection'][1].addLayout(data_collection_vbox)



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
                self.work_thread=read_autoprocessing_results_from_disc(self.visit_list,self.target)
            self.explorer_active=1
            self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.connect(self.work_thread, QtCore.SIGNAL("create_widgets_for_autoprocessing_results"),
                                                     self.create_widgets_for_autoprocessing_results)
            self.work_thread.start()
### -------------------------------------------------------------------
            if self.sender().text()=="Save Files from Autoprocessing in 'inital_model' Folder":
                self.work_thread=save_autoprocessing_results_to_disc(self.dataset_outcome_dict)
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

    def create_widgets_for_autoprocessing_results(self,dict_list):
        data_collection_dict=dict_list[0]
        data_collection_statistics_dict=dict_list[1]

### --- used temporarily to be able to test stuff offline ---
        pickle.dump(data_collection_dict,open('data_collection_dict.p','wb'))
        pickle.dump(data_collection_statistics_dict,open('data_collection_statistics_dict.p','wb'))
#        data_collection_dict = pickle.load( open( "/usr/local/scripts/tobias/XChemExplorer/tmp/data_collection_dict.p", "rb" ) )
#        data_collection_statistics_dict= pickle.load( open( "/usr/local/scripts/tobias/XChemExplorer/tmp/data_collection_statistics_dict.p", "rb" ) )
### ---------------------------------------------------------
        fh = open("imageToSave.png", "wb")
        fh.write(x.decode('base64'))
        fh.close()


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
                                        'Multiplicity\nOuter Shell' ]



#        scroll = QtGui.QScrollArea()
#        scroll.setWidget(mygroupbox)
#        scroll.setWidgetResizable(True)
#        scroll.setFixedHeight(400)
#        layout = QtGui.QVBoxLayout(self)
#        layout.addWidget(scroll)

#        data_collection_vbox_for_table=QtGui.QVBoxLayout()

#        scroll=QtGui.QScrollArea()
#        scroll.setWidgetResizable(True)
#        self.data_collection_vbox_for_table.addWidget(s)

#        w=QtGui.QWidget()
#        vbox=QtGui.QVBoxLayout(w)
#        vbox_table=QtGui.QVBoxLayout()

#        for key in data_collection_dict:
#            if key=='ATAD2A-x383':
#                layout = QtGui.QGridLayout()
#                print len(data_collection_dict[key][1])
#                for run_number,run in enumerate(data_collection_dict[key][0]):
#                    label = QtGui.QLabel(run[0])
#                    layout.addWidget(label,(run_number)*2,0)
#                    # this nifty expression return how often a certain string appears in a list
#                    for column_number,column in enumerate(sorted(filter(lambda x: run[0] in x,data_collection_dict[key][1]))):
#                        print column
#                        pixmap = QtGui.QPixmap(column)
#                        label = QtGui.QLabel()
#                        label.resize(320,200)
#                        label.setPixmap(pixmap.scaled(label.size(), QtCore.Qt.KeepAspectRatio))
#                        layout.addWidget(label, (run_number)*2+1, column_number)
#                self.data_collection_vbox_for_table.addLayout(layout)

#self.myLabel.setPixmap(QtGui.QPixmap(_fromUtf8(directory + '\\' + tempName)).scaled(self.myLabel.size(), QtCore.Qt.KeepAspectRatio))

#        pixmap = QtGui.QPixmap(path)
#        layout = QtGui.QGridLayout(self)
#        for row in range(4):
#            for column in range(4):
#                label = QtGui.QLabel(self)
#                label.setPixmap(pixmap)
#                layout.addWidget(label, row, column)








        # OK this sort of works ------------------------------------------------------------------------
#        scrolled_window=QtGui.QWidget()
#        data_collection_vbox_for_table=QtGui.QVBoxLayout(scrolled_window)
#        i=0
#        for key in sorted(data_collection_statistics_dict):
#            while i<10:
#                table=QtGui.QTableWidget()
#                table.setSortingEnabled(True)
#                table.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
#                table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
#                table.setSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.MinimumExpanding)
#                table.resizeRowsToContents()
#                table.setRowCount(len(data_collection_statistics_dict[key]))
#                table.setColumnCount(len(data_collection_statistics_dict[key][1])-1)
#                for row,line in enumerate(data_collection_statistics_dict[key]):
#                    for column,item in enumerate(line[1:]):
#                        cell_text=QtGui.QTableWidgetItem()
#                        cell_text.setText(str(item))
#                        table.setItem(row, column, cell_text)
#                table.setHorizontalHeaderLabels(diffraction_data_column_name)
#                data_collection_vbox_for_table.addWidget(table)
#                i+=1
#        self.data_collection_scrolled_area.setWidget(scrolled_window)
        #-----------------------------------------------------------------------------------------------

        dataset_outcome = { "success":                      "rgb(0,255,0)",
                            "Failed - centring failed":     "rgb(255,0,0)",
                            "Failed - no diffraction":      "rgb(255,0,0)",
                            "Failed - processing barfs":    "rgb(255,0,0)",
                            "Failed - loop empty":          "rgb(255,0,0)",
                            "Failed - low resolution":      "rgb(255,0,0)",
                            "Failed - no X-rays":           "rgb(255,0,0)",
                            "Failed - unknown":             "rgb(255,0,0)"  }


        # And this is another test ------------------------------------------------------------------------
        table=QtGui.QTableWidget()
        table.setSortingEnabled(True)
#        table.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
#        table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
#        table.setSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.MinimumExpanding)
#        table.resizeColumnsToContents()
        table.setRowCount(len(data_collection_dict))
        table.setColumnCount(3)
        self.dataset_outcome_dict={}        # not sure if wise, but this resets any previous entries
        for row,key in enumerate(sorted(data_collection_statistics_dict)):
            self.dataset_outcome_dict[key]=[]
            # this is the main table
            table.setItem(row, 0, QtGui.QTableWidgetItem(key))
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
                # this nifty expression return how often a certain string appears in a list
                for column_number,column in enumerate(sorted(filter(lambda x: run[0] in x,data_collection_dict[key][1]))):
                    pixmap = QtGui.QPixmap(column)
                    label = QtGui.QLabel()
                    label.resize(320,200)
                    label.setPixmap(pixmap.scaled(label.size(), QtCore.Qt.KeepAspectRatio))
                    layout.addWidget(label, (run_number)*2+1, column_number)
            vbox_cell.addLayout(layout)


            hbox_for_button_and_table=QtGui.QHBoxLayout()
            # dataset outcome buttons
            dataset_outcome_groupbox=QtGui.QGroupBox()
            dataset_outcome_vbox=QtGui.QVBoxLayout()
            for outcome in sorted(dataset_outcome):
                button=QtGui.QPushButton(outcome)
                button.setAutoExclusive(True)
                button.setCheckable(True)
                button.setStyleSheet("background-color: "+dataset_outcome[outcome])
                button.clicked.connect(self.dataset_outcome_button_change_color)
                self.dataset_outcome_dict[key].append(button)
                if outcome=='success':
                    button.setChecked(True)
                dataset_outcome_vbox.addWidget(button)

            dataset_outcome_groupbox.setLayout(dataset_outcome_vbox)
            hbox_for_button_and_table.addWidget(dataset_outcome_groupbox)

            # table for data processing results
            data_collection_table=QtGui.QTableWidget()
            data_collection_table.setRowCount(len(data_collection_statistics_dict[key]))
            data_collection_table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
            data_collection_table.setColumnCount(len(data_collection_statistics_dict[key][0])-1)
            for n,line in enumerate(data_collection_statistics_dict[key]):
                for column,item in enumerate(line[1:]):
                    cell_text=QtGui.QTableWidgetItem()
                    cell_text.setText(str(item))
                    data_collection_table.setItem(n, column, cell_text)
            data_collection_table.setHorizontalHeaderLabels(diffraction_data_column_name)
            data_collection_table.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
            hbox_for_button_and_table.addWidget(data_collection_table)
            vbox_cell.addLayout(hbox_for_button_and_table)

            cell_widget.setLayout(vbox_cell)
            table.setCellWidget(row, 2, cell_widget)
#            table.stretchLastSection()
            table.setColumnWidth(2, 1000)

#            table.setItem(row, 2, vbox)
        table.setHorizontalHeaderLabels(['Sample','Date','Stuff'])
#        table.resizeColumnsToContents()
        table.resizeRowsToContents()
        table.setLineWidth(10)
        self.data_collection_vbox_for_table.addWidget(table)
        #-----------------------------------------------------------------------------------------------



#            vbox_table.addWidget(table)
#            scroll.setWidget(table)
#        s.setWidget(w)
#        s.addWidget(vbox_table)
#        self.data_collection_vbox_for_table.addWidget(scroll)
#        scroll.setWidget(self.data_collection_widget)
#        self.data_collection_widget.setLayout(vbox_table)
#        self.data_collection_scrolled_area.setWidget(data_collection_vbox_for_table)




#                label.addWidget(table)
#                print data_collection_dict[key]
#                print data_collection_statistics_dict[key][1:]
#        hbox_image=QtGui.QHBoxLayout()
#        for image in data_collection_dict['ATAD2A-x367'][1]:
#            print image
#            pixmap = QtGui.QPixmap(image)
#            label = QtGui.QLabel()
##            pixmap_scaled = pixmap.scaled(label.size(), Qt.KeepAspectRatio)
#
#            label.setPixmap(pixmap)
#            hbox_image.addWidget(label)
#        self.data_collection_vbox_for_table.addLayout(hbox_image)

    def dataset_outcome_button_change_color(self):
        print self.sender().text()

class save_autoprocessing_results_to_disc(QtCore.QThread):
    def __init__(self,dataset_outcome_dict):
        QtCore.QThread.__init__(self)
        self.dataset_outcome_dict=dataset_outcome_dict

    def run(self):
        for key in sorted(self.dataset_outcome_dict):
            for button in self.dataset_outcome_dict[key]:
                if button.isChecked():
                    print key
                    print key,button.text()
        self.emit(QtCore.SIGNAL("finished()"))


class read_autoprocessing_results_from_disc(QtCore.QThread):
    def __init__(self,visit_list,target):
        QtCore.QThread.__init__(self)
        self.visit_list=visit_list
        self.target=target
        self.data_collection_dict={}
        self.data_collection_statistics_dict={}

#        self.data_collection_run_dict={}
#        self.data_collection_image_dict={}
#        self.data_collection_logfile_dict={}

    def run(self):
        for visit_directory in sorted(self.visit_list):
            progress_step=100/float(len(glob.glob(os.path.join(visit_directory,'processed',self.target,'*'))))
            progress=0
            visit=visit_directory.split('/')[5]
            for collected_xtals in sorted(glob.glob(os.path.join(visit_directory,'processed',self.target,'*'))):
                xtal=collected_xtals[collected_xtals.rfind('/')+1:]
                self.data_collection_dict[xtal]=[[],[],[],[]]
                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), xtal)
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
                    if image.endswith('thumb.jpeg'):
                        image_list.append(image)
                    if image.endswith('_.png'):
                        image_list.append(image)
                self.data_collection_dict[xtal][1]+=image_list
                self.data_collection_dict[xtal][2]+=logfile_list
                progress += progress_step
                self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

                # convert images to strong and attach to
                tmp=[]
                for image in self.data_collection_dict[xtal][1]:
                    print image[image.rfind('/')+1:]+'   '+image
                    image_file=open(image,"rb")
                    image_string=base64.b64encode(image_file.read())
                    image_string_list.append((image[image.rfind('/')+1:],image_string))
                self.data_collection_dict[xtal][3]+=image_string_list



        progress_step=100/float(len(self.data_collection_dict))
        progress=0
        for key in sorted(self.data_collection_dict):
            self.data_collection_statistics_dict[key]=[]
            if not self.data_collection_dict[key][2]==[]:
                for logfile in self.data_collection_dict[key][2]:
                    self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'parsing aimless logfiles of '+key)
                    aimless_results=parse().GetAimlessLog(logfile)
                    self.data_collection_statistics_dict[key].append([
                                logfile,
                                aimless_results['AutoProc'],
                                aimless_results['Run'],
                                aimless_results['SpaceGroup'],
                                aimless_results['UnitCell'],
                                aimless_results['ResolutionLow']+'-'+aimless_results['ResolutionHigh'],
                                aimless_results['ResolutionLow']+'-'+aimless_results['ResolutionLowInnerShell'],
                                aimless_results['ResolutionHighOuterShell']+'-'+aimless_results['ResolutionHigh'],
                                aimless_results['RmergeOverall'],
                                aimless_results['RmergeLow'],
                                aimless_results['RmergeHigh'],
                                aimless_results['IsigOverall'],
                                aimless_results['IsigLow'],
                                aimless_results['IsigHigh'],
                                aimless_results['CompletenessOverall'],
                                aimless_results['CompletenessLow'],
                                aimless_results['CompletenessHigh'],
                                aimless_results['MultiplicityOverall'],
                                aimless_results['MultiplicityLow'],
                                aimless_results['MultiplicityHigh'] ])
            else:
                self.data_collection_statistics_dict[key]+='###'*20
            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

#            print self.data_collection_dict[key]
#        print self.data_collection_dict['ATAD2A-x367']
#        print ''
#        print self.data_collection_dict['ATAD2A-x367'][0]
#        print ''
#        print self.data_collection_dict['ATAD2A-x367'][1]
#        print ''
#        print self.data_collection_dict['ATAD2A-x367'][2]
#        return


#        hbox_status.addWidget(self.status_bar)
#        hbox_status.addWidget(self.progress_bar)

        self.emit(QtCore.SIGNAL('create_widgets_for_autoprocessing_results'), [self.data_collection_dict,
                                                                               self.data_collection_statistics_dict])
#myPixmap = QtGui.QPixmap(_fromUtf8('image.jpg'))
#myScaledPixmap = myPixmap.scaled(self.label.size(), Qt.KeepAspectRatio)
#self.label.setPixmap(myScaledPixmap)


#def main():
#    app = QApplication(sys.argv)
#    w = MyWindow()
#    w.show()
#    sys.exit(app.exec_())

#class MyWindow(QWidget):
#    def __init__(self, *args):
#        QWidget.__init__(self, *args)
#
#        tablemodel = MyTableModel(my_array, self)
#        tableview = QTableView()
#        tableview.setModel(tablemodel)
#
#        layout = QVBoxLayout(self)
#        layout.addWidget(tableview)
#        self.setLayout(layout)


if __name__ == "__main__":
    app=XChemExplorer(sys.argv)

#app = QtGui.QApplication(sys.argv)
#MainWindow()
##frame.show()
#sys.exit(app.exec_())  
