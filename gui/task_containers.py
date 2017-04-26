from import_modules import *
import datasourcebutton_functions as buttonfunc
import TCdatasets_functions as TCdatafunc

def workflow(self):
    self.workflow =             [   'Overview',         # 0
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
                else:
                    self.need_to_switch_main_tab(task_index)

def datasource_button(self):
    update_from_datasource_button = QtGui.QPushButton("Update Tables\nFrom Datasource")
    update_from_datasource_button.setToolTip(XChemToolTips.update_from_datasource_button_tip())
    update_from_datasource_button.setStyleSheet(
        "QPushButton { padding: 1px; margin: 1px; background: rgb(140,140,140) }")
    self.headlineLabelfont = QtGui.QFont("Arial", 20, QtGui.QFont.Bold)
    update_from_datasource_button.setFont(self.headlineLabelfont)
    update_from_datasource_button.clicked.connect(buttonfunc.datasource_menu_reload_samples)
    update_from_datasource_button.setMaximumWidth((self.screen.width() - 20) / 5)

    return update_from_datasource_button

def datasets_section(self):
    self.dataset_tasks = ['Get New Results from Autoprocessing',
                          'Run DIMPLE on All Autoprocessing MTZ files',
                          'Rescore Datasets',
                          'Read PKL file',
                          'Run xia2 on selected datasets',
                          'Run xia2 on selected datasets - overwrite']

    frame_dataset_task = QtGui.QFrame()
    frame_dataset_task.setFrameShape(QtGui.QFrame.StyledPanel)
    frame_dataset_task.setStyleSheet(
        "QFrame { border: 1px solid black; border-radius: 1px; padding: 0px; margin: 0px }")
    vboxTask = QtGui.QVBoxLayout()
    label = QtGui.QLabel('Datasets')
    label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
    label.setFont(self.headlineLabelfont)
    label.setStyleSheet(
        " QLabel { border: 1px solid black; border-radius: 1px; background: rgb(240,255,140); padding: 0px; margin: 0px }")
    vboxTask.addWidget(label)
    hboxAction = QtGui.QHBoxLayout()
    self.dataset_tasks_combobox = QtGui.QComboBox()
    for task in self.dataset_tasks:
        self.dataset_tasks_combobox.addItem(task)
    self.dataset_tasks_combobox.setToolTip(XChemToolTips.dataset_task_tip())
    self.dataset_tasks_combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")
    hboxAction.addWidget(self.dataset_tasks_combobox)
    vboxButton = QtGui.QVBoxLayout()
    dataset_task_run_button = QtGui.QPushButton("Run")
    dataset_task_run_button.setToolTip(XChemToolTips.dataset_task_run_button_tip())
    dataset_task_run_button.clicked.connect(button_clicked)
    dataset_task_run_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
    vboxButton.addWidget(dataset_task_run_button)
    dataset_task_status_button = QtGui.QPushButton("Status")
    dataset_task_status_button.setToolTip(XChemToolTips.dataset_task_status_button_tip())
    dataset_task_status_button.clicked.connect(button_clicked)
    dataset_task_status_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
    vboxButton.addWidget(dataset_task_status_button)
    hboxAction.addLayout(vboxButton)
    vboxTask.addLayout(hboxAction)
    vboxTask.setSpacing(0)
    vboxTask.setMargin(0)
    frame_dataset_task.setLayout(vboxTask)
    frame_dataset_task.setMaximumWidth((self.screen.width() - 20) / 5)

    self.workflow_widget_dict['Datasets'] = [self.dataset_tasks_combobox, dataset_task_run_button,
                                             dataset_task_status_button]

    return frame_dataset_task

def maps_section(self):
    self.map_cif_file_tasks = ['Run DIMPLE on selected MTZ files',
                               'Remove selected DIMPLE PDB/MTZ files',
                               'Create CIF/PDB/PNG file of ALL compounds',
                               'Create CIF/PDB/PNG file of NEW compounds',
                               'Create CIF/PDB/PNG file of SELECTED compounds']

    frame_map_cif_file_task = QtGui.QFrame()
    frame_map_cif_file_task.setFrameShape(QtGui.QFrame.StyledPanel)
    frame_map_cif_file_task.setStyleSheet(
        "QFrame { border: 1px solid black; border-radius: 1px; padding: 0px; margin: 0px }")
    vboxTask = QtGui.QVBoxLayout()
    #       label=QtGui.QLabel(self.workflow_dict['Maps'])
    label = QtGui.QLabel('Maps & Restraints')
    label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
    label.setFont(self.headlineLabelfont)
    label.setStyleSheet(
        " QLabel { border: 1px solid black; border-radius: 1px; background: rgb(140,255,150); padding: 0px; margin: 0px }")
    vboxTask.addWidget(label)
    hboxAction = QtGui.QHBoxLayout()
    self.map_cif_file_tasks_combobox = QtGui.QComboBox()
    for task in self.map_cif_file_tasks:
        self.map_cif_file_tasks_combobox.addItem(task)
    self.map_cif_file_tasks_combobox.setToolTip(XChemToolTips.map_cif_file_task_tip())
    self.map_cif_file_tasks_combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")
    hboxAction.addWidget(self.map_cif_file_tasks_combobox)
    vboxButton = QtGui.QVBoxLayout()
    map_cif_file_task_run_button = QtGui.QPushButton("Run")
    map_cif_file_task_run_button.setToolTip(XChemToolTips.map_cif_file_task_run_button_tip())
    map_cif_file_task_run_button.clicked.connect(button_clicked)
    map_cif_file_task_run_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
    vboxButton.addWidget(map_cif_file_task_run_button)
    map_cif_file_task_status_button = QtGui.QPushButton("Status")
    map_cif_file_task_status_button.setToolTip(XChemToolTips.map_cif_file_task_status_button_tip())
    map_cif_file_task_status_button.clicked.connect(button_clicked)
    map_cif_file_task_status_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
    vboxButton.addWidget(map_cif_file_task_status_button)
    hboxAction.addLayout(vboxButton)
    vboxTask.addLayout(hboxAction)
    vboxTask.setSpacing(0)
    vboxTask.setMargin(0)
    frame_map_cif_file_task.setLayout(vboxTask)

    self.workflow_widget_dict['Maps'] = [self.map_cif_file_tasks_combobox, map_cif_file_task_run_button,
                                         map_cif_file_task_status_button]
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
    label.setFont(self.headlineLabelfont)
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
    map_cif_file_task_run_button.clicked.connect(button_clicked)
    map_cif_file_task_run_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
    vboxButton.addWidget(map_cif_file_task_run_button)
    map_cif_file_task_status_button=QtGui.QPushButton("Status")
    map_cif_file_task_status_button.setToolTip(XChemToolTips.map_cif_file_task_status_button_tip())
    map_cif_file_task_status_button.clicked.connect(button_clicked)
    map_cif_file_task_status_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
    vboxButton.addWidget(map_cif_file_task_status_button)
    hboxAction.addLayout(vboxButton)
    vboxTask.addLayout(hboxAction)
    vboxTask.setSpacing(0)
    vboxTask.setMargin(0)
    frame_map_cif_file_task.setLayout(vboxTask)
    frame_map_cif_file_task.setMaximumWidth((self.screen.width() - 20) / 5)

    self.workflow_widget_dict['Maps']=[self.map_cif_file_tasks_combobox,map_cif_file_task_run_button,map_cif_file_task_status_button]

    return frame_map_cif_file_task

def panddas_section(self):
    self.panddas_file_tasks = ['pandda.analyse',
                               'pandda.inspect',
                               'run pandda.inspect at home',
                               'Export NEW PANDDA models',
                               'Export ALL PANDDA models',
                               'Show HTML summary',
                               'Update datasource with results from pandda.inspect',
                               'cluster datasets',
                               'Event Map -> SF',
                               'check modelled ligands']

    frame_panddas_file_task = QtGui.QFrame()
    frame_panddas_file_task.setFrameShape(QtGui.QFrame.StyledPanel)
    frame_panddas_file_task.setStyleSheet(
        "QFrame { border: 1px solid black; border-radius: 1px; padding: 0px; margin: 0px }")
    vboxTask = QtGui.QVBoxLayout()
    label = QtGui.QLabel('Hit Identification')
    label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
    label.setFont(self.headlineLabelfont)
    label.setStyleSheet(
        " QLabel { border: 1px solid black; border-radius: 1px; background: rgb(140,200,255); padding: 0px; margin: 0px }")
    vboxTask.addWidget(label)
    hboxAction = QtGui.QHBoxLayout()
    self.panddas_file_tasks_combobox = QtGui.QComboBox()
    for task in self.panddas_file_tasks:
        self.panddas_file_tasks_combobox.addItem(task)
    self.panddas_file_tasks_combobox.setToolTip(XChemToolTips.panddas_file_task_tip())
    self.panddas_file_tasks_combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")
    hboxAction.addWidget(self.panddas_file_tasks_combobox)
    vboxButton = QtGui.QVBoxLayout()
    panddas_file_task_run_button = QtGui.QPushButton("Run")
    panddas_file_task_run_button.setToolTip(XChemToolTips.panddas_file_task_run_button_tip())
    panddas_file_task_run_button.clicked.connect(button_clicked)
    panddas_file_task_run_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
    vboxButton.addWidget(panddas_file_task_run_button)
    panddas_file_task_status_button = QtGui.QPushButton("Status")
    panddas_file_task_status_button.setToolTip(XChemToolTips.panddas_file_task_status_button_tip())
    panddas_file_task_status_button.clicked.connect(button_clicked)
    panddas_file_task_status_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
    vboxButton.addWidget(panddas_file_task_status_button)
    hboxAction.addLayout(vboxButton)
    vboxTask.addLayout(hboxAction)
    vboxTask.setSpacing(0)
    vboxTask.setMargin(0)
    frame_panddas_file_task.setLayout(vboxTask)
    frame_panddas_file_task.setMaximumWidth((self.screen.width() - 20) / 5)

    self.workflow_widget_dict['PANDDAs'] = [self.panddas_file_tasks_combobox, panddas_file_task_run_button,
                                            panddas_file_task_status_button]

    return frame_panddas_file_task

def refine_section(self):
    self.refine_file_tasks = ['Open COOT',
                              'Open COOT - new interface',
                              'Update Deposition Table',
                              'Prepare Group Deposition']

    frame_refine_file_task = QtGui.QFrame()
    frame_refine_file_task.setFrameShape(QtGui.QFrame.StyledPanel)
    frame_refine_file_task.setStyleSheet(
        "QFrame { border: 1px solid black; border-radius: 1px; padding: 0px; margin: 0px }")
    vboxTask = QtGui.QVBoxLayout()
    #        label=QtGui.QLabel(self.workflow_dict['Refinement'])
    label = QtGui.QLabel('Refinement')
    label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
    label.setFont(self.headlineLabelfont)
    label.setStyleSheet(
        " QLabel { border: 1px solid black; border-radius: 1px; background: rgb(245,190,255); padding: 0px; margin: 0px }")
    vboxTask.addWidget(label)
    hboxAction = QtGui.QHBoxLayout()
    self.refine_file_tasks_combobox = QtGui.QComboBox()
    for task in self.refine_file_tasks:
        self.refine_file_tasks_combobox.addItem(task)
    self.refine_file_tasks_combobox.setToolTip(XChemToolTips.refine_file_task_tip())
    self.refine_file_tasks_combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")
    hboxAction.addWidget(self.refine_file_tasks_combobox)
    vboxButton = QtGui.QVBoxLayout()
    refine_file_task_run_button = QtGui.QPushButton("Run")
    refine_file_task_run_button.setToolTip(XChemToolTips.refine_file_task_run_button_tip())
    refine_file_task_run_button.clicked.connect(button_clicked)
    refine_file_task_run_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
    vboxButton.addWidget(refine_file_task_run_button)
    refine_file_task_status_button = QtGui.QPushButton("Status")
    refine_file_task_status_button.setToolTip(XChemToolTips.refine_file_task_status_button_tip())
    refine_file_task_status_button.clicked.connect(button_clicked)
    refine_file_task_status_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
    vboxButton.addWidget(refine_file_task_status_button)
    hboxAction.addLayout(vboxButton)
    vboxTask.addLayout(hboxAction)
    vboxTask.setSpacing(0)
    vboxTask.setMargin(0)
    frame_refine_file_task.setLayout(vboxTask)
    frame_refine_file_task.setMaximumWidth((self.screen.width()) - 20 / 5)

    self.workflow_widget_dict['Refinement'] = [self.refine_file_tasks_combobox, refine_file_task_run_button,
                                               refine_file_task_status_button]

    return frame_refine_file_task

def setup_for_layout(self):
    # build task container elements
    workflow(self)
    datasource = datasource_button(self)
    datasets = datasets_section(self)
    maps = maps_section(self)
    panddas = panddas_section(self)
    refine = refine_section(self)

    self.status_bar=QtGui.QStatusBar()
    self.progress_bar=QtGui.QProgressBar()
    self.progress_bar.setMaximum(100)
    self.status_bar.setMaximumWidth(self.screen.width())
    self.progress_bar.setMaximumWidth(self.screen.width())

    self.hbox_status = QtGui.QHBoxLayout()
    self.hbox_status.addWidget(self.status_bar)
    self.hbox_status.addWidget(self.progress_bar)

    self.hboxTaskFrames = QtGui.QHBoxLayout()


    self.hboxTaskFrames.addWidget(datasource)
    self.hboxTaskFrames.addWidget(datasets)
    self.hboxTaskFrames.addWidget(maps)
    self.hboxTaskFrames.addWidget(panddas)
    self.hboxTaskFrames.addWidget(refine)