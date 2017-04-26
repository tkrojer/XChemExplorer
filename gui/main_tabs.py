from import_modules import *
import dataset_functions as datafunc
import maps_functions as mapfunc
import pandda_functions as panfunc
import settings_functions as settingsfunc
import deposition_functions as depofunc

def add_subtabs(self, tablist, tabwidget):
    tabdict = {}
    for page in tablist:
        tab = QtGui.QWidget()
        vbox = QtGui.QVBoxLayout(tab)
        tabwidget.addTab(tab, page)
        tabdict[page] = [tab, vbox]
    return tabdict

def overview_tab(self):
    ## main tab layout
    # add overview tab to layout
    overview_tab_widget = QtGui.QTabWidget()
    self.tab_dict[self.workflow_dict['Overview']][1].addWidget(overview_tab_widget)

    ## sub tab layout
    # Set names for sub tabs in data source tab
    overview_tab_list = ['Data Source',
                         'Summary']
    # create vbox containers for sub tabs
    self.overview_tab_dict = add_subtabs(self, overview_tab_list, overview_tab_widget)

    ## Data Source tab
    # set columns to automatically display in data source sub tab
    self.data_source_columns_to_display = ['Sample ID',
                                           'Compound ID',
                                           'Smiles',
                                           'Visit',
                                           'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                           'Refinement\nRfree',
                                           'Data Collection\nDate',
                                           'Puck',
                                           'PuckPosition',
                                           'Ligand\nConfidence']

    # setup table to be displayed
    self.mounted_crystal_table = QtGui.QTableWidget()
    self.mounted_crystal_table.setSortingEnabled(True)
    self.mounted_crystal_table.resizeColumnsToContents()

    # add data source tab to layout
    self.overview_tab_dict['Data Source'][1].addWidget(self.mounted_crystal_table)

    ## Overview tab
    # set up and plot graph
    self.overview_figure, self.overview_axes = plt.subplots()
    self.overview_canvas = FigureCanvas(self.overview_figure)
    self.update_summary_plot()

    # add overview tab to layout
    self.overview_tab_dict['Summary'][1].addWidget(self.overview_canvas)

def dataset_tab(self):
    ## Dataset tab
    # add dataset tab to layout
    dls_tab_widget = QtGui.QTabWidget()

    # set names for subtabs in dls tab
    dls_tab_list = ['Summary',
                    'Reprocess']
    # create vbox containers for sub tabs
    self.dls_tab_dict = add_subtabs(self, dls_tab_list, dls_tab_widget)

    ## Data Collection tab (which is actually the Datasets main tab, but OK)
    self.dls_data_collection_vbox = QtGui.QVBoxLayout()
    self.tab_dict[self.workflow_dict['Datasets']][1].addLayout(self.dls_data_collection_vbox)

    # add checkbox to check for new data collection
    hbox = QtGui.QHBoxLayout()
    check_for_new_data_collection = QtGui.QCheckBox('Check for new data collection every two minutes')
    check_for_new_data_collection.toggle()
    check_for_new_data_collection.setChecked(False)
    check_for_new_data_collection.stateChanged.connect(datafunc.collection_autoupdate)
    hbox.addWidget(check_for_new_data_collection)

    # add dropdown for target selection
    hbox.addWidget(QtGui.QLabel('                                             '))
    hbox.addWidget(QtGui.QLabel('Select Target: '))
    self.target_selection_combobox = QtGui.QComboBox()
    self.populate_target_selection_combobox(self.target_selection_combobox)
    self.target_selection_combobox.activated[str].connect(datafunc.target_selection_combobox_activated)
    hbox.addWidget(self.target_selection_combobox)
    self.target = str(self.target_selection_combobox.currentText())
    self.dls_data_collection_vbox.addLayout(hbox)

    ## Summary sub tab
    # columns to display in table
    self.data_collection_summary_column_name = ['Sample ID',
                                                # 'Date',
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
                                                'Show Diffraction\nImage']

    # setup summary table
    self.data_collection_summary_table = QtGui.QTableWidget()
    self.data_collection_summary_table.setColumnCount(len(self.data_collection_summary_column_name))
    self.data_collection_summary_table.setSortingEnabled(True)
    self.data_collection_summary_table.setHorizontalHeaderLabels(self.data_collection_summary_column_name)

    # add vbox to insert table into to layout
    self.data_collection_summarys_vbox_for_table = QtGui.QVBoxLayout()
    self.dls_tab_dict['Summary'][1].addLayout(self.data_collection_summarys_vbox_for_table)

    # insert table into vbox and add to layout
    self.data_collection_summarys_vbox_for_details = QtGui.QVBoxLayout()
    self.data_collection_details_currently_on_display = None
    self.dls_tab_dict['Summary'][1].addLayout(self.data_collection_summarys_vbox_for_details)
    self.data_collection_summarys_vbox_for_table.addWidget(self.data_collection_summary_table)
    self.dls_data_collection_vbox.addWidget(dls_tab_widget)

    ## Dewar sub-tab
    ## WARNING: CURRENTLY UNUSED

    self.dewar_configuration_dict = {}
    self.dewar_sample_configuration_dict = {}
    self.dewar_label_active = ''
    self.dewar_configuration_layout = QtGui.QGridLayout()

    # create context menu
    self.popMenu = QtGui.QMenu()
    recollect = QtGui.QAction("recollect", self.window)
    recollect.triggered.connect(datafunc.flag_sample_for_recollection)
    undo_recollect = QtGui.QAction("undo", self.window)
    undo_recollect.triggered.connect(datafunc.undo_flag_sample_for_recollection)
    self.popMenu.addAction(recollect)
    self.popMenu.addAction(undo_recollect)

    for puck in range(38):
        for position in range(17):
            frame = QtGui.QFrame()
            frame.setFrameShape(QtGui.QFrame.StyledPanel)
            vbox_for_frame = QtGui.QVBoxLayout()
            if puck == 0 and position == 0:
                label = QtGui.QLabel('')
                vbox_for_frame.addWidget(label)
                frame.setLayout(vbox_for_frame)
            elif puck == 0 and position != 0:
                label = QtGui.QLabel(str(position))
                vbox_for_frame.addWidget(label)
                frame.setLayout(vbox_for_frame)
            elif position == 0 and puck != 0:
                label = QtGui.QLabel(str(puck))
                vbox_for_frame.addWidget(label)
                frame.setLayout(vbox_for_frame)
            else:
                frame = QtGui.QPushButton('')
                frame.setStyleSheet("font-size:5px;border-width: 0px")
                frame.clicked.connect(self.show_html_summary_in_firefox)
                # how to right click on button
                self.dewar_configuration_dict[str(puck) + '-' + str(position)] = frame
                self.dewar_sample_configuration_dict[str(puck) + '-' + str(position)] = []
                frame.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
                frame.customContextMenuRequested.connect(datafunc.on_context_menu)
            self.dewar_configuration_layout.addWidget(frame, position, puck)

    ## END WARNING!

    ## Reprocess sub tab
    reprocess_vbox = QtGui.QVBoxLayout()

    # frame for all options
    options_frame = QtGui.QFrame()
    options_frame.setFrameShape(QtGui.QFrame.StyledPanel)
    hbox = QtGui.QHBoxLayout()
    frame_select = QtGui.QFrame()
    frame_select.setFrameShape(QtGui.QFrame.StyledPanel)

    # layout for data collection directory functions
    hbox_select = QtGui.QHBoxLayout()
    label = QtGui.QLabel('Data collection directory:')
    hbox_select.addWidget(label)

    # frame around directory path
    dir_frame = QtGui.QFrame()
    dir_frame.setFrameShape(QtGui.QFrame.StyledPanel)
    dir_label_box = QtGui.QVBoxLayout()

    # data collection path
    self.diffraction_data_dir_label = QtGui.QLabel(self.diffraction_data_directory) # diffraction_data_directory in:
    dir_label_box.addWidget(self.diffraction_data_dir_label)                        # settings_directories.py
    dir_frame.setLayout(dir_label_box)
    hbox_select.addWidget(dir_frame)

    # 'select' button for data collection directory
    select_button = QtGui.QPushButton("Select")
    select_button.clicked.connect(datafunc.select_diffraction_data_directory)
    hbox_select.addWidget(select_button)
    frame_select.setLayout(hbox_select)
    hbox.addWidget(frame_select)

    # frame for dataset search
    frame_search = QtGui.QFrame()
    frame_search.setFrameShape(QtGui.QFrame.StyledPanel)
    hbox_search = QtGui.QHBoxLayout()

    # search datasets button
    button = QtGui.QPushButton("Search Datasets")
    button.clicked.connect(datafunc.search_for_datasets)
    hbox_search.addWidget(button)

    # frame for number of datasets found
    frame_search_info = QtGui.QFrame()
    frame_search_info.setFrameShape(QtGui.QFrame.StyledPanel)
    hbox_search_info = QtGui.QHBoxLayout()
    self.diffraction_data_search_label = QtGui.QLabel(self.diffraction_data_search_info)
    hbox_search_info.addWidget(self.diffraction_data_search_label)
    frame_search_info.setLayout(hbox_search_info)
    hbox_search.addWidget(frame_search_info)
    frame_search.setLayout(hbox_search)
    hbox.addWidget(frame_search)

    # frame for translation datasetID -> sampleID
    frame_translate = QtGui.QFrame()
    frame_translate.setFrameShape(QtGui.QFrame.StyledPanel)
    vbox_translate = QtGui.QVBoxLayout()

    # label for translate box
    translate_label = QtGui.QLabel('translate:\ndatasetID -> sampleID')
    translate_label.setAlignment(QtCore.Qt.AlignCenter)
    vbox_translate.addWidget(translate_label)

    # button to perform translation
    translate_button = QtGui.QPushButton('Open CSV')
    translate_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
    translate_button.clicked.connect(datafunc.translate_datasetID_to_sampleID)
    vbox_translate.addWidget(translate_button)
    frame_translate.setLayout(vbox_translate)
    hbox.addWidget(frame_translate)
    hbox.addStretch(0)

    # add all frames to main options frame
    options_frame.setLayout(hbox)
    reprocess_vbox.addWidget(frame)

    # columns to  display in reprocess table
    self.reprocess_datasets_column_list = ['Dataset ID',
                                           'Sample ID',
                                           'Run\nxia2',
                                           'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                           'Rmerge\nLow',
                                           'Dimple\nRfree',
                                           'DataProcessing\nSpaceGroup',
                                           'DataProcessing\nUnitCell',
                                           'DataProcessing\nStatus']

    # initiate table and add table to layout
    self.reprocess_datasets_table = QtGui.QTableWidget()
    self.reprocess_datasets_table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
    self.reprocess_datasets_table.setSortingEnabled(True)
    self.reprocess_datasets_table.setColumnCount(len(self.reprocess_datasets_column_list))
    self.reprocess_datasets_table.setHorizontalHeaderLabels(self.reprocess_datasets_column_list)
    reprocess_vbox.addWidget(self.reprocess_datasets_table)

    # create context menu
    self.popMenu_for_reprocess_datasets_table = QtGui.QMenu()
    run_xia2_on_selected = QtGui.QAction("mark selected for reprocessing", self.window)
    run_xia2_on_selected.triggered.connect(self.select_sample_for_xia2)
    self.popMenu_for_reprocess_datasets_table.addAction(run_xia2_on_selected)
    self.reprocess_datasets_table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
    self.reprocess_datasets_table.customContextMenuRequested.connect(datafunc.on_context_menu_reprocess_data)

    ######## GAVE UP ON COMMENTING HERE
    frame = QtGui.QFrame()
    frame.setFrameShape(QtGui.QFrame.StyledPanel)
    hbox = QtGui.QHBoxLayout()

    frame_options = QtGui.QFrame()
    frame_options.setFrameShape(QtGui.QFrame.StyledPanel)
    hbox_options = QtGui.QHBoxLayout()
    label = QtGui.QLabel('Data processing protocol:')
    hbox_options.addWidget(label)
    self.xia2_3d_checkbox = QtGui.QCheckBox(' xia2 3d')
    hbox_options.addWidget(self.xia2_3d_checkbox)
    self.xia2_3dii_checkbox = QtGui.QCheckBox('xia2 3dii')
    hbox_options.addWidget(self.xia2_3dii_checkbox)
    self.xia2_dials_checkbox = QtGui.QCheckBox('Dials')
    hbox_options.addWidget(self.xia2_dials_checkbox)
    frame_options.setLayout(hbox_options)
    hbox.addWidget(frame_options)

    frame_sg = QtGui.QFrame()
    frame_sg.setFrameShape(QtGui.QFrame.StyledPanel)
    hbox_sg = QtGui.QHBoxLayout()
    label = QtGui.QLabel('Space Group:')
    hbox_sg.addWidget(label)
    self.reprocess_space_group_comboxbox = QtGui.QComboBox()
    self.reprocess_space_group_comboxbox.addItem('ignore')
    for sg in XChemMain.space_group_list():
        self.reprocess_space_group_comboxbox.addItem(sg)
    hbox_sg.addWidget(self.reprocess_space_group_comboxbox)
    frame_sg.setLayout(hbox_sg)
    hbox.addWidget(frame_sg)

    frame_ref = QtGui.QFrame()
    frame_ref.setFrameShape(QtGui.QFrame.StyledPanel)
    hbox_ref = QtGui.QHBoxLayout()
    label = QtGui.QLabel('Reference MTZ:')
    hbox_ref.addWidget(label)
    frame_ref_info = QtGui.QFrame()
    frame_ref_info.setFrameShape(QtGui.QFrame.StyledPanel)
    hbox_frame_ref_info = QtGui.QHBoxLayout()
    self.reprocess_reference_mtz_file_label = QtGui.QLabel(self.diffraction_data_reference_mtz)
    hbox_frame_ref_info.addWidget(self.reprocess_reference_mtz_file_label)
    frame_ref_info.setLayout(hbox_frame_ref_info)
    hbox_ref.addWidget(frame_ref_info)
    button = QtGui.QPushButton("Select")
    button.clicked.connect(datafunc.select_reprocess_reference_mtz)
    hbox_ref.addWidget(button)
    frame_ref.setLayout(hbox_ref)
    hbox.addWidget(frame_ref)

    frame_isigma = QtGui.QFrame()
    frame_isigma.setFrameShape(QtGui.QFrame.StyledPanel)
    vbox_isigma = QtGui.QVBoxLayout()
    label = QtGui.QLabel('Resolution\nLimit:\nMn<I/sig(I)>')
    label.setAlignment(QtCore.Qt.AlignCenter)
    vbox_isigma.addWidget(label)
    self.reprocess_isigma_combobox = QtGui.QComboBox()
    misigma = ['default', '4', '3', '2.5', '2', '1.5', '1', '0.5']
    for item in misigma:
        self.reprocess_isigma_combobox.addItem(item)
    self.reprocess_isigma_combobox.setCurrentIndex(0)
    self.reprocess_isigma_combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")
    vbox_isigma.addWidget(self.reprocess_isigma_combobox)
    frame_isigma.setLayout(vbox_isigma)
    hbox.addWidget(frame_isigma)

    frame_cc_half = QtGui.QFrame()
    frame_cc_half.setFrameShape(QtGui.QFrame.StyledPanel)
    vbox_cc_half = QtGui.QVBoxLayout()
    label = QtGui.QLabel('Resolution\nLimit:\nCC 1/2')
    label.setAlignment(QtCore.Qt.AlignCenter)
    vbox_cc_half.addWidget(label)
    self.reprocess_cc_half_combobox = QtGui.QComboBox()
    cc_half = ['default', '0.9', '0.8', '0.7', '0.6', '0.5', '0.4', '0.3', '0.2', '0.1']
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

def maps_tab(self):
    # select all samples in maps tab checkbox
    initial_model_checkbutton_hbox = QtGui.QHBoxLayout()
    select_sample_for_dimple = QtGui.QCheckBox('(de-)select all samples for DIMPLE')
    select_sample_for_dimple.toggle()
    select_sample_for_dimple.setChecked(False)
    select_sample_for_dimple.stateChanged.connect(mapfunc.set_run_dimple_flag)
    initial_model_checkbutton_hbox.addWidget(select_sample_for_dimple)

    set_new_reference_button = QtGui.QPushButton("Set New Reference (if applicable)")
    set_new_reference_button.clicked.connect(mapfunc.set_new_reference_if_applicable)
    initial_model_checkbutton_hbox.addWidget(set_new_reference_button)

    self.reference_file_list = self.get_reference_file_list(' ')
    self.reference_file_selection_combobox = QtGui.QComboBox()
    self.populate_reference_combobox(self.reference_file_selection_combobox)
    initial_model_checkbutton_hbox.addWidget(self.reference_file_selection_combobox)

    self.tab_dict[self.workflow_dict['Maps']][1].addLayout(initial_model_checkbutton_hbox)
    self.initial_model_vbox_for_table = QtGui.QVBoxLayout()
    self.inital_model_column_list = ['Sample ID',
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
                                     'LastUpdated']

    self.initial_model_table = QtGui.QTableWidget()
    self.initial_model_table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
    self.initial_model_table.setSortingEnabled(True)
    self.initial_model_table.setColumnCount(len(self.inital_model_column_list))
    self.initial_model_table.setHorizontalHeaderLabels(self.inital_model_column_list)
    self.initial_model_vbox_for_table.addWidget(self.initial_model_table)
    self.tab_dict[self.workflow_dict['Maps']][1].addLayout(self.initial_model_vbox_for_table)

    # create context menu
    self.popMenu_for_initial_model_table = QtGui.QMenu()
    run_dimple = QtGui.QAction("mark selected for dimple run", self.window)
    run_dimple.triggered.connect(self.select_sample_for_dimple)
    self.popMenu_for_initial_model_table.addAction(run_dimple)
    self.initial_model_table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
    self.initial_model_table.customContextMenuRequested.connect(mapfunc.on_context_menu_initial_model)

def panddas_tab(self):
    pandda_tab_widget = QtGui.QTabWidget()
    # set names for subtabs in panddas tab
    pandda_tab_list = ['pandda.analyse',
                       'Dataset Summary',
                       'Results Summary',
                       'Inspect Summary']
    # create vbox containers for sub tabs
    self.pandda_tab_dict = add_subtabs(self, pandda_tab_list, pandda_tab_widget)

    self.panddas_results_vbox = QtGui.QVBoxLayout()
    self.tab_dict[self.workflow_dict['PANDDAs']][1].addLayout(self.panddas_results_vbox)

    self.pandda_analyse_hbox = QtGui.QHBoxLayout()
    self.pandda_tab_dict['pandda.analyse'][1].addLayout(self.pandda_analyse_hbox)

    grid_pandda = QtGui.QGridLayout()
    grid_pandda.setColumnStretch(0, 20)
    grid_pandda.setRowStretch(0, 20)
    # left hand side: table with information about available datasets
    self.pandda_column_name = ['Sample ID',
                               'Refinement\nSpace Group',
                               'Resolution\n[Mn<I/sig(I)> = 1.5]',
                               'Dimple\nRcryst',
                               'Dimple\nRfree',
                               'Crystal Form\nName',
                               'PanDDA\nlaunched?',
                               'PanDDA\nhit?',
                               'PanDDA\nreject?',
                               'PanDDA\nStatus']

    self.pandda_analyse_data_table = QtGui.QTableWidget()
    self.pandda_analyse_data_table.setSortingEnabled(True)
    self.pandda_analyse_data_table.resizeColumnsToContents()
    self.pandda_analyse_data_table.setColumnCount(len(self.pandda_column_name))
    self.pandda_analyse_data_table.setHorizontalHeaderLabels(self.pandda_column_name)

    frame_pandda = QtGui.QFrame()
    grid_pandda.addWidget(self.pandda_analyse_data_table, 0, 0)


    # right hand side: input parameters for PANDDAs run
    frame_right = QtGui.QFrame()
    frame_right.setFrameShape(QtGui.QFrame.StyledPanel)

    self.pandda_analyse_input_params_vbox = QtGui.QVBoxLayout()

    pandda_input_dir_hbox = QtGui.QHBoxLayout()
    label = QtGui.QLabel('data directory')
    #        label.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Minimum)
    self.pandda_analyse_input_params_vbox.addWidget(label)
    self.pandda_input_data_dir_entry = QtGui.QLineEdit()
    self.pandda_input_data_dir_entry.setText(os.path.join(self.initial_model_directory, '*'))
    self.pandda_input_data_dir_entry.setFixedWidth(300)
    pandda_input_dir_hbox.addWidget(self.pandda_input_data_dir_entry)
    self.select_pandda_input_dir_button = QtGui.QPushButton("Select Input Template")
    #        self.select_pandda_input_dir_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
    self.select_pandda_input_dir_button.setMaximumWidth(200)
    self.select_pandda_input_dir_button.clicked.connect(panfunc.select_pandda_input_template)
    pandda_input_dir_hbox.addWidget(self.select_pandda_input_dir_button)
    self.pandda_analyse_input_params_vbox.addLayout(pandda_input_dir_hbox)

    pandda_pdb_style_hbox = QtGui.QHBoxLayout()
    label = QtGui.QLabel('pdb style')
    pandda_pdb_style_hbox.addWidget(label)
    self.pandda_pdb_style_entry = QtGui.QLineEdit()
    self.pandda_pdb_style_entry.setText('dimple.pdb')
    self.pandda_pdb_style_entry.setFixedWidth(200)
    pandda_pdb_style_hbox.addWidget(self.pandda_pdb_style_entry)
    self.pandda_analyse_input_params_vbox.addLayout(pandda_pdb_style_hbox)

    pandda_mtz_style_hbox = QtGui.QHBoxLayout()
    label = QtGui.QLabel('mtz style')
    pandda_mtz_style_hbox.addWidget(label)
    self.pandda_mtz_style_entry = QtGui.QLineEdit()
    self.pandda_mtz_style_entry.setText('dimple.mtz')
    self.pandda_mtz_style_entry.setFixedWidth(200)
    pandda_mtz_style_hbox.addWidget(self.pandda_mtz_style_entry)
    self.pandda_analyse_input_params_vbox.addLayout(pandda_mtz_style_hbox)

    pandda_output_dir_hbox = QtGui.QHBoxLayout()
    label = QtGui.QLabel('output directory')
    self.pandda_analyse_input_params_vbox.addWidget(label)
    self.pandda_output_data_dir_entry = QtGui.QLineEdit()
    self.pandda_output_data_dir_entry.setText(self.panddas_directory)
    self.pandda_output_data_dir_entry.setFixedWidth(300)
    pandda_output_dir_hbox.addWidget(self.pandda_output_data_dir_entry)
    self.select_pandda_output_dir_button = QtGui.QPushButton("Select PANNDAs Directory")
    self.select_pandda_output_dir_button.setMaximumWidth(200)
    self.select_pandda_output_dir_button.clicked.connect(settingsfunc.settings_button_clicked)
    pandda_output_dir_hbox.addWidget(self.select_pandda_output_dir_button)
    self.pandda_analyse_input_params_vbox.addLayout(pandda_output_dir_hbox)

    # qstat or local machine
    label = QtGui.QLabel('submit')
    #        label.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Minimum)
    self.pandda_analyse_input_params_vbox.addWidget(label)
    self.pandda_submission_mode_selection_combobox = QtGui.QComboBox()
    if self.external_software['qsub']:
        self.pandda_submission_mode_selection_combobox.addItem('qsub')
    self.pandda_submission_mode_selection_combobox.addItem('local machine')
    self.pandda_submission_mode_selection_combobox.setMaximumWidth(200)
    self.pandda_analyse_input_params_vbox.addWidget(self.pandda_submission_mode_selection_combobox)

    label = QtGui.QLabel('number of processors')
    #        label.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Minimum)
    self.pandda_analyse_input_params_vbox.addWidget(label)
    self.pandda_nproc = multiprocessing.cpu_count() - 1
    self.pandda_nproc_entry = QtGui.QLineEdit()
    self.pandda_nproc_entry.setText(str(self.pandda_nproc).replace(' ', ''))
    self.pandda_nproc_entry.setFixedWidth(200)
    self.pandda_analyse_input_params_vbox.addWidget(self.pandda_nproc_entry)

    label = QtGui.QLabel('order events by:')
    #        label.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Minimum)
    self.pandda_analyse_input_params_vbox.addWidget(label)
    self.pandda_sort_event_combobox = QtGui.QComboBox()
    self.pandda_sort_event_combobox.addItem('cluster_size')
    self.pandda_sort_event_combobox.addItem('z_peak')
    self.pandda_sort_event_combobox.setMaximumWidth(200)
    self.pandda_analyse_input_params_vbox.addWidget(self.pandda_sort_event_combobox)

    # crystal form option
    label = QtGui.QLabel('Use space group of reference file as filter:')
    #        label.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Minimum)
    self.pandda_analyse_input_params_vbox.addWidget(label)
    # reference file combobox, label with spg display
    hbox = QtGui.QHBoxLayout()
    self.reference_file_list = self.get_reference_file_list(' ')
    self.pandda_reference_file_selection_combobox = QtGui.QComboBox()
    self.populate_reference_combobox(self.pandda_reference_file_selection_combobox)
    self.pandda_reference_file_selection_combobox.activated[str].connect(panfunc.change_pandda_spg_label)
    hbox.addWidget(self.pandda_reference_file_selection_combobox)
    self.pandda_reference_file_spg_label = QtGui.QLabel()
    hbox.addWidget(self.pandda_reference_file_spg_label)
    self.pandda_analyse_input_params_vbox.addLayout(hbox)

    label = QtGui.QLabel('\n\n\nExpert Parameters (only change if you know what you are doing!):\n')
    #        label.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Minimum)
    self.pandda_analyse_input_params_vbox.addWidget(label)

    # minimum number of datasets
    label = QtGui.QLabel('min_build_datasets')
    #        label.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Minimum)
    self.pandda_analyse_input_params_vbox.addWidget(label)
    self.pandda_min_build_dataset_entry = QtGui.QLineEdit()
    self.pandda_min_build_dataset_entry.setText('40')
    self.pandda_min_build_dataset_entry.setFixedWidth(200)
    self.pandda_analyse_input_params_vbox.addWidget(self.pandda_min_build_dataset_entry)

    label = QtGui.QLabel('max_new_datasets')
    #        label.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Minimum)
    self.pandda_analyse_input_params_vbox.addWidget(label)
    self.pandda_max_new_datasets_entry = QtGui.QLineEdit()
    self.pandda_max_new_datasets_entry.setText('200')
    self.pandda_max_new_datasets_entry.setFixedWidth(200)
    self.pandda_analyse_input_params_vbox.addWidget(self.pandda_max_new_datasets_entry)

    label = QtGui.QLabel(
        'grid_spacing (default=0.6)\nNote: higher values speed up calculations, but maps might be less pretty)')
    #        label.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Minimum)
    self.pandda_analyse_input_params_vbox.addWidget(label)
    self.pandda_grid_spacing_entry = QtGui.QLineEdit()
    self.pandda_grid_spacing_entry.setText('0.6')
    self.pandda_grid_spacing_entry.setFixedWidth(200)
    self.pandda_analyse_input_params_vbox.addWidget(self.pandda_grid_spacing_entry)

    self.pandda_analyse_input_params_vbox.addStretch(1)

    frame_right.setSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
    frame_right.setLayout(self.pandda_analyse_input_params_vbox)

    grid_pandda.addWidget(frame_right, 0, 1, 5, 5)
    frame_pandda.setLayout(grid_pandda)
    self.pandda_analyse_hbox.addWidget(frame_pandda)

    # next three blocks display html documents created by pandda.analyse
    self.pandda_initial_html_file = os.path.join(self.panddas_directory, 'results_summaries', 'pandda_initial.html')
    self.pandda_analyse_html_file = os.path.join(self.panddas_directory, 'results_summaries', 'pandda_analyse.html')
    self.pandda_inspect_html_file = os.path.join(self.panddas_directory, 'results_summaries', 'pandda_inspect.html')

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

    self.panddas_results_vbox.addWidget(pandda_tab_widget)
    self.show_pandda_html_summary()

def refinement_tab(self):
    self.summary_vbox_for_table = QtGui.QVBoxLayout()
    self.summary_column_name = ['Sample ID',
                                'Compound ID',
                                'Refinement\nSpace Group',
                                'Refinement\nResolution',
                                'Refinement\nRcryst',
                                'Refinement\nRfree',
                                'Refinement\nOutcome',
                                'PanDDA site details',
                                'Refinement\nStatus']
    self.summary_table = QtGui.QTableWidget()
    self.summary_table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
    self.summary_table.setSortingEnabled(True)
    self.summary_table.setColumnCount(len(self.summary_column_name))
    self.summary_table.setHorizontalHeaderLabels(self.summary_column_name)
    self.summary_vbox_for_table.addWidget(self.summary_table)
    self.tab_dict[self.workflow_dict['Refinement']][1].addLayout(self.summary_vbox_for_table)

def deposition_tab(self):
    self.deposition_vbox = QtGui.QVBoxLayout()

    scroll = QtGui.QScrollArea()
    self.deposition_vbox.addWidget(scroll)
    scrollContent = QtGui.QWidget(scroll)

    scrollLayout = QtGui.QVBoxLayout(scrollContent)
    scrollContent.setLayout(scrollLayout)

    label_title = QtGui.QLabel('HTML export & ZENODO upload')
    label_title.setStyleSheet("font: 30pt Comic Sans MS")
    scrollLayout.addWidget(label_title)
    scrollLayout.addWidget(QtGui.QLabel(''))
    label_text = QtGui.QLabel(XChemToolTips.html_summary_introduction())
    label_text.setStyleSheet("font: 17pt Arial")
    scrollLayout.addWidget(label_text)
    scrollLayout.addWidget(QtGui.QLabel(''))
    image = QtGui.QLabel()
    pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'), 'image', 'html_summary_page.png'))

    image.setPixmap(pixmap)
    scrollLayout.addWidget(image)
    scrollLayout.addWidget(QtGui.QLabel(''))
    #        label_link=QtGui.QLabel('''<a href='ftp://ftp.sgc.ox.ac.uk/pub/tkrojer/XChemExplorer/Export_HTML_summary.pdf'>Click here for more information!</a>''')
    #        label_link.setOpenExternalLinks(True)
    #        scrollLayout.addWidget(label_link)
    #        scrollLayout.addWidget(QtGui.QLabel(''))
    label_heading = QtGui.QLabel('1. Specify HTML export directory in the settings tab')
    label_heading.setStyleSheet("font: bold 20pt Arial")
    scrollLayout.addWidget(label_heading)
    label_text = QtGui.QLabel(XChemToolTips.html_export_directory_background())
    label_text.setStyleSheet("font: 17pt Arial")
    scrollLayout.addWidget(label_text)
    if os.getcwd().startswith('/dls'):
        label_text = QtGui.QLabel('Note: default for labxchem project at DLS is <labxchem_directory>/processing/html.')
        label_text.setStyleSheet("font: 17pt Arial")
        scrollLayout.addWidget(label_text)
        label_text = QtGui.QLabel('In your case: ' + self.html_export_directory)
        label_text.setStyleSheet("font: 17pt Arial")
        scrollLayout.addWidget(label_text)
    scrollLayout.addWidget(QtGui.QLabel('\n'))
    label_heading = QtGui.QLabel("2. Prepare files for HTML export")
    label_heading.setStyleSheet("font: bold 20pt Arial")
    scrollLayout.addWidget(label_heading)
    label_text = QtGui.QLabel(XChemToolTips.html_export_step())
    label_text.setStyleSheet("font: 17pt Arial")
    scrollLayout.addWidget(label_text)
    button = QtGui.QPushButton('Export to HTML')
    button.clicked.connect(self.export_to_html)
    button.setMaximumWidth(200)
    scrollLayout.addWidget(button)
    scrollLayout.addWidget(QtGui.QLabel('\n'))
    label_heading = QtGui.QLabel("3. Prepare ICB files")
    label_heading.setStyleSheet("font: bold 20pt Arial")
    scrollLayout.addWidget(label_heading)

    label_text = QtGui.QLabel(XChemToolTips.icb_file_background())
    label_text.setStyleSheet("font: 17pt Arial")
    scrollLayout.addWidget(label_text)

    label_text = QtGui.QLabel(XChemToolTips.prepare_ICB_files())
    label_text.setStyleSheet("font: 17pt Arial")
    scrollLayout.addWidget(label_text)
    button = QtGui.QPushButton('Open ICM-pro')
    button.clicked.connect(self.open_icm)
    button.setMaximumWidth(200)
    scrollLayout.addWidget(button)
    scrollLayout.addWidget(QtGui.QLabel(''))
    image = QtGui.QLabel()
    pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'), 'image', 'drag_and_drop_icb_file.png'))
    image.setPixmap(pixmap)
    scrollLayout.addWidget(image)
    scrollLayout.addWidget(QtGui.QLabel(''))
    image = QtGui.QLabel()
    pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'), 'image', 'run_icm_script.png'))
    image.setPixmap(pixmap)
    scrollLayout.addWidget(image)
    scrollLayout.addWidget(QtGui.QLabel('\n'))
    label_heading = QtGui.QLabel("4. Prepare files for ZENODO upload")
    label_heading.setStyleSheet("font: bold 20pt Arial")
    scrollLayout.addWidget(label_heading)
    label_text = QtGui.QLabel(XChemToolTips.zenodo_upload_start(self.html_export_directory))
    label_text.setStyleSheet("font: 17pt Arial")
    scrollLayout.addWidget(label_text)
    button = QtGui.QPushButton('prepare files')
    button.clicked.connect(depofunc.prepare_files_for_zenodo_upload)
    button.setMaximumWidth(200)
    scrollLayout.addWidget(button)
    scrollLayout.addWidget(QtGui.QLabel('\n'))
    label_heading = QtGui.QLabel("5. ZENODO")
    label_heading.setStyleSheet("font: bold 20pt Arial")
    scrollLayout.addWidget(label_heading)
    label_text = QtGui.QLabel(XChemToolTips.zenodo_upload_part_one(self.html_export_directory))
    label_text.setStyleSheet("font: 17pt Arial")
    scrollLayout.addWidget(label_text)
    image = QtGui.QLabel()
    pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'), 'image', 'new_zenodo_upload.png'))
    image.setPixmap(pixmap)
    scrollLayout.addWidget(image)
    scrollLayout.addWidget(QtGui.QLabel('\n'))
    label_heading = QtGui.QLabel("5. ZENODO upload ID")
    label_heading.setStyleSheet("font: bold 20pt Arial")
    scrollLayout.addWidget(label_heading)
    label_text = QtGui.QLabel(XChemToolTips.zenodo_upload_part_two())
    label_text.setStyleSheet("font: 17pt Arial")
    scrollLayout.addWidget(label_text)
    image = QtGui.QLabel()
    pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'), 'image', 'zenodo_upload_id.png'))
    image.setPixmap(pixmap)
    scrollLayout.addWidget(image)
    label_text = QtGui.QLabel(XChemToolTips.zenodo_upload_part_three())
    label_text.setStyleSheet("font: 17pt Arial")
    scrollLayout.addWidget(label_text)
    hbox_zenodo_upload_id = QtGui.QHBoxLayout()
    hbox_zenodo_upload_id.addWidget(QtGui.QLabel('upload ID:'))
    self.zenodo_upload_id_entry = QtGui.QLineEdit()
    self.zenodo_upload_id_entry.setFixedWidth(200)
    hbox_zenodo_upload_id.addWidget(self.zenodo_upload_id_entry)
    hbox_zenodo_upload_id.addStretch(1)
    scrollLayout.addLayout(hbox_zenodo_upload_id)
    button = QtGui.QPushButton('update html files with upload ID')
    button.clicked.connect(self.update_html_for_zenodo_upload)
    button.setMaximumWidth(300)
    scrollLayout.addWidget(button)
    scrollLayout.addWidget(QtGui.QLabel('\n'))
    label_heading = QtGui.QLabel("6. ZENODO upload HTML files")
    label_heading.setStyleSheet("font: bold 20pt Arial")
    scrollLayout.addWidget(label_heading)
    label_text = QtGui.QLabel(XChemToolTips.zenodo_upload_part_four(self.html_export_directory))
    label_text.setStyleSheet("font: 17pt Arial")
    scrollLayout.addWidget(label_text)
    scrollLayout.addWidget(QtGui.QLabel('\n'))

    scrollLayout.addStretch(1)
    #        label_title.setStyleSheet("font: 30pt Comic Sans MS")
    #        label_heading.setStyleSheet("font: bold 20pt Arial")
    #        label_text.setStyleSheet("font: 17pt Arial")

    scroll.setWidget(scrollContent)

    self.tab_dict[self.workflow_dict['Deposition']][1].addLayout(self.deposition_vbox)

def settings_tab(self):
    self.settings_container = QtGui.QWidget()
    self.buttons_etc = QtGui.QWidget()
    self.settings_vbox = QtGui.QVBoxLayout()

    self.scroll = QtGui.QScrollArea(self.settings_container)
    self.settings_vbox.addWidget(self.scroll)
    scrollContent_settings = QtGui.QWidget(self.scroll)

    scrollLayout_settings = QtGui.QVBoxLayout(scrollContent_settings)
    scrollContent_settings.setLayout(scrollLayout_settings)

    # Settings Tab
    self.data_collection_vbox_for_settings = QtGui.QVBoxLayout()

    self.buttons_etc.setLayout(self.data_collection_vbox_for_settings)
    self.scroll.setWidget(self.buttons_etc)

    self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nProject Directory: - REQUIRED -'))
    settings_hbox_initial_model_directory = QtGui.QHBoxLayout()
    self.initial_model_directory_label = QtGui.QLabel(self.initial_model_directory)
    settings_hbox_initial_model_directory.addWidget(self.initial_model_directory_label)
    settings_buttoon_initial_model_directory = QtGui.QPushButton('Select Project Directory')
    settings_buttoon_initial_model_directory.clicked.connect(settingsfunc.settings_button_clicked)
    settings_hbox_initial_model_directory.addWidget(settings_buttoon_initial_model_directory)
    self.data_collection_vbox_for_settings.addLayout(settings_hbox_initial_model_directory)

    self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nReference Structure Directory: - OPTIONAL -'))
    settings_hbox_reference_directory = QtGui.QHBoxLayout()
    self.reference_directory_label = QtGui.QLabel(self.reference_directory)
    settings_hbox_reference_directory.addWidget(self.reference_directory_label)
    settings_buttoon_reference_directory = QtGui.QPushButton('Select Reference Structure Directory')
    settings_buttoon_reference_directory.clicked.connect(settingsfunc.settings_button_clicked)
    settings_hbox_reference_directory.addWidget(settings_buttoon_reference_directory)
    self.data_collection_vbox_for_settings.addLayout(settings_hbox_reference_directory)

    self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nData Source: - REQUIRED -'))
    settings_hbox_data_source_file = QtGui.QHBoxLayout()
    if self.data_source_file != '':
        self.data_source_file_label = QtGui.QLabel(os.path.join(self.database_directory, self.data_source_file))
    else:
        self.data_source_file_label = QtGui.QLabel('')
    settings_hbox_data_source_file.addWidget(self.data_source_file_label)
    settings_buttoon_data_source_file = QtGui.QPushButton('Select Data Source File')
    settings_buttoon_data_source_file.clicked.connect(settingsfunc.settings_button_clicked)
    settings_hbox_data_source_file.addWidget(settings_buttoon_data_source_file)
    self.data_collection_vbox_for_settings.addLayout(settings_hbox_data_source_file)

    # Data Collection
    self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nData Collection Directory: - OPTIONAL -'))

    settings_beamline_frame = QtGui.QFrame()
    settings_beamline_frame.setFrameShape(QtGui.QFrame.StyledPanel)
    settings_beamline_vbox = QtGui.QVBoxLayout()

    settings_hbox_beamline_directory = QtGui.QHBoxLayout()
    self.beamline_directory_label = QtGui.QLabel(self.beamline_directory)
    settings_hbox_beamline_directory.addWidget(self.beamline_directory_label)
    settings_buttoon_beamline_directory = QtGui.QPushButton('Select Data Collection Directory')
    settings_buttoon_beamline_directory.clicked.connect(settingsfunc.settings_button_clicked)
    settings_hbox_beamline_directory.addWidget(settings_buttoon_beamline_directory)
    settings_beamline_vbox.addLayout(settings_hbox_beamline_directory)

    settings_hbox_data_collection_summary_file = QtGui.QHBoxLayout()
    self.data_collection_summary_file_label = QtGui.QLabel(self.data_collection_summary_file)
    settings_hbox_data_collection_summary_file.addWidget(self.data_collection_summary_file_label)
    settings_button_data_collection_summary_file = QtGui.QPushButton('Select Existing\nCollection Summary File')
    settings_button_data_collection_summary_file.clicked.connect(settingsfunc.settings_button_clicked)
    settings_hbox_data_collection_summary_file.addWidget(settings_button_data_collection_summary_file)

    settings_button_new_data_collection_summary_file = QtGui.QPushButton('Assign New\nCollection Summary File')
    settings_button_new_data_collection_summary_file.clicked.connect(settingsfunc.settings_button_clicked)
    settings_hbox_data_collection_summary_file.addWidget(settings_button_new_data_collection_summary_file)

    settings_beamline_vbox.addLayout(settings_hbox_data_collection_summary_file)

    settings_beamline_frame.setLayout(settings_beamline_vbox)
    self.data_collection_vbox_for_settings.addWidget(settings_beamline_frame)

    self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nCCP4_SCR Directory: - OPTIONAL -'))
    settings_hbox_ccp4_scratch_directory = QtGui.QHBoxLayout()
    self.ccp4_scratch_directory_label = QtGui.QLabel(self.ccp4_scratch_directory)
    settings_hbox_ccp4_scratch_directory.addWidget(self.ccp4_scratch_directory_label)
    settings_buttoon_ccp4_scratch_directory = QtGui.QPushButton('Select CCP4_SCR Directory')
    settings_buttoon_ccp4_scratch_directory.clicked.connect(settingsfunc.settings_button_clicked)
    settings_hbox_ccp4_scratch_directory.addWidget(settings_buttoon_ccp4_scratch_directory)
    self.data_collection_vbox_for_settings.addLayout(settings_hbox_ccp4_scratch_directory)

    self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nPANDDAs directory: - OPTIONAL -'))
    settings_hbox_panddas_directory = QtGui.QHBoxLayout()
    self.panddas_directory_label = QtGui.QLabel(self.panddas_directory)
    settings_hbox_panddas_directory.addWidget(self.panddas_directory_label)
    settings_button_panddas_directory = QtGui.QPushButton('Select PANNDAs Directory')
    settings_button_panddas_directory.clicked.connect(settingsfunc.settings_button_clicked)
    settings_hbox_panddas_directory.addWidget(settings_button_panddas_directory)
    self.data_collection_vbox_for_settings.addLayout(settings_hbox_panddas_directory)

    self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nHTML export directory: - OPTIONAL -'))
    settings_hbox_html_export_directory = QtGui.QHBoxLayout()
    self.html_export_directory_label = QtGui.QLabel(self.html_export_directory)
    settings_hbox_html_export_directory.addWidget(self.html_export_directory_label)
    settings_button_html_export_directory = QtGui.QPushButton('Select HTML Export Directory')
    settings_button_html_export_directory.clicked.connect(settingsfunc.settings_button_clicked)
    settings_hbox_html_export_directory.addWidget(settings_button_html_export_directory)
    self.data_collection_vbox_for_settings.addLayout(settings_hbox_html_export_directory)

    self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nGroup deposition directory: - OPTIONAL -'))
    settings_hbox_group_deposition_directory = QtGui.QHBoxLayout()
    self.group_deposition_directory_label = QtGui.QLabel(self.group_deposit_directory)
    settings_hbox_group_deposition_directory.addWidget(self.group_deposition_directory_label)
    settings_button_group_deposition_directory = QtGui.QPushButton('Select Group deposition Directory')
    settings_button_group_deposition_directory.clicked.connect(settingsfunc.settings_button_clicked)
    settings_hbox_group_deposition_directory.addWidget(settings_button_group_deposition_directory)
    self.data_collection_vbox_for_settings.addLayout(settings_hbox_group_deposition_directory)

    self.data_collection_vbox_for_settings.setSpacing(0)
    self.data_collection_vbox_for_settings.setContentsMargins(30, 0, 0, 0)

    self.buttons_etc.resize(self.screen.width() - 100, self.buttons_etc.sizeHint().height())
    self.tab_dict[self.workflow_dict['Settings']][1].addLayout(self.settings_vbox)

def tabs_setup(self):
    # setup main tabs
    self.main_tab_widget = QtGui.QTabWidget()
    self.main_tab_widget.setMaximumSize(self.screen.width(),self.screen.height()-245)
    self.tab_dict = add_subtabs(self, self.workflow, self.main_tab_widget)
    overview_tab(self)
    dataset_tab(self)
    maps_tab(self)
    panddas_tab(self)
    refinement_tab(self)
    deposition_tab(self)
    settings_tab(self)


















