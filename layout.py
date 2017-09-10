import sys, os, subprocess
from PyQt4 import QtGui, QtCore

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'lib'))

import XChemToolTips


class LayoutObjects():
    def __init__(self):
        pass

    # function to initialise the top menu bar
    def initialise_menu_bar(self):
        ################################################################################################################
        #                                                                                                              #
        #                                              MENU BAR - TOP OF GUI                                           #
        #                                                                                                              #
        ################################################################################################################

        # initiate menu widget
        menu_bar = QtGui.QMenuBar()

        # dictionary containing config for top menu setup
        # menu dict = { 'order letter: menu_item_name': ["name_in_menu",
        #                                       [
        #                                           ['text_in_menu', 'shortcut', 'trigger function']],
        #                                           [...],
        #                                       ]]
        #             }
        menu_dict = {'A: file': ["&File",
                                 [
                                     ['Open Config File', 'Ctrl+O', 'object.open_config_file'],
                                     ['Save Config File', 'Ctrl+S', 'object.save_config_file'],
                                     ['Quit', 'Ctrl+Q', 'object.quit_xce']
                                 ]],
                     'B: datasource': ["&Datasource",
                                       [
                                           ['Reload Samples From Datasource', '',
                                            'object.datasource_menu_reload_samples'],
                                           ['Save Samples to Datasource', '', 'object.datasource_menu_save_samples'],
                                           ['Import CSV file into Datasource', '',
                                            'object.datasource_menu_import_csv_file'],
                                           ['Export CSV file from Datasource', '',
                                            'object.datasource_menu_export_csv_file'],
                                           ['Select columns to show', '', 'object.select_datasource_columns_to_display'],
                                           ['Create New Datasource (SQLite)', '', 'object.create_new_data_source'],
                                           ['Export CSV for wonka', '', 'object.export_data_for_WONKA']
                                       ]],
                     'C: preferences': ["&Preferences",
                                        [
                                            ['Edit preferences', '', 'object.show_preferences']
                                        ]],
                     'D: deposition': ["&Deposition",
                                       [
                                           ['Edit information', '', 'object.deposition_data'],
                                           ['Export to HTML', '', 'object.export_to_html'],
                                           ['Find PanDDA apo structures', '',
                                            'object.create_missing_apo_records_in_depositTable'],
                                           ['Update file info of apo structures', '',
                                            'object.update_file_information_of_apo_records'],
                                           ['Prepare mmcif for apo structures', '',
                                            'object.prepare_models_for_deposition'],
                                           # ^ this is added to a dictionary somewhere, so need to check what it interferes
                                           #  with when code changed
                                           ['Prepare mmcif for ligand bound structures', '',
                                            'object.prepare_models_for_deposition'],
                                           # ^ this is added to a dictionary somewhere, so need to check what it interferes
                                           #  with when code changed
                                           ['Copy files to group deposition directory', '',
                                            'object.prepare_for_group_deposition_upload'],
                                           ['Update DB with PDB codes', '', 'object.enter_pdb_codes'],
                                           ['Check SMILES', '', 'object.check_smiles_in_db_and_pdb']
                                       ]],
                     'E: help': ["&Help",
                                 [
                                     ['Open XCE tutorial', '', str(
                                         'lambda: object.layout_funcs.openFile("/dls/science/groups/i04-1/software/docs/'
                                         'XChemExplorer.pdf")')],
                                     ['Troubleshooting', '', str(
                                         'lambda: object.layout_funcs.openFile("/dls/science/groups/i04-1/software/'
                                         'xce_troubleshooting.pdf")')]
                                 ]]

                     }

        # create menu from menu dictionary
        menu_bar = self.layout_funcs.setup_menubar(self, menu_bar, menu_dict)

        # END OF MENU BAR - CODE BELOW: stuff removed from apo structure stuff that appears might have a funky
        # consequence - work out later.

        # self.prepare_mmcif_files_dict = {}
        # self.prepare_mmcif_files_dict['apo'] = prepare_mmcif_files_for_apo_structures
        # self.prepare_mmcif_files_dict['ligand_bound'] = prepare_mmcif_files_for_ligand_bound_structures

        return menu_bar


    # function containing setup for bottom boxes
    def initialise_bottom_boxes(self):
        ################################################################################################################
        #                                                                                                              #
        #                                             DATASOURCE BUTTON                                                #
        #                                                                                                              #
        ################################################################################################################
        # config settings fot the 'update datasource' button
        datasource_button_dict = {'datasource_button': [r"Update Tables\nFrom Datasource",
                                                        [
                                                            ['XChemToolTips.update_from_datasource_button_tip()',
                                                             # tooltip
                                                             'QPushButton { padding: 1px; margin: 1px; '
                                                             'background: rgb(140,140,140) }',
                                                             # stylesheet
                                                             'object.headlineLabelfont',  # font
                                                             'object.datasource_menu_reload_samples']  # action
                                                        ]]}
        # setup datasource button
        update_from_datasource_button = self.layout_funcs.setup_push_button(self, datasource_button_dict)

        ################################################################################################################
        #                                                                                                              #
        #                                               DATASETS BOX                                                   #
        #                                                                                                              #
        ################################################################################################################
        # settings for the run button
        dataset_task_run_button_dict = {'dataset_run_button':
                                            [r"Run",
                                             [
                                                 ['XChemToolTips.dataset_task_run_button_tip()',  # tooltip
                                                  'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                  '',  # font
                                                  'object.button_clicked']  # action
                                             ]]}
        # setup the run button with push button function
        self.dataset_task_run_button = self.layout_funcs.setup_push_button(self, dataset_task_run_button_dict)

        # settings for the status button
        dataset_task_status_button_dict = {'dataset_status_button':
                                               [r"Status",
                                                [
                                                    ['XChemToolTips.dataset_task_status_button_tip()',  # tooltip
                                                     'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                     '',  # font
                                                     'object.button_clicked']  # action
                                                ]]}
        # setup the task button with push button function
        self.dataset_task_status_button = self.layout_funcs.setup_push_button(self, dataset_task_status_button_dict)

        # array of both button objects to apply to bottom box layout
        dataset_buttons = [self.dataset_task_run_button, self.dataset_task_status_button]

        # items that will appear it the dropdown (combobox)
        dataset_tasks = ['Get New Results from Autoprocessing',
                         'Run DIMPLE on All Autoprocessing MTZ files',
                         'Rescore Datasets',
                         'Read PKL file',
                         'Run xia2 on selected datasets',
                         'Run xia2 on selected datasets - overwrite']

        # label for the bottom box layout
        dataset_label = "Datasets"

        # return the frame and combobox from the bottom box setup function
        frame_dataset_task, self.dataset_tasks_combobox = self.layout_funcs.bottom_box_setup(self, dataset_label,
                                                                                             dataset_tasks,
                                                                                             'XChemToolTips.dataset_task_tip()',
                                                                                             dataset_buttons,
                                                                                             'background: rgb(240, 255, 140); ')

        # define the combobox and buttons in dictionary key to determine behaviour
        self.workflow_widget_dict['Datasets'] = [self.dataset_tasks_combobox, self.dataset_task_run_button,
                                                 self.dataset_task_status_button]

        ################################################################################################################
        #                                                                                                              #
        #                                            MAPS & RESTRAINTS BOX                                             #
        #                                                                                                              #
        ################################################################################################################
        # settings for the run button
        map_cif_file_task_run_button_dict = {'dataset_run_button':
                                                 [r"Run",
                                                  [
                                                      ['XChemToolTips.map_cif_file_task_run_button_tip()',  # tooltip
                                                       'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                       '',  # font
                                                       'object.button_clicked']  # action
                                                  ]]}
        # setup the run button with push button function
        self.map_cif_file_task_run_button = self.layout_funcs.setup_push_button(self, map_cif_file_task_run_button_dict)

        # settings for the status button
        map_cif_file_task_status_button_dict = {'dataset_status_button':
                                                    [r"Status",
                                                     [
                                                         ['XChemToolTips.map_cif_file_task_status_button_tip()',
                                                          # tooltip
                                                          'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                          '',  # font
                                                          'object.button_clicked']  # action
                                                     ]]}
        # setup the task button with push button function
        self.map_cif_file_task_status_button = self.layout_funcs.setup_push_button(self,
                                                                                   map_cif_file_task_status_button_dict)

        # array of both button objects to apply to bottom box layout
        map_cif_file_buttons = [self.map_cif_file_task_run_button, self.map_cif_file_task_status_button]

        # items that will appear it the dropdown (combobox)
        map_cif_file_tasks = ['Run DIMPLE on selected MTZ files',
                              'Remove selected DIMPLE PDB/MTZ files',
                              'Create CIF/PDB/PNG file of ALL compounds',
                              'Create CIF/PDB/PNG file of NEW compounds',
                              'Create CIF/PDB/PNG file of SELECTED compounds']

        # label for the bottom box layout
        map_cif_file_label = "Maps & Restraints"

        # return the frame and combobox from the bottom box setup function
        frame_map_cif_file_task, self.map_cif_file_tasks_combobox = self.layout_funcs.bottom_box_setup(self,
                                                                                                       map_cif_file_label,
                                                                                                       map_cif_file_tasks,
                                                                                                       'XChemToolTips.map_cif_file_'
                                                                                                       'task_tip()',
                                                                                                       map_cif_file_buttons,
                                                                                                       'background: rgb(140, 255, '
                                                                                                       '150); ')

        # define the combobox and buttons in dictionary key to determine behaviour
        self.workflow_widget_dict['Maps'] = [self.map_cif_file_tasks_combobox, self.map_cif_file_task_run_button,
                                             self.map_cif_file_task_status_button]

        ################################################################################################################
        #                                                                                                              #
        #                                               HIT IDENTIFICATION BOX                                         #
        #                                                                                                              #
        ################################################################################################################
        # settings for the run button
        panddas_file_task_run_button_dict = {'dataset_run_button':
                                                 [r"Run",
                                                  [
                                                      ['XChemToolTips.panddas_file_task_run_button_tip()',  # tooltip
                                                       'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                       '',  # font
                                                       'object.button_clicked']  # action
                                                  ]]}
        # setup the run button with push button function
        self.panddas_file_task_run_button = self.layout_funcs.setup_push_button(self, panddas_file_task_run_button_dict)

        # settings for the status button
        panddas_file_task_status_button_dict = {'dataset_status_button':
                                                    [r"Status",
                                                     [
                                                         ['XChemToolTips.panddas_file_task_status_button_tip()',
                                                          # tooltip
                                                          'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                          '',  # font
                                                          'object.button_clicked']  # action
                                                     ]]}
        # setup the task button with push button function
        self.panddas_file_task_status_button = self.layout_funcs.setup_push_button(self,
                                                                                   panddas_file_task_status_button_dict)

        # array of both button objects to apply to bottom box layout
        panddas_file_buttons = [self.panddas_file_task_run_button, self.panddas_file_task_status_button]

        # items that will appear it the dropdown (combobox)
        panddas_file_tasks = ['pandda.analyse',
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
                              'Build ground state model']

        # label for the bottom box layout
        panddas_file_label = "Hit Identification"

        # return the frame and combobox from the bottom box setup function
        frame_panddas_file_task, self.panddas_file_tasks_combobox = self.layout_funcs.bottom_box_setup(self,
                                                                                                       panddas_file_label,
                                                                                                       panddas_file_tasks,
                                                                                                       'XChemToolTips.panddas_file_'
                                                                                                       'task_tip()',
                                                                                                       panddas_file_buttons,
                                                                                                       'background: rgb(140,200,255)'
                                                                                                       '; ')

        # define the combobox and buttons in dictionary key to determine behaviour
        self.workflow_widget_dict['PANDDAs'] = [self.panddas_file_tasks_combobox, self.panddas_file_task_run_button,
                                                self.panddas_file_task_status_button]

        ################################################################################################################
        #                                                                                                              #
        #                                                      REFINEMENT BOX                                          #
        #                                                                                                              #
        ################################################################################################################
        # settings for the run button
        refine_file_task_run_button_dict = {'dataset_run_button':
                                                [r"Run",
                                                 [
                                                     ['XChemToolTips.refine_file_task_run_button_tip()',  # tooltip
                                                      'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                      '',  # font
                                                      'object.button_clicked']  # action
                                                 ]]}
        # setup the run button with push button function
        self.refine_file_task_run_button = self.layout_funcs.setup_push_button(self, refine_file_task_run_button_dict)

        # settings for the status button
        refine_file_task_status_button_dict = {'dataset_status_button':
                                                   [r"Status",
                                                    [
                                                        ['XChemToolTips.refine_file_task_status_button_tip()',
                                                         # tooltip
                                                         'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                         '',  # font
                                                         'object.button_clicked']  # action
                                                    ]]}
        # setup the task button with push button function
        self.refine_file_task_status_button = self.layout_funcs.setup_push_button(self, refine_file_task_status_button_dict)

        # array of both button objects to apply to bottom box layout
        refine_file_buttons = [self.refine_file_task_run_button, self.refine_file_task_status_button]

        # items that will appear it the dropdown (combobox)
        refine_file_tasks = ['Open COOT',
                             'Open COOT - new interface',
                             'Open COOT for old PanDDA',
                             'Update Deposition Table',
                             'Prepare Group Deposition']

        # label for the bottom box layout
        refine_file_label = "Refinement"

        # return the frame and combobox from the bottom box setup function
        frame_refine_file_task, self.refine_file_tasks_combobox = self.layout_funcs.bottom_box_setup(self,
                                                                                                     refine_file_label,
                                                                                                     refine_file_tasks,
                                                                                                     'XChemToolTips.refine_file_task'
                                                                                                     '_tip()',
                                                                                                     refine_file_buttons,
                                                                                                     'background: rgb(245, 190, 255)'
                                                                                                     ';')

        # define the combobox and buttons in dictionary key to determine behaviour
        self.workflow_widget_dict['Refinement'] = [self.refine_file_tasks_combobox, self.refine_file_task_run_button,
                                                   self.refine_file_task_status_button]

        return update_from_datasource_button, frame_dataset_task, frame_map_cif_file_task, frame_panddas_file_task, \
               frame_refine_file_task

    def main_layout(self):
        # initialise menu bar
        menu_bar = self.initialise_menu_bar()

        # initialise bottom boxes
        update_from_datasource_button, frame_dataset_task, frame_map_cif_file_task, frame_panddas_file_task, \
        frame_refine_file_task = self.initialise_bottom_boxes()

        # set all table columns
        #self.define_all_tables()

        # Tab layout & content
        # --------------------
        #
        # Overview
        # |- datasource - TABLE
        # |- summary - GRAPH
        #
        # Datasets
        # |- summary - TABLE
        # |- reprocess - TABLE
        #
        # Maps - TABLE
        #
        # PANDDAS
        # |- pandda.analyse - TABLE
        # |- Dataset Summary ------------------
        # |- Processing Output                  |   HTML
        # |- pandda.inspect                     |
        # |- Statistical Map Summaries --------
        #
        # Refinement - TABLE
        #
        # Deposition
        #
        # Settings
        #
        # NOTES:-
        # -------
        # 1. all table columns are defined by self.define_all_tables()
        # 2. all table column definitions are named by <tab>_<subtab>_table_columns (subtab omitted if not relevant)
        # 3. subtabs are created by self.layout_funcs.make_tab_dict()

        ################################################################################################################
        #                                                                                                              #
        #                                                 OVERVIEW TAB                                                 #
        #                                                                                                              #
        ################################################################################################################
        # define subtab list, widget and dict
        overview_tab_list = ['Data Source', 'Summary']
        self.overview_tab_widget = QtGui.QTabWidget()
        self.overview_tab_dict = {}

        # make subtabs
        self.layout_funcs.make_tab_dict(overview_tab_list, self.overview_tab_widget, self.overview_tab_dict)

        # add overview subtabs to overview tab
        self.tab_dict[self.workflow_dict['Overview']][1].addWidget(self.overview_tab_widget)

        # initiate the table in overview/datasource
        self.overview_datasource_table = QtGui.QTableWidget()
        self.overview_datasource_table.setSortingEnabled(True)
        self.overview_datasource_table.resizeColumnsToContents()

        # add table to datasource subtab
        self.overview_tab_dict['Data Source'][1].addWidget(self.overview_datasource_table)

        # initiate the graph in overview/summary
        self.overview_figure, self.overview_axes = plt.subplots()
        self.overview_canvas = FigureCanvas(self.overview_figure)
        self.update_summary_plot()

        # add graph to summary tab
        self.overview_tab_dict['Summary'][1].addWidget(self.overview_canvas)

        ################################################################################################################
        #                                                                                                              #
        #                                                 DATASETS TAB                                                 #
        #                                                                                                              #
        ################################################################################################################
        # define subtab list, widget and dict
        datasets_tab_list = ['Summary', 'Reprocess']
        self.datasets_tab_widget = QtGui.QTabWidget()
        self.datasets_tab_dict = {}

        # make subtabs
        self.layout_funcs.make_tab_dict(datasets_tab_list, self.datasets_tab_widget, self.datasets_tab_dict)

        # main body - things that are always displayed

        # add a container to hold everythting and add to main tab layout
        self.datasets_data_collection_vbox = QtGui.QVBoxLayout()
        self.tab_dict[self.workflow_dict['Datasets']][1].addLayout(self.datasets_data_collection_vbox)

        # add a horizontal box to hold option to autocheck for new data
        self.autocheck_hbox = QtGui.QHBoxLayout()

        # checkbox for autocollect
        self.check_for_new_data_collection = QtGui.QCheckBox('Check for new data collection every two minutes')
        self.layout_funcs.add_checkbox(self, self.check_for_new_data_collection, 'object.continously_check_for_new_data_collection')

        # select target dropdown
        select_target_label = QtGui.QLabel('Select Target: ')
        select_target_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.target_selection_combobox = QtGui.QComboBox()
        self.populate_target_selection_combobox(self.target_selection_combobox)
        self.target_selection_combobox.activated[str].connect(self.target_selection_combobox_activated)
        self.target = str(self.target_selection_combobox.currentText())

        self.autocheck_hbox_widgets = [self.check_for_new_data_collection, select_target_label,
                                       self.target_selection_combobox]

        self.layout_funcs.add_to_box(self.autocheck_hbox, self.autocheck_hbox_widgets)

        # add target dropdown to top bar
        self.datasets_data_collection_vbox.addLayout(self.autocheck_hbox)

        # summary sub-tab

        # table
        self.datasets_summary_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(self.datasets_summary_table, self.datasets_summary_table_columns)  # setup table
        self.datasets_summarys_vbox_for_table = QtGui.QVBoxLayout()  # setup layout to hold table
        self.datasets_summarys_vbox_for_table.addWidget(self.datasets_summary_table)  # add table to layout
        self.datasets_tab_dict['Summary'][1].addLayout(self.datasets_summarys_vbox_for_table)  # add layout to tab
        self.datasets_summarys_vbox_for_details = QtGui.QVBoxLayout()  # vbox for details
        self.data_collection_details_currently_on_display = None  # switch for displaying/updating table
        self.datasets_tab_dict['Summary'][1].addLayout(self.datasets_summarys_vbox_for_details)  # add details

        self.datasets_data_collection_vbox.addWidget(self.datasets_tab_widget)  # add subtabs to main tab

        # reprocessing sub-tab
        # top options
        reprocess_vbox = QtGui.QVBoxLayout()  # box to hold reprocessing subtab content

        frame_select = QtGui.QFrame()  # frame to hold top options box
        self.hbox_select = QtGui.QHBoxLayout()  # top options box

        dc_label = QtGui.QLabel('Data collection directory: ')
        dc_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.diffraction_data_dir_label = QtGui.QLabel(self.diffraction_data_directory)
        self.diffraction_data_dir_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)

        select_button = QtGui.QPushButton("Select")
        select_button.clicked.connect(self.select_diffraction_data_directory)

        search_button = QtGui.QPushButton("Search Datasets")
        search_button.clicked.connect(self.search_for_datasets)

        self.diffraction_data_search_label = QtGui.QLabel(self.diffraction_data_search_info)

        translate_label = QtGui.QLabel('translate: datasetID -> sampleID')
        translate_label.setAlignment(QtCore.Qt.AlignCenter)
        csv_button = QtGui.QPushButton('Open CSV')
        csv_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
        csv_button.clicked.connect(self.translate_datasetID_to_sampleID)

        self.hbox_select_widgets = [dc_label, select_button, search_button, self.diffraction_data_search_label,
                                    translate_label, csv_button]

        self.layout_funcs.add_to_box(self.hbox_select, self.hbox_select_widgets)

        frame_select.setLayout(self.hbox_select)

        # table
        self.datasets_reprocess_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(self.datasets_reprocess_table, self.datasets_reprocess_columns)

        # create context menu
        self.popMenu_for_datasets_reprocess_table = QtGui.QMenu()
        run_xia2_on_selected = QtGui.QAction("mark selected for reprocessing", self.window)
        run_xia2_on_selected.triggered.connect(self.select_sample_for_xia2)
        self.popMenu_for_datasets_reprocess_table.addAction(run_xia2_on_selected)
        self.datasets_reprocess_table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.datasets_reprocess_table.customContextMenuRequested.connect(self.on_context_menu_reprocess_data)

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)
        data_protocol_hbox = QtGui.QHBoxLayout()

        frame_options = QtGui.QFrame()
        frame_options.setFrameShape(QtGui.QFrame.StyledPanel)
        hbox_options = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Data processing protocol: ')
        hbox_options.addWidget(label)
        self.xia2_3d_checkbox = QtGui.QCheckBox('xia2 3d')
        hbox_options.addWidget(self.xia2_3d_checkbox)
        self.xia2_3dii_checkbox = QtGui.QCheckBox('xia2 3dii')
        hbox_options.addWidget(self.xia2_3dii_checkbox)
        self.xia2_dials_checkbox = QtGui.QCheckBox('Dials')
        hbox_options.addWidget(self.xia2_dials_checkbox)
        frame_options.setLayout(hbox_options)

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
        button.clicked.connect(self.select_reprocess_reference_mtz)
        hbox_ref.addWidget(button)
        frame_ref.setLayout(hbox_ref)

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

        data_ptotocol_hbox_widgets = [frame_options, frame_sg, frame_ref, frame_isigma, frame_cc_half]
        self.layout_funcs.add_to_box(data_protocol_hbox, data_ptotocol_hbox_widgets)

        frame.setLayout(data_protocol_hbox)

        self.reprocess_hbox_widgets = [frame_select, self.datasets_reprocess_table, frame]
        self.layout_funcs.add_to_box(reprocess_vbox, self.reprocess_hbox_widgets)

        self.datasets_tab_dict['Reprocess'][1].addLayout(reprocess_vbox)

        ################################################################################################################
        #                                                                                                              #
        #                                                    MAPS TAB                                                  #
        #                                                                                                              #
        ################################################################################################################
        initial_model_checkbutton_hbox = QtGui.QHBoxLayout()

        self.select_sample_for_dimple_box = QtGui.QCheckBox('(de-)select all samples for DIMPLE')
        self.layout_funcs.add_checkbox(self, self.select_sample_for_dimple_box, 'object.set_run_dimple_flag')
        initial_model_checkbutton_hbox.addWidget(self.select_sample_for_dimple_box)

        set_new_reference_button = QtGui.QPushButton("Set New Reference (if applicable)")
        set_new_reference_button.clicked.connect(self.set_new_reference_if_applicable)
        initial_model_checkbutton_hbox.addWidget(set_new_reference_button)

        refresh_reference_file_list_button = QtGui.QPushButton("Refresh reference file list")
        refresh_reference_file_list_button.clicked.connect(self.refresh_reference_file_list)
        initial_model_checkbutton_hbox.addWidget(refresh_reference_file_list_button)

        self.reference_file_list = self.get_reference_file_list(' ')
        self.reference_file_selection_combobox = QtGui.QComboBox()
        self.populate_reference_combobox(self.reference_file_selection_combobox)
        initial_model_checkbutton_hbox.addWidget(self.reference_file_selection_combobox)

        self.tab_dict[self.workflow_dict['Maps']][1].addLayout(initial_model_checkbutton_hbox)
        self.initial_model_vbox_for_table = QtGui.QVBoxLayout()

        # table
        self.maps_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(self.maps_table, self.maps_table_columns)
        self.initial_model_vbox_for_table.addWidget(self.maps_table)
        self.tab_dict[self.workflow_dict['Maps']][1].addLayout(self.initial_model_vbox_for_table)

        # create context menu
        self.popMenu_for_maps_table = QtGui.QMenu()
        run_dimple = QtGui.QAction("mark selected for dimple run", self.window)
        run_dimple.triggered.connect(self.select_sample_for_dimple)
        self.popMenu_for_maps_table.addAction(run_dimple)
        self.maps_table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.maps_table.customContextMenuRequested.connect(self.on_context_menu_initial_model)

        ################################################################################################################
        #                                                                                                              #
        #                                                   PANDDA TAB                                                 #
        #                                                                                                              #
        ################################################################################################################
        pandda_tab_list = ['pandda.analyse',
                           'Dataset Summary',
                           'Processing Output',
                           'pandda.inspect',
                           'Statistical Map Summaries']

        self.pandda_tab_widget = QtGui.QTabWidget()
        self.pandda_tab_dict = {}
        self.layout_funcs.make_tab_dict(pandda_tab_list, self.pandda_tab_widget, self.pandda_tab_dict)

        ## PanDDA tab
        self.panddas_results_vbox = QtGui.QVBoxLayout()
        self.tab_dict[self.workflow_dict['PANDDAs']][1].addLayout(self.panddas_results_vbox)

        self.pandda_analyse_hbox = QtGui.QHBoxLayout()
        self.pandda_tab_dict['pandda.analyse'][1].addLayout(self.pandda_analyse_hbox)
        self.pandda_map_layout = QtGui.QVBoxLayout()
        self.pandda_map_list = QtGui.QComboBox()
        self.pandda_maps_html = QtWebKit.QWebView()
        self.pandda_map_layout.addWidget(self.pandda_map_list)
        self.pandda_map_layout.addWidget(self.pandda_maps_html)

        self.pandda_tab_dict['Statistical Map Summaries'][1].addLayout(self.pandda_map_layout)
        self.pandda_maps_html.show()

        grid_pandda = QtGui.QGridLayout()
        grid_pandda.setColumnStretch(0, 20)
        grid_pandda.setRowStretch(0, 20)
        # left hand side: table with information about available datasets

        # table
        self.pandda_analyse_data_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(self.pandda_analyse_data_table, self.pandda_table_columns)

        frame_pandda = QtGui.QFrame()
        grid_pandda.addWidget(self.pandda_analyse_data_table, 0, 0)

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
        self.pandda_status_label.setFont(QtGui.QFont("Arial", 25, QtGui.QFont.Bold))
        grid_pandda.addWidget(self.pandda_status_label, 3, 0)

        # right hand side: input parameters for PANDDAs run
        frame_right = QtGui.QFrame()
        frame_right.setFrameShape(QtGui.QFrame.StyledPanel)

        self.pandda_analyse_input_params_vbox = QtGui.QVBoxLayout()

        pandda_input_dir_hbox = QtGui.QHBoxLayout()
        label = QtGui.QLabel('data directory')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        self.pandda_input_data_dir_entry = QtGui.QLineEdit()
        self.pandda_input_data_dir_entry.setText(os.path.join(self.initial_model_directory, '*'))
        self.pandda_input_data_dir_entry.setFixedWidth(300)
        pandda_input_dir_hbox.addWidget(self.pandda_input_data_dir_entry)
        self.select_pandda_input_dir_button = QtGui.QPushButton("Select Input Template")
        self.select_pandda_input_dir_button.setMaximumWidth(200)
        self.select_pandda_input_dir_button.clicked.connect(self.select_pandda_input_template)
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
        self.select_pandda_output_dir_button.clicked.connect(self.settings_button_clicked)
        pandda_output_dir_hbox.addWidget(self.select_pandda_output_dir_button)
        self.pandda_analyse_input_params_vbox.addLayout(pandda_output_dir_hbox)

        # qstat or local machine
        label = QtGui.QLabel('submit')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        self.pandda_submission_mode_selection_combobox = QtGui.QComboBox()
        if self.external_software['qsub']:
            self.pandda_submission_mode_selection_combobox.addItem('qsub')
        self.pandda_submission_mode_selection_combobox.addItem('local machine')
        self.pandda_submission_mode_selection_combobox.setMaximumWidth(200)
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_submission_mode_selection_combobox)

        label = QtGui.QLabel('number of processors')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        self.pandda_nproc = multiprocessing.cpu_count() - 1
        self.pandda_nproc_entry = QtGui.QLineEdit()
        self.pandda_nproc_entry.setText(str(self.pandda_nproc).replace(' ', ''))
        self.pandda_nproc_entry.setFixedWidth(200)
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_nproc_entry)

        label = QtGui.QLabel('order events by:')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        self.pandda_sort_event_combobox = QtGui.QComboBox()
        self.pandda_sort_event_combobox.addItem('cluster_size')
        self.pandda_sort_event_combobox.addItem('z_peak')
        self.pandda_sort_event_combobox.setMaximumWidth(200)
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_sort_event_combobox)

        # crystal form option
        label = QtGui.QLabel('Use space group of reference file as filter:')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        # reference file combobox, label with spg display
        hbox = QtGui.QHBoxLayout()
        self.reference_file_list = self.get_reference_file_list(' ')
        self.pandda_reference_file_selection_combobox = QtGui.QComboBox()
        self.populate_reference_combobox(self.pandda_reference_file_selection_combobox)
        self.pandda_reference_file_selection_combobox.activated[str].connect(self.change_pandda_spg_label)
        hbox.addWidget(self.pandda_reference_file_selection_combobox)
        self.pandda_reference_file_spg_label = QtGui.QLabel()
        hbox.addWidget(self.pandda_reference_file_spg_label)
        self.pandda_analyse_input_params_vbox.addLayout(hbox)

        label = QtGui.QLabel('\nExpert Parameters (only change if you know what you are doing!):')
        self.pandda_analyse_input_params_vbox.addWidget(label)

        self.wilson_checkbox = QtGui.QCheckBox('Wilson B-factor Scaling')
        self.layout_funcs.add_checkbox(self, self.wilson_checkbox, 'object.set_run_dimple_flag')

        # minimum number of datasets
        label = QtGui.QLabel('min_build_datasets')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        self.pandda_min_build_dataset_entry = QtGui.QLineEdit()
        self.pandda_min_build_dataset_entry.setText('40')
        self.pandda_min_build_dataset_entry.setFixedWidth(200)
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_min_build_dataset_entry)

        label = QtGui.QLabel('max_new_datasets')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        self.pandda_max_new_datasets_entry = QtGui.QLineEdit()
        self.pandda_max_new_datasets_entry.setText('500')
        self.pandda_max_new_datasets_entry.setFixedWidth(200)
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_max_new_datasets_entry)

        label = QtGui.QLabel(
            'grid_spacing (default=0.5)\nNote: higher values speed up calculations, but maps might be less pretty)')
        self.pandda_analyse_input_params_vbox.addWidget(label)
        self.pandda_grid_spacing_entry = QtGui.QLineEdit()
        self.pandda_grid_spacing_entry.setText('0.5')
        self.pandda_grid_spacing_entry.setFixedWidth(200)
        self.pandda_analyse_input_params_vbox.addWidget(self.pandda_grid_spacing_entry)

        self.pandda_analyse_input_params_vbox.addStretch(1)

        frame_right.setSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        frame_right.setLayout(self.pandda_analyse_input_params_vbox)

        grid_pandda.addWidget(frame_right, 0, 1, 5, 5)
        frame_pandda.setLayout(grid_pandda)
        self.pandda_analyse_hbox.addWidget(frame_pandda)

        # next three blocks display html documents created by pandda.analyse
        self.layout_funcs.pandda_html(self)

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

        self.panddas_results_vbox.addWidget(self.pandda_tab_widget)
        self.show_pandda_html_summary()

        ################################################################################################################
        #                                                                                                              #
        #                                                REFINEMENT TAB                                                #
        #                                                                                                              #
        ################################################################################################################
        self.summary_vbox_for_table = QtGui.QVBoxLayout()

        # table
        self.refinement_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(self.refinement_table, self.refinement_table_columns)
        self.summary_vbox_for_table.addWidget(self.refinement_table)
        self.tab_dict[self.workflow_dict['Refinement']][1].addLayout(self.summary_vbox_for_table)

        ################################################################################################################
        #                                                                                                              #
        #                                                DEPOSITION TAB                                                #
        #                                                                                                              #
        ################################################################################################################
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
        label_heading = QtGui.QLabel('1. Specify HTML export directory in the settings tab')
        label_heading.setStyleSheet("font: bold 20pt Arial")
        scrollLayout.addWidget(label_heading)
        label_text = QtGui.QLabel(XChemToolTips.html_export_directory_background())
        label_text.setStyleSheet("font: 17pt Arial")
        scrollLayout.addWidget(label_text)
        if os.getcwd().startswith('/dls'):
            label_text = QtGui.QLabel(
                'Note: default for labxchem project at DLS is <labxchem_directory>/processing/html.')
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
        button.clicked.connect(self.prepare_files_for_zenodo_upload)
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
        scroll.setWidget(scrollContent)
        self.tab_dict[self.workflow_dict['Deposition']][1].addLayout(self.deposition_vbox)

        ################################################################################################################
        #                                                                                                              #
        #                                                 SETTINGS TAB                                                 #
        #                                                                                                              #
        ################################################################################################################
        self.settings_container = QtGui.QWidget()
        self.buttons_etc = QtGui.QWidget()
        self.settings_vbox = QtGui.QVBoxLayout()

        self.scroll = QtGui.QScrollArea(self.settings_container)
        self.settings_vbox.addWidget(self.scroll)
        scrollContent_settings = QtGui.QWidget(scroll)

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
        settings_buttoon_initial_model_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_initial_model_directory.addWidget(settings_buttoon_initial_model_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_initial_model_directory)

        self.data_collection_vbox_for_settings.addWidget(
            QtGui.QLabel('\n\nReference Structure Directory: - OPTIONAL -'))
        settings_hbox_reference_directory = QtGui.QHBoxLayout()
        self.reference_directory_label = QtGui.QLabel(self.reference_directory)
        settings_hbox_reference_directory.addWidget(self.reference_directory_label)
        settings_buttoon_reference_directory = QtGui.QPushButton('Select Reference Structure Directory')
        settings_buttoon_reference_directory.clicked.connect(self.settings_button_clicked)
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
        settings_buttoon_data_source_file.clicked.connect(self.settings_button_clicked)
        settings_hbox_data_source_file.addWidget(settings_buttoon_data_source_file)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_data_source_file)

        ## Data Collection
        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nData Collection Directory: - OPTIONAL -'))

        settings_beamline_frame = QtGui.QFrame()
        settings_beamline_frame.setFrameShape(QtGui.QFrame.StyledPanel)
        settings_beamline_vbox = QtGui.QVBoxLayout()

        settings_hbox_beamline_directory = QtGui.QHBoxLayout()
        self.beamline_directory_label = QtGui.QLabel(self.beamline_directory)
        settings_hbox_beamline_directory.addWidget(self.beamline_directory_label)
        settings_buttoon_beamline_directory = QtGui.QPushButton('Select Data Collection Directory')
        settings_buttoon_beamline_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_beamline_directory.addWidget(settings_buttoon_beamline_directory)
        settings_beamline_vbox.addLayout(settings_hbox_beamline_directory)

        settings_hbox_datasets_summary_file = QtGui.QHBoxLayout()
        self.datasets_summary_file_label = QtGui.QLabel(self.datasets_summary_file)
        settings_hbox_datasets_summary_file.addWidget(self.datasets_summary_file_label)
        settings_button_datasets_summary_file = QtGui.QPushButton('Select Existing\nCollection Summary File')
        settings_button_datasets_summary_file.clicked.connect(self.settings_button_clicked)
        settings_hbox_datasets_summary_file.addWidget(settings_button_datasets_summary_file)

        settings_button_new_datasets_summary_file = QtGui.QPushButton('Assign New\nCollection Summary File')
        settings_button_new_datasets_summary_file.clicked.connect(self.settings_button_clicked)
        settings_hbox_datasets_summary_file.addWidget(settings_button_new_datasets_summary_file)

        settings_beamline_vbox.addLayout(settings_hbox_datasets_summary_file)

        settings_beamline_frame.setLayout(settings_beamline_vbox)
        self.data_collection_vbox_for_settings.addWidget(settings_beamline_frame)

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nCCP4_SCR Directory: - OPTIONAL -'))
        settings_hbox_ccp4_scratch_directory = QtGui.QHBoxLayout()
        self.ccp4_scratch_directory_label = QtGui.QLabel(self.ccp4_scratch_directory)
        settings_hbox_ccp4_scratch_directory.addWidget(self.ccp4_scratch_directory_label)
        settings_buttoon_ccp4_scratch_directory = QtGui.QPushButton('Select CCP4_SCR Directory')
        settings_buttoon_ccp4_scratch_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_ccp4_scratch_directory.addWidget(settings_buttoon_ccp4_scratch_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_ccp4_scratch_directory)

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nPANDDAs directory: - OPTIONAL -'))
        settings_hbox_panddas_directory = QtGui.QHBoxLayout()
        self.panddas_directory_label = QtGui.QLabel(self.panddas_directory)
        settings_hbox_panddas_directory.addWidget(self.panddas_directory_label)
        settings_button_panddas_directory = QtGui.QPushButton('Select PANNDAs Directory')
        settings_button_panddas_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_panddas_directory.addWidget(settings_button_panddas_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_panddas_directory)

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nHTML export directory: - OPTIONAL -'))
        settings_hbox_html_export_directory = QtGui.QHBoxLayout()
        self.html_export_directory_label = QtGui.QLabel(self.html_export_directory)
        settings_hbox_html_export_directory.addWidget(self.html_export_directory_label)
        settings_button_html_export_directory = QtGui.QPushButton('Select HTML Export Directory')
        settings_button_html_export_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_html_export_directory.addWidget(settings_button_html_export_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_html_export_directory)

        self.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nGroup deposition directory: - OPTIONAL -'))
        settings_hbox_group_deposition_directory = QtGui.QHBoxLayout()
        self.group_deposition_directory_label = QtGui.QLabel(self.group_deposit_directory)
        settings_hbox_group_deposition_directory.addWidget(self.group_deposition_directory_label)
        settings_button_group_deposition_directory = QtGui.QPushButton('Select Group deposition Directory')
        settings_button_group_deposition_directory.clicked.connect(self.settings_button_clicked)
        settings_hbox_group_deposition_directory.addWidget(settings_button_group_deposition_directory)
        self.data_collection_vbox_for_settings.addLayout(settings_hbox_group_deposition_directory)

        self.data_collection_vbox_for_settings.setSpacing(0)
        self.data_collection_vbox_for_settings.setContentsMargins(30, 0, 0, 0)

        self.buttons_etc.resize(self.screen.width() - 100, self.buttons_etc.sizeHint().height())
        self.tab_dict[self.workflow_dict['Settings']][1].addLayout(self.settings_vbox)

        ################################################################################################################
        #                                                                                                              #
        #                                                  STATUS BAR                                                  #
        #                                                                                                              #
        ################################################################################################################
        self.status_bar = QtGui.QStatusBar()
        self.progress_bar = QtGui.QProgressBar()
        self.progress_bar.setMaximum(100)
        self.status_bar.setMaximumWidth(self.screen.width())
        self.progress_bar.setMaximumWidth(self.screen.width())
        hbox_status = QtGui.QHBoxLayout()
        hbox_status.addWidget(self.status_bar)
        hbox_status.addWidget(self.progress_bar)

        vbox_main = QtGui.QVBoxLayout()
        menu_bar.setMaximumWidth(self.screen.width())
        vbox_main.addWidget(menu_bar)
        self.main_tab_widget.setMaximumSize(self.screen.width(), self.screen.height() - 245)
        vbox_main.addWidget(self.main_tab_widget)

        hboxTaskFrames = QtGui.QHBoxLayout()

        hboxTaskFrames.addWidget(update_from_datasource_button)
        hboxTaskFrames.addWidget(frame_dataset_task)
        hboxTaskFrames.addWidget(frame_map_cif_file_task)
        hboxTaskFrames.addWidget(frame_panddas_file_task)
        hboxTaskFrames.addWidget(frame_refine_file_task)

        vbox_main.addLayout(hboxTaskFrames)

        vbox_main.addLayout(hbox_status)

        self.window.setLayout(vbox_main)

        self.status_bar.showMessage('Ready')
        self.window.show()

        if self.data_source_file != '':
            write_enabled = self.check_write_permissions_of_data_source()
            if not write_enabled:
                self.data_source_set = False

    def workflow(self):
        ################################################################################################################
        #                                                                                                              #
        # ========================================== WORKFLOW TASK CONTAINER ========================================= #
        #                                                                                                              #
        ################################################################################################################
        self.workflow_widget_dict = {}

        # workflow task container - order of tabs as they appear for the main window
        self.workflow = ['Overview',  # 0
                         'Datasets',  # 1
                         'Maps',  # 2
                         'PANDDAs',  # 3
                         'Refinement',  # 4
                         'Deposition',  # 6
                         'Settings']  # 5

        # dictionary with keys corresponding to each stage in the workflow
        self.workflow_dict = {self.workflow[0]: 'Overview',
                              self.workflow[1]: 'Datasets',
                              self.workflow[2]: 'Maps',
                              self.workflow[3]: 'PANDDAs',
                              self.workflow[4]: 'Refinement',
                              self.workflow[6]: 'Settings',
                              self.workflow[5]: 'Deposition'}

        # tab widget
        self.main_tab_widget = QtGui.QTabWidget()
        self.tab_dict = {}
        self.layout_funcs.make_tab_dict(self.workflow, self.main_tab_widget, self.tab_dict)

class LayoutFuncs():
    def __init__(self):
        pass

    def make_tab_dict(self, tab_list, tab_widget, tab_dict):
        for page in tab_list:
            tab = QtGui.QWidget()
            vbox = QtGui.QVBoxLayout(tab)
            tab_widget.addTab(tab, page)
            tab_dict[page] = [tab, vbox]

    def add_checkbox(self, object, checkbox, function, checkopt=False):
        checkbox.toggle()
        checkbox.setChecked(checkopt)
        eval(str('checkbox.stateChanged.connect(' + function + ')'))

    def table_setup(self, table, table_columns, sortingopt=True):
        table.setColumnCount(len(table_columns))
        table.setSortingEnabled(sortingopt)
        table.setHorizontalHeaderLabels(table_columns)
        table.resizeRowsToContents()
        table.resizeColumnsToContents()

    def pandda_html(self, object):
        if os.path.exists(str(object.panddas_directory + '/interesting_datasets')):
            print('WARNING: USING RESULTS FROM OLD PANDDA ANALYSE! THIS IS NOT FULLY SUPPORTED IN XCE2')
            print('PLEASE CHANGE YOUR PANDDA DIRECTORY TO A NEW RUN, OR USE THE OLD VERSION OF XCE!')
            object.pandda_initial_html_file = str(object.panddas_directory + '/results_summareis/pandda_initial.html')
            object.pandda_analyse_html_file = str(object.panddas_directory + '/results_summaries/pandda_analyse.html')
        object.pandda_initial_html_file = str(
            object.panddas_directory + '/analyses/html_summaries/' + 'pandda_initial.html')
        object.pandda_analyse_html_file = str(
            object.panddas_directory + '/analyses/html_summaries/' + 'pandda_analyse.html')
        object.pandda_inspect_html_file = str(
            object.panddas_directory + '/analyses/html_summaries/' + 'pandda_inspect.html')

    # function for datasource, run and status button setup
    def setup_push_button(self, object, button_dict):
        # use iterkeys to determine order of key by letter
        for name in sorted(button_dict.iterkeys()):
            # add current item to menu bar
            button = eval('QtGui.QPushButton("' + str(button_dict[name][0]) + '")')
            # for each configuration item
            for button_config in button_dict[name][1]:
                eval(str('button.setToolTip(' + str(button_config[0]) + ')'))
                eval(str('button.setStyleSheet("' + str(button_config[1] + '")')))
                if len(button_config[2]) > 1:
                    eval(str('button.setFont(' + str(button_config[2]) + ')'))
                eval(str('button.clicked.connect(' + str(button_config[3]) + ')'))

        return button

    # function to setup one of the bottom boxes
    def bottom_box_setup(self, object, label, dropdown_options, dropdown_tooltip, buttons, colour):

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)
        frame.setStyleSheet("QFrame { border: 1px solid black; border-radius: 1px; padding: 0px; margin: 0px }")

        vbox = QtGui.QVBoxLayout()
        label = QtGui.QLabel(label)
        label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        label.setFont(object.headlineLabelfont)
        label.setStyleSheet(str(" QLabel { border: 1px solid black; border-radius: 1px;" + str(colour) +
                                "padding: 0px; margin: 0px }"))
        vbox.addWidget(label)

        hboxAction = QtGui.QHBoxLayout()
        combobox = QtGui.QComboBox()
        for task in dropdown_options:
            combobox.addItem(task)
        eval('combobox.setToolTip(' + str(dropdown_tooltip) + ')')
        combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")
        hboxAction.addWidget(combobox)

        vboxButton = QtGui.QVBoxLayout()
        for button in buttons:
            vboxButton.addWidget(button)
        hboxAction.addLayout(vboxButton)
        vbox.addLayout(hboxAction)
        vbox.setSpacing(0)
        vbox.setMargin(0)
        frame.setLayout(vbox)
        frame.setMaximumWidth((object.screen.width() - 20) / 5)

        return frame, combobox

    # function to add items to top menu bar
    def setup_menubar(self, object, menu_bar, menu_items_dict):
        # use iterkeys to determine order of key by letter
        for config in sorted(menu_items_dict.iterkeys()):
            # add current item to menu bar
            menu = eval('menu_bar.addMenu("' + str(menu_items_dict[config][0]) + '")')
            # for each configuration item
            for menu_item in menu_items_dict[config][1]:
                # add the drop down option
                action = eval(str('QtGui.QAction("' + str(menu_item[0]) + '", object.window)'))
                # add a shortcut if defined
                if len(menu_item[1]) > 1:
                    eval(str('action.setShortcut("' + str(menu_item[1]) + '")'))
                # connect the relevant function and add as an action
                eval(str('action.triggered.connect(' + menu_item[2] + ')'))
                menu.addAction(action)

        return menu_bar

    # function for opening help and tutorial files
    def openFile(self, file):
        if sys.platform == 'linux2':
            subprocess.call(["xdg-open", file])
        else:
            os.startfile(file)

    def add_to_box(self, frame, widgets_list):
        for widget in widgets_list:
            frame.addWidget(widget)
