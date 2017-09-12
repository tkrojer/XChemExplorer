import multiprocessing
import subprocess
import sys, os
from PyQt4 import QtGui, QtCore, QtWebKit

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'lib'))
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'web'))
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'gui_scripts'))

from settings_preferences import *

import XChemToolTips
import XChemMain

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas


class LayoutObjects():
    def __init__(self, object):
        self.layout_funcs = LayoutFuncs()

    # function to initialise the top menu bar
    def initialise_menu_bar(self, object):
        ################################################################################################################
        #                                                                                                              #
        #                                              MENU BAR - TOP OF GUI                                           #
        #                                                                                                              #
        ################################################################################################################

        # initiate menu widget
        menu_bar = QtGui.QMenuBar()

        # import menu bar dictionary
        setup().top_menu_dict(self)

        # create menu from menu dictionary
        menu_bar = self.layout_funcs.setup_menubar(object, menu_bar, self.menu_dict)

        # END OF MENU BAR - CODE BELOW: stuff removed from apo structure stuff that appears might have a funky
        # consequence - work out later.

        # object.prepare_mmcif_files_dict = {}
        # object.prepare_mmcif_files_dict['apo'] = prepare_mmcif_files_for_apo_structures
        # object.prepare_mmcif_files_dict['ligand_bound'] = prepare_mmcif_files_for_ligand_bound_structures

        return menu_bar

    # function containing setup for bottom boxes
    def initialise_bottom_boxes(self, object):

        # import all buttons
        setup().bottom_box_buttons(object)

        # setup datasource button
        update_from_datasource_button = self.layout_funcs.setup_push_button(object, object.datasource_button_dict)

        ################################################################################################################
        #                                                                                                              #
        #                                               DATASETS BOX                                                   #
        #                                                                                                              #
        ################################################################################################################

        # setup the run button with push button function
        object.dataset_task_run_button = self.layout_funcs.setup_push_button(object,
                                                                             object.dataset_task_run_button_dict)

        # setup the task button with push button function
        object.dataset_task_status_button = self.layout_funcs.setup_push_button(object,
                                                                                object.dataset_task_status_button_dict)

        # array of both button objects to apply to bottom box layout
        dataset_buttons = [object.dataset_task_run_button, object.dataset_task_status_button]

        # label for the bottom box layout
        dataset_label = "Datasets"

        # return the frame and combobox from the bottom box setup function
        frame_dataset_task, object.dataset_tasks_combobox = self.layout_funcs.bottom_box_setup(object, dataset_label,
                                                                                               object.dataset_tasks,
                                                                                               'XChemToolTips.'
                                                                                               'dataset_task_tip()',
                                                                                               dataset_buttons,
                                                                                               'background: '
                                                                                               'rgb(240, 255, 140); ')

        # define the combobox and buttons in dictionary key to determine behaviour
        object.workflow_widget_dict['Datasets'] = [object.dataset_tasks_combobox, object.dataset_task_run_button,
                                                   object.dataset_task_status_button]

        ################################################################################################################
        #                                                                                                              #
        #                                            MAPS & RESTRAINTS BOX                                             #
        #                                                                                                              #
        ################################################################################################################
        # settings for the run button

        # setup the run button with push button function
        object.map_cif_file_task_run_button = \
            self.layout_funcs.setup_push_button(object,
                                                object.map_cif_file_task_run_button_dict)

        # setup the task button with push button function
        object.map_cif_file_task_status_button = \
            self.layout_funcs.setup_push_button(object,
                                                object.map_cif_file_task_status_button_dict)

        # array of both button objects to apply to bottom box layout
        map_cif_file_buttons = [object.map_cif_file_task_run_button, object.map_cif_file_task_status_button]

        # label for the bottom box layout
        map_cif_file_label = "Maps & Restraints"

        # return the frame and combobox from the bottom box setup function
        frame_map_cif_file_task, object.map_cif_file_tasks_combobox = \
            self.layout_funcs.bottom_box_setup(object,
                                               map_cif_file_label,
                                               object.map_cif_file_tasks,
                                               'XChemToolTips.map_cif_file_'
                                               'task_tip()',
                                               map_cif_file_buttons,
                                               'background: rgb(140, 255, '
                                               '150); ')

        # define the combobox and buttons in dictionary key to determine behaviour
        object.workflow_widget_dict['Maps'] = [object.map_cif_file_tasks_combobox, object.map_cif_file_task_run_button,
                                               object.map_cif_file_task_status_button]

        ################################################################################################################
        #                                                                                                              #
        #                                               HIT IDENTIFICATION BOX                                         #
        #                                                                                                              #
        ################################################################################################################
        # settings for the run button

        # setup the run button with push button function
        object.panddas_file_task_run_button = \
            self.layout_funcs.setup_push_button(object,
                                                object.panddas_file_task_run_button_dict)

        # setup the task button with push button function
        object.panddas_file_task_status_button = \
            self.layout_funcs.setup_push_button(object,
                                                object.panddas_file_task_status_button_dict)

        # array of both button objects to apply to bottom box layout
        panddas_file_buttons = [object.panddas_file_task_run_button, object.panddas_file_task_status_button]

        # label for the bottom box layout
        panddas_file_label = "Hit Identification"

        # return the frame and combobox from the bottom box setup function
        frame_panddas_file_task, object.panddas_file_tasks_combobox = \
            self.layout_funcs.bottom_box_setup(object,
                                               panddas_file_label,
                                               object.panddas_file_tasks,
                                               'XChemToolTips.panddas_file_'
                                               'task_tip()',
                                               panddas_file_buttons,
                                               'background: rgb(140,200,255)'
                                               '; ')

        # define the combobox and buttons in dictionary key to determine behaviour
        object.workflow_widget_dict['PANDDAs'] = [object.panddas_file_tasks_combobox,
                                                  object.panddas_file_task_run_button,
                                                  object.panddas_file_task_status_button]

        ################################################################################################################
        #                                                                                                              #
        #                                                      REFINEMENT BOX                                          #
        #                                                                                                              #
        ################################################################################################################
        # settings for the run button

        # setup the run button with push button function
        object.refine_file_task_run_button = \
            self.layout_funcs.setup_push_button(object, object.refine_file_task_run_button_dict)

        # setup the task button with push button function
        object.refine_file_task_status_button = \
            self.layout_funcs.setup_push_button(object, object.refine_file_task_status_button_dict)

        # array of both button objects to apply to bottom box layout
        refine_file_buttons = [object.refine_file_task_run_button, object.refine_file_task_status_button]

        # label for the bottom box layout
        refine_file_label = "Refinement"

        # return the frame and combobox from the bottom box setup function
        frame_refine_file_task, object.refine_file_tasks_combobox = \
            self.layout_funcs.bottom_box_setup(object,
                                               refine_file_label,
                                               object.refine_file_tasks,
                                               'XChemToolTips.refine_file_task'
                                               '_tip()',
                                               refine_file_buttons,
                                               'background: rgb(245, 190, 255)'
                                               ';')

        # define the combobox and buttons in dictionary key to determine behaviour
        object.workflow_widget_dict['Refinement'] = [object.refine_file_tasks_combobox,
                                                     object.refine_file_task_run_button,
                                                     object.refine_file_task_status_button]

        return update_from_datasource_button, frame_dataset_task, frame_map_cif_file_task, frame_panddas_file_task, \
               frame_refine_file_task

    def main_layout(self, object):
        # initialise menu bar
        menu_bar = self.initialise_menu_bar(object)

        # initialise bottom boxes
        update_from_datasource_button, frame_dataset_task, frame_map_cif_file_task, frame_panddas_file_task, \
        frame_refine_file_task = self.initialise_bottom_boxes(object)

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

        ################################################################################################################
        #                                                                                                              #
        #                                                 OVERVIEW TAB                                                 #
        #                                                                                                              #
        ################################################################################################################
        # define subtab list, widget and dict
        overview_tab_list = ['Data Source', 'Summary']
        object.overview_tab_widget = QtGui.QTabWidget()
        object.overview_tab_dict = {}

        # make subtabs
        self.layout_funcs.make_tab_dict(overview_tab_list, object.overview_tab_widget, object.overview_tab_dict)

        # add overview subtabs to overview tab
        object.tab_dict[object.workflow_dict['Overview']][1].addWidget(object.overview_tab_widget)

        # initiate the table in overview/datasource
        object.overview_datasource_table = QtGui.QTableWidget()
        object.overview_datasource_table.setSortingEnabled(True)
        object.overview_datasource_table.resizeColumnsToContents()

        # add table to datasource subtab
        object.overview_tab_dict['Data Source'][1].addWidget(object.overview_datasource_table)

        # initiate the graph in overview/summary
        object.overview_figure, object.overview_axes = plt.subplots()
        object.overview_canvas = FigureCanvas(object.overview_figure)
        object.update_summary_plot()

        # add graph to summary tab
        object.overview_tab_dict['Summary'][1].addWidget(object.overview_canvas)

        ################################################################################################################
        #                                                                                                              #
        #                                                 DATASETS TAB                                                 #
        #                                                                                                              #
        ################################################################################################################
        # define subtab list, widget and dict
        datasets_tab_list = ['Summary', 'Reprocess']
        object.datasets_tab_widget = QtGui.QTabWidget()
        object.datasets_tab_dict = {}

        # make subtabs
        self.layout_funcs.make_tab_dict(datasets_tab_list, object.datasets_tab_widget, object.datasets_tab_dict)

        # main body - things that are always displayed
        # add a container to hold everythting and add to main tab layout
        object.datasets_data_collection_vbox = QtGui.QVBoxLayout()
        object.tab_dict[object.workflow_dict['Datasets']][1].addLayout(object.datasets_data_collection_vbox)

        # add a horizontal box to hold option to autocheck for new data
        object.autocheck_hbox = QtGui.QHBoxLayout()

        # checkbox for autocollect
        object.check_for_new_data_collection = QtGui.QCheckBox('Check for new data collection every two minutes')
        self.layout_funcs.add_checkbox(object, object.check_for_new_data_collection,
                                       'object.continously_check_for_new_data_collection')

        # select target dropdown
        select_target_label = QtGui.QLabel('Select Target: ')
        select_target_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        object.target_selection_combobox = QtGui.QComboBox()
        object.populate_target_selection_combobox(object.target_selection_combobox)
        object.target_selection_combobox.activated[str].connect(object.target_selection_combobox_activated)
        object.target = str(object.target_selection_combobox.currentText())

        object.autocheck_hbox_widgets = [object.check_for_new_data_collection, select_target_label,
                                         object.target_selection_combobox]  # array defining order of objects to add

        self.layout_funcs.add_to_box(object.autocheck_hbox, object.autocheck_hbox_widgets)  # add objects in order

        # add target dropdown to top bar
        object.datasets_data_collection_vbox.addLayout(object.autocheck_hbox)

        # summary sub-tab
        # table
        object.datasets_summary_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(object.datasets_summary_table, object.datasets_summary_table_columns)
        object.datasets_summarys_vbox_for_table = QtGui.QVBoxLayout()  # setup layout to hold table
        object.datasets_summarys_vbox_for_table.addWidget(object.datasets_summary_table)  # add table to layout
        object.datasets_tab_dict['Summary'][1].addLayout(object.datasets_summarys_vbox_for_table)  # add layout to tab
        object.datasets_summarys_vbox_for_details = QtGui.QVBoxLayout()  # vbox for details
        object.data_collection_details_currently_on_display = None  # switch for displaying/updating table
        object.datasets_tab_dict['Summary'][1].addLayout(object.datasets_summarys_vbox_for_details)  # add details

        object.datasets_data_collection_vbox.addWidget(object.datasets_tab_widget)  # add subtab to main tab

        # reprocessing sub-tab
        # top options
        # data collection label
        dc_label = QtGui.QLabel('Data collection directory: ')
        dc_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)  # align left and centre of container
        object.diffraction_data_dir_label = QtGui.QLabel(object.diffraction_data_directory)  # add directory as text
        object.diffraction_data_dir_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)  # align as above

        # select label
        select_button = QtGui.QPushButton("Select")
        select_button.clicked.connect(object.select_diffraction_data_directory)  # attach file open dialogue

        # search button
        search_button = QtGui.QPushButton("Search Datasets")
        search_button.clicked.connect(object.search_for_datasets)  # search for datasets in the selected directory

        # search info
        object.diffraction_data_search_label = QtGui.QLabel(object.diffraction_data_search_info)

        # translate label
        translate_label = QtGui.QLabel('translate: datasetID -> sampleID')
        translate_label.setAlignment(QtCore.Qt.AlignCenter)  # align in centre of container

        # CSV button
        csv_button = QtGui.QPushButton('Open CSV')
        csv_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
        csv_button.clicked.connect(object.translate_datasetID_to_sampleID)  # open the relevant csv file

        # create hbox to hold everything and add widgets to it
        object.hbox_select = QtGui.QHBoxLayout()  # top options box
        object.hbox_select_widgets = [dc_label, select_button, search_button, object.diffraction_data_search_label,
                                      translate_label, csv_button]  # array defining order of objects to be added
        self.layout_funcs.add_to_box(object.hbox_select, object.hbox_select_widgets)  # add objects in order

        # frame to hold everything
        frame_select = QtGui.QFrame()
        frame_select.setLayout(object.hbox_select)  # apply to containing frame

        # table - main body
        object.datasets_reprocess_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(object.datasets_reprocess_table, object.datasets_reprocess_columns)  # setup

        # create context menu - no idea where this lives...
        object.popMenu_for_datasets_reprocess_table = QtGui.QMenu()
        run_xia2_on_selected = QtGui.QAction("mark selected for reprocessing", object.window)
        run_xia2_on_selected.triggered.connect(object.select_sample_for_xia2)
        object.popMenu_for_datasets_reprocess_table.addAction(run_xia2_on_selected)
        object.datasets_reprocess_table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        object.datasets_reprocess_table.customContextMenuRequested.connect(object.on_context_menu_reprocess_data)

        # options at bottom of tab
        # data processing label
        label = QtGui.QLabel('Data processing protocol: ')
        label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)

        # option checkboxes
        object.xia2_3d_checkbox = QtGui.QCheckBox('xia2 3d')
        object.xia2_3dii_checkbox = QtGui.QCheckBox('xia2 3dii')
        object.xia2_dials_checkbox = QtGui.QCheckBox('Dials')

        # spacegroup label
        sg_label = QtGui.QLabel('Space Group:')

        # spacegroup dropdown menu
        object.reprocess_space_group_comboxbox = QtGui.QComboBox()
        object.reprocess_space_group_comboxbox.addItem('ignore')
        for sg in XChemMain.space_group_list():
            object.reprocess_space_group_comboxbox.addItem(sg)

        # mtz label
        mtz_label = QtGui.QLabel('Reference MTZ:')

        # file label
        object.reprocess_reference_mtz_file_label = QtGui.QLabel(object.diffraction_data_reference_mtz)

        # select button
        select_button = QtGui.QPushButton("Select")
        select_button.clicked.connect(object.select_reprocess_reference_mtz)

        # define order of widgets to be added to options hbox
        hbox_options_widgets = [label, object.xia2_3d_checkbox, object.xia2_3dii_checkbox,
                                object.xia2_dials_checkbox, sg_label, object.reprocess_space_group_comboxbox,
                                mtz_label, object.reprocess_reference_mtz_file_label, select_button]

        # create hbox, add everything to it and then put it in a frame
        hbox_options = QtGui.QHBoxLayout()
        self.layout_funcs.add_to_box(hbox_options, hbox_options_widgets)

        frame_options = QtGui.QFrame()
        frame_options.setLayout(hbox_options)

        # following are contained in vboxes
        # res limit isig label
        label = QtGui.QLabel('Resolution\nLimit:\nMn<I/sig(I)>')
        label.setAlignment(QtCore.Qt.AlignCenter)

        # res limit isig dropdown menu
        object.reprocess_isigma_combobox = QtGui.QComboBox()
        misigma = ['default', '4', '3', '2.5', '2', '1.5', '1', '0.5']
        for item in misigma:
            object.reprocess_isigma_combobox.addItem(item)
        object.reprocess_isigma_combobox.setCurrentIndex(0)
        object.reprocess_isigma_combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")

        # create vertical box to add labels and dropdowns to, create box and put in frame
        vbox_isigma = QtGui.QVBoxLayout()
        vbox_isigma_widgets = [label, object.reprocess_isigma_combobox]
        self.layout_funcs.add_to_box(vbox_isigma, vbox_isigma_widgets)
        frame_isigma = QtGui.QFrame()
        frame_isigma.setLayout(vbox_isigma)

        # res limit cc half label
        res_cc_label = QtGui.QLabel('Resolution\nLimit:\nCC 1/2')
        label.setAlignment(QtCore.Qt.AlignCenter)

        # res limit cc half dropdown
        object.reprocess_cc_half_combobox = QtGui.QComboBox()
        cc_half = ['default', '0.9', '0.8', '0.7', '0.6', '0.5', '0.4', '0.3', '0.2', '0.1']
        for item in cc_half:
            object.reprocess_cc_half_combobox.addItem(item)
        object.reprocess_cc_half_combobox.setCurrentIndex(0)
        object.reprocess_cc_half_combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")

        # create a vbox for label and dropdown, and add items to it
        vbox_cc_half = QtGui.QVBoxLayout()
        vbox_cc_half_widgets = [res_cc_label, object.reprocess_cc_half_combobox]
        self.layout_funcs.add_to_box(vbox_cc_half, vbox_cc_half_widgets)

        # create frame to hold everything and add vbox
        frame_cc_half = QtGui.QFrame()
        frame_cc_half.setLayout(vbox_cc_half)

        # create a hbox to hold the bottom frames and add everything
        data_protocol_hbox = QtGui.QHBoxLayout()
        data_protocol_hbox_widgets = [frame_options, frame_isigma, frame_cc_half]
        self.layout_funcs.add_to_box(data_protocol_hbox, data_protocol_hbox_widgets)

        bottom_options_frame = QtGui.QFrame()  # create frame to hold everything (horizontal)
        bottom_options_frame.setLayout(data_protocol_hbox)

        # code below sets final layout for whole subtab
        reprocess_vbox = QtGui.QVBoxLayout()  # box to hold reprocessing subtab content
        object.reprocess_hbox_widgets = [frame_select, object.datasets_reprocess_table, bottom_options_frame]
        self.layout_funcs.add_to_box(reprocess_vbox, object.reprocess_hbox_widgets)

        object.datasets_tab_dict['Reprocess'][1].addLayout(reprocess_vbox)

        ################################################################################################################
        #                                                                                                              #
        #                                                    MAPS TAB                                                  #
        #                                                                                                              #
        ################################################################################################################
        initial_model_checkbutton_hbox = QtGui.QHBoxLayout()

        object.select_sample_for_dimple_box = QtGui.QCheckBox('(de-)select all samples for DIMPLE')
        self.layout_funcs.add_checkbox(object, object.select_sample_for_dimple_box, 'object.set_run_dimple_flag')
        initial_model_checkbutton_hbox.addWidget(object.select_sample_for_dimple_box)

        set_new_reference_button = QtGui.QPushButton("Set New Reference (if applicable)")
        set_new_reference_button.clicked.connect(object.set_new_reference_if_applicable)
        initial_model_checkbutton_hbox.addWidget(set_new_reference_button)

        refresh_reference_file_list_button = QtGui.QPushButton("Refresh reference file list")
        refresh_reference_file_list_button.clicked.connect(object.refresh_reference_file_list)
        initial_model_checkbutton_hbox.addWidget(refresh_reference_file_list_button)

        object.reference_file_list = object.get_reference_file_list(' ')
        object.reference_file_selection_combobox = QtGui.QComboBox()
        object.populate_reference_combobox(object.reference_file_selection_combobox)
        initial_model_checkbutton_hbox.addWidget(object.reference_file_selection_combobox)

        object.tab_dict[object.workflow_dict['Maps']][1].addLayout(initial_model_checkbutton_hbox)
        object.initial_model_vbox_for_table = QtGui.QVBoxLayout()

        # table
        object.maps_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(object.maps_table, object.maps_table_columns)
        object.initial_model_vbox_for_table.addWidget(object.maps_table)
        object.tab_dict[object.workflow_dict['Maps']][1].addLayout(object.initial_model_vbox_for_table)

        # create context menu
        object.popMenu_for_maps_table = QtGui.QMenu()
        run_dimple = QtGui.QAction("mark selected for dimple run", object.window)
        run_dimple.triggered.connect(object.select_sample_for_dimple)
        object.popMenu_for_maps_table.addAction(run_dimple)
        object.maps_table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        object.maps_table.customContextMenuRequested.connect(object.on_context_menu_initial_model)

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

        object.pandda_tab_widget = QtGui.QTabWidget()
        object.pandda_tab_dict = {}
        self.layout_funcs.make_tab_dict(pandda_tab_list, object.pandda_tab_widget, object.pandda_tab_dict)

        ## PanDDA tab
        object.panddas_results_vbox = QtGui.QVBoxLayout()
        object.tab_dict[object.workflow_dict['PANDDAs']][1].addLayout(object.panddas_results_vbox)

        object.pandda_analyse_hbox = QtGui.QHBoxLayout()
        object.pandda_tab_dict['pandda.analyse'][1].addLayout(object.pandda_analyse_hbox)
        object.pandda_map_layout = QtGui.QVBoxLayout()
        object.pandda_map_list = QtGui.QComboBox()
        object.pandda_maps_html = QtWebKit.QWebView()
        object.pandda_map_layout.addWidget(object.pandda_map_list)
        object.pandda_map_layout.addWidget(object.pandda_maps_html)

        object.pandda_tab_dict['Statistical Map Summaries'][1].addLayout(object.pandda_map_layout)
        object.pandda_maps_html.show()

        grid_pandda = QtGui.QGridLayout()
        grid_pandda.setColumnStretch(0, 20)
        grid_pandda.setRowStretch(0, 20)
        # left hand side: table with information about available datasets

        # table
        object.pandda_analyse_data_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(object.pandda_analyse_data_table, object.pandda_table_columns)

        frame_pandda = QtGui.QFrame()
        grid_pandda.addWidget(object.pandda_analyse_data_table, 0, 0)

        object.pandda_status = 'UNKNOWN'
        object.pandda_status_label = QtGui.QLabel()

        pandda_status_options = [['/pandda.done', 'Finished!', 'color: green'],
                                 ['/pandda.running', 'Running...', 'color: orange'],
                                 ['/pandda.errored', 'Error encountered... please check the log files for pandda!',
                                  'color: red']]

        for option in pandda_status_options:
            if os.path.exists(str(object.panddas_directory + option[0])):
                object.pandda_status = option[1]
                object.pandda_status_label.setStyleSheet(option[2])

        object.pandda_status_label.setText(str('STATUS: ' + object.pandda_status))
        object.pandda_status_label.setFont(QtGui.QFont("Arial", 25, QtGui.QFont.Bold))
        grid_pandda.addWidget(object.pandda_status_label, 3, 0)

        # right hand side: input parameters for PANDDAs run
        frame_right = QtGui.QFrame()
        frame_right.setFrameShape(QtGui.QFrame.StyledPanel)

        object.pandda_analyse_input_params_vbox = QtGui.QVBoxLayout()

        pandda_input_dir_hbox = QtGui.QHBoxLayout()
        label = QtGui.QLabel('data directory')
        object.pandda_analyse_input_params_vbox.addWidget(label)
        object.pandda_input_data_dir_entry = QtGui.QLineEdit()
        object.pandda_input_data_dir_entry.setText(os.path.join(object.initial_model_directory, '*'))
        object.pandda_input_data_dir_entry.setFixedWidth(300)
        pandda_input_dir_hbox.addWidget(object.pandda_input_data_dir_entry)
        object.select_pandda_input_dir_button = QtGui.QPushButton("Select Input Template")
        object.select_pandda_input_dir_button.setMaximumWidth(200)
        object.select_pandda_input_dir_button.clicked.connect(object.select_pandda_input_template)
        pandda_input_dir_hbox.addWidget(object.select_pandda_input_dir_button)
        object.pandda_analyse_input_params_vbox.addLayout(pandda_input_dir_hbox)

        pandda_pdb_style_hbox = QtGui.QHBoxLayout()
        label = QtGui.QLabel('pdb style')
        pandda_pdb_style_hbox.addWidget(label)
        object.pandda_pdb_style_entry = QtGui.QLineEdit()
        object.pandda_pdb_style_entry.setText('dimple.pdb')
        object.pandda_pdb_style_entry.setFixedWidth(200)
        pandda_pdb_style_hbox.addWidget(object.pandda_pdb_style_entry)
        object.pandda_analyse_input_params_vbox.addLayout(pandda_pdb_style_hbox)

        pandda_mtz_style_hbox = QtGui.QHBoxLayout()
        label = QtGui.QLabel('mtz style')
        pandda_mtz_style_hbox.addWidget(label)
        object.pandda_mtz_style_entry = QtGui.QLineEdit()
        object.pandda_mtz_style_entry.setText('dimple.mtz')
        object.pandda_mtz_style_entry.setFixedWidth(200)
        pandda_mtz_style_hbox.addWidget(object.pandda_mtz_style_entry)
        object.pandda_analyse_input_params_vbox.addLayout(pandda_mtz_style_hbox)

        pandda_output_dir_hbox = QtGui.QHBoxLayout()
        label = QtGui.QLabel('output directory')
        object.pandda_analyse_input_params_vbox.addWidget(label)
        object.pandda_output_data_dir_entry = QtGui.QLineEdit()
        object.pandda_output_data_dir_entry.setText(object.panddas_directory)
        object.pandda_output_data_dir_entry.setFixedWidth(300)
        pandda_output_dir_hbox.addWidget(object.pandda_output_data_dir_entry)
        object.select_pandda_output_dir_button = QtGui.QPushButton("Select PANNDAs Directory")
        object.select_pandda_output_dir_button.setMaximumWidth(200)
        object.select_pandda_output_dir_button.clicked.connect(object.settings_button_clicked)
        pandda_output_dir_hbox.addWidget(object.select_pandda_output_dir_button)
        object.pandda_analyse_input_params_vbox.addLayout(pandda_output_dir_hbox)

        # qstat or local machine
        label = QtGui.QLabel('submit')
        object.pandda_analyse_input_params_vbox.addWidget(label)
        object.pandda_submission_mode_selection_combobox = QtGui.QComboBox()
        if object.external_software['qsub']:
            object.pandda_submission_mode_selection_combobox.addItem('qsub')
        object.pandda_submission_mode_selection_combobox.addItem('local machine')
        object.pandda_submission_mode_selection_combobox.setMaximumWidth(200)
        object.pandda_analyse_input_params_vbox.addWidget(object.pandda_submission_mode_selection_combobox)

        label = QtGui.QLabel('number of processors')
        object.pandda_analyse_input_params_vbox.addWidget(label)
        object.pandda_nproc = multiprocessing.cpu_count() - 1
        object.pandda_nproc_entry = QtGui.QLineEdit()
        object.pandda_nproc_entry.setText(str(object.pandda_nproc).replace(' ', ''))
        object.pandda_nproc_entry.setFixedWidth(200)
        object.pandda_analyse_input_params_vbox.addWidget(object.pandda_nproc_entry)

        label = QtGui.QLabel('order events by:')
        object.pandda_analyse_input_params_vbox.addWidget(label)
        object.pandda_sort_event_combobox = QtGui.QComboBox()
        object.pandda_sort_event_combobox.addItem('cluster_size')
        object.pandda_sort_event_combobox.addItem('z_peak')
        object.pandda_sort_event_combobox.setMaximumWidth(200)
        object.pandda_analyse_input_params_vbox.addWidget(object.pandda_sort_event_combobox)

        # crystal form option
        label = QtGui.QLabel('Use space group of reference file as filter:')
        object.pandda_analyse_input_params_vbox.addWidget(label)
        # reference file combobox, label with spg display
        hbox = QtGui.QHBoxLayout()
        object.reference_file_list = object.get_reference_file_list(' ')
        object.pandda_reference_file_selection_combobox = QtGui.QComboBox()
        object.populate_reference_combobox(object.pandda_reference_file_selection_combobox)
        object.pandda_reference_file_selection_combobox.activated[str].connect(object.change_pandda_spg_label)
        hbox.addWidget(object.pandda_reference_file_selection_combobox)
        object.pandda_reference_file_spg_label = QtGui.QLabel()
        hbox.addWidget(object.pandda_reference_file_spg_label)
        object.pandda_analyse_input_params_vbox.addLayout(hbox)

        label = QtGui.QLabel('\nExpert Parameters (only change if you know what you are doing!):')
        object.pandda_analyse_input_params_vbox.addWidget(label)

        object.wilson_checkbox = QtGui.QCheckBox('Wilson B-factor Scaling')
        self.layout_funcs.add_checkbox(object, object.wilson_checkbox, 'object.set_run_dimple_flag')

        # minimum number of datasets
        label = QtGui.QLabel('min_build_datasets')
        object.pandda_analyse_input_params_vbox.addWidget(label)
        object.pandda_min_build_dataset_entry = QtGui.QLineEdit()
        object.pandda_min_build_dataset_entry.setText('40')
        object.pandda_min_build_dataset_entry.setFixedWidth(200)
        object.pandda_analyse_input_params_vbox.addWidget(object.pandda_min_build_dataset_entry)

        label = QtGui.QLabel('max_new_datasets')
        object.pandda_analyse_input_params_vbox.addWidget(label)
        object.pandda_max_new_datasets_entry = QtGui.QLineEdit()
        object.pandda_max_new_datasets_entry.setText('500')
        object.pandda_max_new_datasets_entry.setFixedWidth(200)
        object.pandda_analyse_input_params_vbox.addWidget(object.pandda_max_new_datasets_entry)

        label = QtGui.QLabel(
            'grid_spacing (default=0.5)\nNote: higher values speed up calculations, but maps might be less pretty)')
        object.pandda_analyse_input_params_vbox.addWidget(label)
        object.pandda_grid_spacing_entry = QtGui.QLineEdit()
        object.pandda_grid_spacing_entry.setText('0.5')
        object.pandda_grid_spacing_entry.setFixedWidth(200)
        object.pandda_analyse_input_params_vbox.addWidget(object.pandda_grid_spacing_entry)

        object.pandda_analyse_input_params_vbox.addStretch(1)

        frame_right.setSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        frame_right.setLayout(object.pandda_analyse_input_params_vbox)

        grid_pandda.addWidget(frame_right, 0, 1, 5, 5)
        frame_pandda.setLayout(grid_pandda)
        object.pandda_analyse_hbox.addWidget(frame_pandda)

        # next three blocks display html documents created by pandda.analyse
        self.layout_funcs.pandda_html(object)

        object.pandda_initial_html = QtWebKit.QWebView()
        object.pandda_tab_dict['Dataset Summary'][1].addWidget(object.pandda_initial_html)
        object.pandda_initial_html.load(QtCore.QUrl(object.pandda_initial_html_file))
        object.pandda_initial_html.show()

        object.pandda_analyse_html = QtWebKit.QWebView()
        object.pandda_tab_dict['Processing Output'][1].addWidget(object.pandda_analyse_html)
        object.pandda_analyse_html.load(QtCore.QUrl(object.pandda_analyse_html_file))
        object.pandda_analyse_html.show()

        object.pandda_inspect_html = QtWebKit.QWebView()
        object.pandda_tab_dict['pandda.inspect'][1].addWidget(object.pandda_inspect_html)
        object.pandda_analyse_html.load(QtCore.QUrl(object.pandda_inspect_html_file))
        object.pandda_analyse_html.show()

        object.panddas_results_vbox.addWidget(object.pandda_tab_widget)
        object.show_pandda_html_summary()

        ################################################################################################################
        #                                                                                                              #
        #                                                REFINEMENT TAB                                                #
        #                                                                                                              #
        ################################################################################################################
        object.summary_vbox_for_table = QtGui.QVBoxLayout()

        # table
        object.refinement_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(object.refinement_table, object.refinement_table_columns)
        object.summary_vbox_for_table.addWidget(object.refinement_table)
        object.tab_dict[object.workflow_dict['Refinement']][1].addLayout(object.summary_vbox_for_table)

        ################################################################################################################
        #                                                                                                              #
        #                                                DEPOSITION TAB                                                #
        #                                                                                                              #
        ################################################################################################################
        object.deposition_vbox = QtGui.QVBoxLayout()

        scroll = QtGui.QScrollArea()
        object.deposition_vbox.addWidget(scroll)
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
            label_text = QtGui.QLabel('In your case: ' + object.html_export_directory)
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
        button.clicked.connect(object.export_to_html)
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
        button.clicked.connect(object.open_icm)
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
        label_text = QtGui.QLabel(XChemToolTips.zenodo_upload_start(object.html_export_directory))
        label_text.setStyleSheet("font: 17pt Arial")
        scrollLayout.addWidget(label_text)
        button = QtGui.QPushButton('prepare files')
        button.clicked.connect(object.prepare_files_for_zenodo_upload)
        button.setMaximumWidth(200)
        scrollLayout.addWidget(button)
        scrollLayout.addWidget(QtGui.QLabel('\n'))
        label_heading = QtGui.QLabel("5. ZENODO")
        label_heading.setStyleSheet("font: bold 20pt Arial")
        scrollLayout.addWidget(label_heading)
        label_text = QtGui.QLabel(XChemToolTips.zenodo_upload_part_one(object.html_export_directory))
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
        object.zenodo_upload_id_entry = QtGui.QLineEdit()
        object.zenodo_upload_id_entry.setFixedWidth(200)
        hbox_zenodo_upload_id.addWidget(object.zenodo_upload_id_entry)
        hbox_zenodo_upload_id.addStretch(1)
        scrollLayout.addLayout(hbox_zenodo_upload_id)
        button = QtGui.QPushButton('update html files with upload ID')
        button.clicked.connect(object.update_html_for_zenodo_upload)
        button.setMaximumWidth(300)
        scrollLayout.addWidget(button)
        scrollLayout.addWidget(QtGui.QLabel('\n'))
        label_heading = QtGui.QLabel("6. ZENODO upload HTML files")
        label_heading.setStyleSheet("font: bold 20pt Arial")
        scrollLayout.addWidget(label_heading)
        label_text = QtGui.QLabel(XChemToolTips.zenodo_upload_part_four(object.html_export_directory))
        label_text.setStyleSheet("font: 17pt Arial")
        scrollLayout.addWidget(label_text)
        scrollLayout.addWidget(QtGui.QLabel('\n'))
        scrollLayout.addStretch(1)
        scroll.setWidget(scrollContent)
        object.tab_dict[object.workflow_dict['Deposition']][1].addLayout(object.deposition_vbox)

        ################################################################################################################
        #                                                                                                              #
        #                                                 SETTINGS TAB                                                 #
        #                                                                                                              #
        ################################################################################################################
        object.settings_container = QtGui.QWidget()
        object.buttons_etc = QtGui.QWidget()
        object.settings_vbox = QtGui.QVBoxLayout()

        object.scroll = QtGui.QScrollArea(object.settings_container)
        object.settings_vbox.addWidget(object.scroll)
        scrollContent_settings = QtGui.QWidget(scroll)

        scrollLayout_settings = QtGui.QVBoxLayout(scrollContent_settings)
        scrollContent_settings.setLayout(scrollLayout_settings)

        # Settings Tab
        object.data_collection_vbox_for_settings = QtGui.QVBoxLayout()

        object.buttons_etc.setLayout(object.data_collection_vbox_for_settings)
        object.scroll.setWidget(object.buttons_etc)

        object.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nProject Directory: - REQUIRED -'))
        settings_hbox_initial_model_directory = QtGui.QHBoxLayout()
        object.initial_model_directory_label = QtGui.QLabel(object.initial_model_directory)
        settings_hbox_initial_model_directory.addWidget(object.initial_model_directory_label)
        settings_buttoon_initial_model_directory = QtGui.QPushButton('Select Project Directory')
        settings_buttoon_initial_model_directory.clicked.connect(object.settings_button_clicked)
        settings_hbox_initial_model_directory.addWidget(settings_buttoon_initial_model_directory)
        object.data_collection_vbox_for_settings.addLayout(settings_hbox_initial_model_directory)

        object.data_collection_vbox_for_settings.addWidget(
            QtGui.QLabel('\n\nReference Structure Directory: - OPTIONAL -'))
        settings_hbox_reference_directory = QtGui.QHBoxLayout()
        object.reference_directory_label = QtGui.QLabel(object.reference_directory)
        settings_hbox_reference_directory.addWidget(object.reference_directory_label)
        settings_buttoon_reference_directory = QtGui.QPushButton('Select Reference Structure Directory')
        settings_buttoon_reference_directory.clicked.connect(object.settings_button_clicked)
        settings_hbox_reference_directory.addWidget(settings_buttoon_reference_directory)
        object.data_collection_vbox_for_settings.addLayout(settings_hbox_reference_directory)

        object.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nData Source: - REQUIRED -'))
        settings_hbox_data_source_file = QtGui.QHBoxLayout()
        if object.data_source_file != '':
            object.data_source_file_label = QtGui.QLabel(
                os.path.join(object.database_directory, object.data_source_file))
        else:
            object.data_source_file_label = QtGui.QLabel('')
        settings_hbox_data_source_file.addWidget(object.data_source_file_label)
        settings_buttoon_data_source_file = QtGui.QPushButton('Select Data Source File')
        settings_buttoon_data_source_file.clicked.connect(object.settings_button_clicked)
        settings_hbox_data_source_file.addWidget(settings_buttoon_data_source_file)
        object.data_collection_vbox_for_settings.addLayout(settings_hbox_data_source_file)

        ## Data Collection
        object.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nData Collection Directory: - OPTIONAL -'))

        settings_beamline_frame = QtGui.QFrame()
        settings_beamline_frame.setFrameShape(QtGui.QFrame.StyledPanel)
        settings_beamline_vbox = QtGui.QVBoxLayout()

        settings_hbox_beamline_directory = QtGui.QHBoxLayout()
        object.beamline_directory_label = QtGui.QLabel(object.beamline_directory)
        settings_hbox_beamline_directory.addWidget(object.beamline_directory_label)
        settings_buttoon_beamline_directory = QtGui.QPushButton('Select Data Collection Directory')
        settings_buttoon_beamline_directory.clicked.connect(object.settings_button_clicked)
        settings_hbox_beamline_directory.addWidget(settings_buttoon_beamline_directory)
        settings_beamline_vbox.addLayout(settings_hbox_beamline_directory)

        settings_hbox_datasets_summary_file = QtGui.QHBoxLayout()
        object.datasets_summary_file_label = QtGui.QLabel(object.datasets_summary_file)
        settings_hbox_datasets_summary_file.addWidget(object.datasets_summary_file_label)
        settings_button_datasets_summary_file = QtGui.QPushButton('Select Existing\nCollection Summary File')
        settings_button_datasets_summary_file.clicked.connect(object.settings_button_clicked)
        settings_hbox_datasets_summary_file.addWidget(settings_button_datasets_summary_file)

        settings_button_new_datasets_summary_file = QtGui.QPushButton('Assign New\nCollection Summary File')
        settings_button_new_datasets_summary_file.clicked.connect(object.settings_button_clicked)
        settings_hbox_datasets_summary_file.addWidget(settings_button_new_datasets_summary_file)

        settings_beamline_vbox.addLayout(settings_hbox_datasets_summary_file)

        settings_beamline_frame.setLayout(settings_beamline_vbox)
        object.data_collection_vbox_for_settings.addWidget(settings_beamline_frame)

        object.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nCCP4_SCR Directory: - OPTIONAL -'))
        settings_hbox_ccp4_scratch_directory = QtGui.QHBoxLayout()
        object.ccp4_scratch_directory_label = QtGui.QLabel(object.ccp4_scratch_directory)
        settings_hbox_ccp4_scratch_directory.addWidget(object.ccp4_scratch_directory_label)
        settings_buttoon_ccp4_scratch_directory = QtGui.QPushButton('Select CCP4_SCR Directory')
        settings_buttoon_ccp4_scratch_directory.clicked.connect(object.settings_button_clicked)
        settings_hbox_ccp4_scratch_directory.addWidget(settings_buttoon_ccp4_scratch_directory)
        object.data_collection_vbox_for_settings.addLayout(settings_hbox_ccp4_scratch_directory)

        object.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nPANDDAs directory: - OPTIONAL -'))
        settings_hbox_panddas_directory = QtGui.QHBoxLayout()
        object.panddas_directory_label = QtGui.QLabel(object.panddas_directory)
        settings_hbox_panddas_directory.addWidget(object.panddas_directory_label)
        settings_button_panddas_directory = QtGui.QPushButton('Select PANNDAs Directory')
        settings_button_panddas_directory.clicked.connect(object.settings_button_clicked)
        settings_hbox_panddas_directory.addWidget(settings_button_panddas_directory)
        object.data_collection_vbox_for_settings.addLayout(settings_hbox_panddas_directory)

        object.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nHTML export directory: - OPTIONAL -'))
        settings_hbox_html_export_directory = QtGui.QHBoxLayout()
        object.html_export_directory_label = QtGui.QLabel(object.html_export_directory)
        settings_hbox_html_export_directory.addWidget(object.html_export_directory_label)
        settings_button_html_export_directory = QtGui.QPushButton('Select HTML Export Directory')
        settings_button_html_export_directory.clicked.connect(object.settings_button_clicked)
        settings_hbox_html_export_directory.addWidget(settings_button_html_export_directory)
        object.data_collection_vbox_for_settings.addLayout(settings_hbox_html_export_directory)

        object.data_collection_vbox_for_settings.addWidget(QtGui.QLabel('\n\nGroup deposition directory: - OPTIONAL -'))
        settings_hbox_group_deposition_directory = QtGui.QHBoxLayout()
        object.group_deposition_directory_label = QtGui.QLabel(object.group_deposit_directory)
        settings_hbox_group_deposition_directory.addWidget(object.group_deposition_directory_label)
        settings_button_group_deposition_directory = QtGui.QPushButton('Select Group deposition Directory')
        settings_button_group_deposition_directory.clicked.connect(object.settings_button_clicked)
        settings_hbox_group_deposition_directory.addWidget(settings_button_group_deposition_directory)
        object.data_collection_vbox_for_settings.addLayout(settings_hbox_group_deposition_directory)

        object.data_collection_vbox_for_settings.setSpacing(0)
        object.data_collection_vbox_for_settings.setContentsMargins(30, 0, 0, 0)

        object.buttons_etc.resize(object.screen.width() - 100, object.buttons_etc.sizeHint().height())
        object.tab_dict[object.workflow_dict['Settings']][1].addLayout(object.settings_vbox)

        ################################################################################################################
        #                                                                                                              #
        #                                                  STATUS BAR                                                  #
        #                                                                                                              #
        ################################################################################################################
        object.status_bar = QtGui.QStatusBar()
        object.progress_bar = QtGui.QProgressBar()
        object.progress_bar.setMaximum(100)
        object.status_bar.setMaximumWidth(object.screen.width())
        object.progress_bar.setMaximumWidth(object.screen.width())
        hbox_status = QtGui.QHBoxLayout()
        hbox_status.addWidget(object.status_bar)
        hbox_status.addWidget(object.progress_bar)

        vbox_main = QtGui.QVBoxLayout()
        menu_bar.setMaximumWidth(object.screen.width())
        vbox_main.addWidget(menu_bar)
        object.main_tab_widget.setMaximumSize(object.screen.width(), object.screen.height() - 245)
        vbox_main.addWidget(object.main_tab_widget)

        hboxTaskFrames = QtGui.QHBoxLayout()

        hboxTaskFrames.addWidget(update_from_datasource_button)
        hboxTaskFrames.addWidget(frame_dataset_task)
        hboxTaskFrames.addWidget(frame_map_cif_file_task)
        hboxTaskFrames.addWidget(frame_panddas_file_task)
        hboxTaskFrames.addWidget(frame_refine_file_task)

        vbox_main.addLayout(hboxTaskFrames)

        vbox_main.addLayout(hbox_status)

        object.window.setLayout(vbox_main)

        object.status_bar.showMessage('Ready')
        object.window.show()

        if object.data_source_file != '':
            write_enabled = object.check_write_permissions_of_data_source()
            if not write_enabled:
                object.data_source_set = False

    def workflow(self, object):
        ################################################################################################################
        #                                                                                                              #
        # ========================================== WORKFLOW TASK CONTAINER ========================================= #
        #                                                                                                              #
        ################################################################################################################
        object.workflow_widget_dict = {}

        # workflow task container - order of tabs as they appear for the main window
        object.workflow = ['Overview',  # 0
                           'Datasets',  # 1
                           'Maps',  # 2
                           'PANDDAs',  # 3
                           'Refinement',  # 4
                           'Deposition',  # 6
                           'Settings']  # 5

        # dictionary with keys corresponding to each stage in the workflow
        object.workflow_dict = {object.workflow[0]: 'Overview',
                                object.workflow[1]: 'Datasets',
                                object.workflow[2]: 'Maps',
                                object.workflow[3]: 'PANDDAs',
                                object.workflow[4]: 'Refinement',
                                object.workflow[6]: 'Settings',
                                object.workflow[5]: 'Deposition'}

        # tab widget
        object.main_tab_widget = QtGui.QTabWidget()
        object.tab_dict = {}
        self.layout_funcs.make_tab_dict(object.workflow, object.main_tab_widget, object.tab_dict)


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

    def add_to_box(self, frame, widgets_list):
        for widget in widgets_list:
            frame.addWidget(widget)
