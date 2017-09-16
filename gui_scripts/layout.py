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
    def __init__(self, xce_object):
        self.layout_funcs = LayoutFuncs()

    # function to initialise the top menu bar
    def initialise_menu_bar(self, xce_object):
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
        menu_bar = self.layout_funcs.setup_menubar(xce_object, menu_bar, self.menu_dict)

        # END OF MENU BAR - CODE BELOW: stuff removed from apo structure stuff that appears might have a funky
        # consequence - work out later.

        # xce_object.prepare_mmcif_files_dict = {}
        # xce_object.prepare_mmcif_files_dict['apo'] = prepare_mmcif_files_for_apo_structures
        # xce_object.prepare_mmcif_files_dict['ligand_bound'] = prepare_mmcif_files_for_ligand_bound_structures

        return menu_bar

    # function containing setup for bottom boxes
    def initialise_bottom_boxes(self, xce_object):

        # import all buttons
        setup().bottom_box_buttons(xce_object)

        # setup datasource button
        update_from_datasource_button = self.layout_funcs.setup_push_button(xce_object,
                                                                            xce_object.datasource_button_dict)

        ################################################################################################################
        #                                                                                                              #
        #                                               DATASETS BOX                                                   #
        #                                                                                                              #
        ################################################################################################################

        # setup the run button with push button function
        xce_object.dataset_task_run_button = self.layout_funcs.setup_push_button(xce_object,
                                                                                 xce_object.dataset_task_run_button_dict)

        # setup the task button with push button function
        xce_object.dataset_task_status_button = self.layout_funcs.setup_push_button(xce_object,
                                                                                    xce_object.dataset_task_status_button_dict)

        # array of both button xce_objects to apply to bottom box layout
        dataset_buttons = [xce_object.dataset_task_run_button, xce_object.dataset_task_status_button]

        # label for the bottom box layout
        dataset_label = "Datasets"

        # return the frame and combobox from the bottom box setup function
        frame_dataset_task, xce_object.dataset_tasks_combobox = self.layout_funcs.bottom_box_setup(xce_object,
                                                                                                   dataset_label,
                                                                                                   xce_object.dataset_tasks,
                                                                                                   'XChemToolTips.'
                                                                                                   'dataset_task_tip()',
                                                                                                   dataset_buttons,
                                                                                                   'background: '
                                                                                                   'rgb(240, 255, 140); ')

        # define the combobox and buttons in dictionary key to determine behaviour
        xce_object.workflow_widget_dict['Datasets'] = [xce_object.dataset_tasks_combobox,
                                                       xce_object.dataset_task_run_button,
                                                       xce_object.dataset_task_status_button]

        ################################################################################################################
        #                                                                                                              #
        #                                            MAPS & RESTRAINTS BOX                                             #
        #                                                                                                              #
        ################################################################################################################
        # settings for the run button

        # setup the run button with push button function
        xce_object.map_cif_file_task_run_button = \
            self.layout_funcs.setup_push_button(xce_object,
                                                xce_object.map_cif_file_task_run_button_dict)

        # setup the task button with push button function
        xce_object.map_cif_file_task_status_button = \
            self.layout_funcs.setup_push_button(xce_object,
                                                xce_object.map_cif_file_task_status_button_dict)

        # array of both button xce_objects to apply to bottom box layout
        map_cif_file_buttons = [xce_object.map_cif_file_task_run_button, xce_object.map_cif_file_task_status_button]

        # label for the bottom box layout
        map_cif_file_label = "Maps & Restraints"

        # return the frame and combobox from the bottom box setup function
        frame_map_cif_file_task, xce_object.map_cif_file_tasks_combobox = \
            self.layout_funcs.bottom_box_setup(xce_object,
                                               map_cif_file_label,
                                               xce_object.map_cif_file_tasks,
                                               'XChemToolTips.map_cif_file_'
                                               'task_tip()',
                                               map_cif_file_buttons,
                                               'background: rgb(140, 255, '
                                               '150); ')

        # define the combobox and buttons in dictionary key to determine behaviour
        xce_object.workflow_widget_dict['Maps'] = [xce_object.map_cif_file_tasks_combobox,
                                                   xce_object.map_cif_file_task_run_button,
                                                   xce_object.map_cif_file_task_status_button]

        ################################################################################################################
        #                                                                                                              #
        #                                               HIT IDENTIFICATION BOX                                         #
        #                                                                                                              #
        ################################################################################################################
        # settings for the run button

        # setup the run button with push button function
        xce_object.panddas_file_task_run_button = \
            self.layout_funcs.setup_push_button(xce_object,
                                                xce_object.panddas_file_task_run_button_dict)

        # setup the task button with push button function
        xce_object.panddas_file_task_status_button = \
            self.layout_funcs.setup_push_button(xce_object,
                                                xce_object.panddas_file_task_status_button_dict)

        # array of both button xce_objects to apply to bottom box layout
        panddas_file_buttons = [xce_object.panddas_file_task_run_button, xce_object.panddas_file_task_status_button]

        # label for the bottom box layout
        panddas_file_label = "Hit Identification"

        # return the frame and combobox from the bottom box setup function
        frame_panddas_file_task, xce_object.panddas_file_tasks_combobox = \
            self.layout_funcs.bottom_box_setup(xce_object,
                                               panddas_file_label,
                                               xce_object.panddas_file_tasks,
                                               'XChemToolTips.panddas_file_'
                                               'task_tip()',
                                               panddas_file_buttons,
                                               'background: rgb(140,200,255)'
                                               '; ')

        # define the combobox and buttons in dictionary key to determine behaviour
        xce_object.workflow_widget_dict['PANDDAs'] = [xce_object.panddas_file_tasks_combobox,
                                                      xce_object.panddas_file_task_run_button,
                                                      xce_object.panddas_file_task_status_button]

        ################################################################################################################
        #                                                                                                              #
        #                                                      REFINEMENT BOX                                          #
        #                                                                                                              #
        ################################################################################################################
        # settings for the run button

        # setup the run button with push button function
        xce_object.refine_file_task_run_button = \
            self.layout_funcs.setup_push_button(xce_object, xce_object.refine_file_task_run_button_dict)

        # setup the task button with push button function
        xce_object.refine_file_task_status_button = \
            self.layout_funcs.setup_push_button(xce_object, xce_object.refine_file_task_status_button_dict)

        # array of both button xce_objects to apply to bottom box layout
        refine_file_buttons = [xce_object.refine_file_task_run_button, xce_object.refine_file_task_status_button]

        # label for the bottom box layout
        refine_file_label = "Refinement"

        # return the frame and combobox from the bottom box setup function
        frame_refine_file_task, xce_object.refine_file_tasks_combobox = \
            self.layout_funcs.bottom_box_setup(xce_object,
                                               refine_file_label,
                                               xce_object.refine_file_tasks,
                                               'XChemToolTips.refine_file_task'
                                               '_tip()',
                                               refine_file_buttons,
                                               'background: rgb(245, 190, 255)'
                                               ';')

        # define the combobox and buttons in dictionary key to determine behaviour
        xce_object.workflow_widget_dict['Refinement'] = [xce_object.refine_file_tasks_combobox,
                                                         xce_object.refine_file_task_run_button,
                                                         xce_object.refine_file_task_status_button]

        return update_from_datasource_button, frame_dataset_task, frame_map_cif_file_task, frame_panddas_file_task, \
               frame_refine_file_task

    def main_layout(self, xce_object):
        # initialise menu bar
        menu_bar = self.initialise_menu_bar(xce_object)

        # initialise bottom boxes
        update_from_datasource_button, frame_dataset_task, frame_map_cif_file_task, frame_panddas_file_task, \
        frame_refine_file_task = self.initialise_bottom_boxes(xce_object)

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
        xce_object.overview_tab_widget = QtGui.QTabWidget()
        xce_object.overview_tab_dict = {}

        # make subtabs
        self.layout_funcs.make_tab_dict(overview_tab_list, xce_object.overview_tab_widget, xce_object.overview_tab_dict)

        # initiate the table in overview/datasource
        xce_object.overview_datasource_table = QtGui.QTableWidget()
        xce_object.overview_datasource_table.setSortingEnabled(True)
        xce_object.overview_datasource_table.resizeColumnsToContents()

        # initiate the graph in overview/summary
        xce_object.overview_figure, xce_object.overview_axes = plt.subplots()
        xce_object.overview_canvas = FigureCanvas(xce_object.overview_figure)
        xce_object.update_summary_plot()

        ################################################################################################################
        #                                                                                                              #
        #                                                 DATASETS TAB                                                 #
        #                                                                                                              #
        ################################################################################################################
        # define subtab list, widget and dict
        datasets_tab_list = ['Summary', 'Reprocess']
        xce_object.datasets_tab_widget = QtGui.QTabWidget()
        xce_object.datasets_tab_dict = {}

        # make subtabs
        self.layout_funcs.make_tab_dict(datasets_tab_list, xce_object.datasets_tab_widget, xce_object.datasets_tab_dict)

        # main body - things that are always displayed
        # add a container to hold everythting and add to main tab layout
        xce_object.datasets_data_collection_vbox = QtGui.QVBoxLayout()

        # add a horizontal box to hold option to autocheck for new data
        xce_object.autocheck_hbox = QtGui.QHBoxLayout()

        # checkbox for autocollect
        xce_object.check_for_new_data_collection = QtGui.QCheckBox('Check for new data collection every two minutes')
        self.layout_funcs.add_checkbox(xce_object, xce_object.check_for_new_data_collection,
                                       'xce_object.continously_check_for_new_data_collection')

        # select target dropdown
        select_target_label = QtGui.QLabel('Select Target: ')
        select_target_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        xce_object.target_selection_combobox = QtGui.QComboBox()
        xce_object.populate_target_selection_combobox(xce_object.target_selection_combobox)
        xce_object.target_selection_combobox.activated[str].connect(xce_object.target_selection_combobox_activated)
        xce_object.target = str(xce_object.target_selection_combobox.currentText())

        xce_object.autocheck_hbox_widgets = [xce_object.check_for_new_data_collection, select_target_label,
                                             xce_object.target_selection_combobox]  # array defining order of xce_objects to add

        self.layout_funcs.add_to_box(xce_object.autocheck_hbox,
                                     xce_object.autocheck_hbox_widgets)  # add xce_objects in order

        # add target dropdown to top bar
        xce_object.datasets_data_collection_vbox.addLayout(xce_object.autocheck_hbox)

        # summary sub-tab
        # table
        xce_object.datasets_summary_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(xce_object.datasets_summary_table, xce_object.datasets_summary_table_columns)
        xce_object.datasets_summarys_vbox_for_table = QtGui.QVBoxLayout()  # setup layout to hold table
        xce_object.datasets_summarys_vbox_for_table.addWidget(xce_object.datasets_summary_table)  # add table to layout
        xce_object.datasets_summarys_vbox_for_details = QtGui.QVBoxLayout()  # vbox for details
        xce_object.data_collection_details_currently_on_display = None  # switch for displaying/updating table

        xce_object.datasets_data_collection_vbox.addWidget(xce_object.datasets_tab_widget)  # add subtab to main tab

        # reprocessing sub-tab
        # top options
        # data collection label
        dc_label = QtGui.QLabel('Data collection directory: ')
        dc_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)  # align left and centre of container
        xce_object.diffraction_data_dir_label = QtGui.QLabel(
            xce_object.diffraction_data_directory)  # add directory as text
        xce_object.diffraction_data_dir_label.setAlignment(
            QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)  # align as above

        # select label
        select_button = QtGui.QPushButton("Select")
        select_button.clicked.connect(xce_object.select_diffraction_data_directory)  # attach file open dialogue

        # search button
        search_button = QtGui.QPushButton("Search Datasets")
        search_button.clicked.connect(xce_object.search_for_datasets)  # search for datasets in the selected directory

        # search info
        xce_object.diffraction_data_search_label = QtGui.QLabel(xce_object.diffraction_data_search_info)

        # translate label
        translate_label = QtGui.QLabel('translate: datasetID -> sampleID')
        translate_label.setAlignment(QtCore.Qt.AlignCenter)  # align in centre of container

        # CSV button
        csv_button = QtGui.QPushButton('Open CSV')
        csv_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
        csv_button.clicked.connect(xce_object.translate_datasetID_to_sampleID)  # open the relevant csv file

        # create hbox to hold everything and add widgets to it
        xce_object.hbox_select = QtGui.QHBoxLayout()  # top options box
        xce_object.hbox_select_widgets = [dc_label, select_button, search_button,
                                          xce_object.diffraction_data_search_label,
                                          translate_label,
                                          csv_button]  # array defining order of xce_objects to be added
        self.layout_funcs.add_to_box(xce_object.hbox_select, xce_object.hbox_select_widgets)  # add xce_objects in order

        # frame to hold everything
        frame_select = QtGui.QFrame()
        frame_select.setLayout(xce_object.hbox_select)  # apply to containing frame

        # table - main body
        xce_object.datasets_reprocess_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(xce_object.datasets_reprocess_table,
                                      xce_object.datasets_reprocess_columns)  # setup

        # create context menu - no idea where this lives...
        xce_object.popMenu_for_datasets_reprocess_table = QtGui.QMenu()
        run_xia2_on_selected = QtGui.QAction("mark selected for reprocessing", xce_object.window)
        run_xia2_on_selected.triggered.connect(xce_object.select_sample_for_xia2)
        xce_object.popMenu_for_datasets_reprocess_table.addAction(run_xia2_on_selected)
        xce_object.datasets_reprocess_table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        xce_object.datasets_reprocess_table.customContextMenuRequested.connect(
            xce_object.on_context_menu_reprocess_data)

        # options at bottom of tab
        # data processing label
        label = QtGui.QLabel('Data processing protocol: ')
        label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)

        # option checkboxes
        xce_object.xia2_3d_checkbox = QtGui.QCheckBox('xia2 3d')
        xce_object.xia2_3dii_checkbox = QtGui.QCheckBox('xia2 3dii')
        xce_object.xia2_dials_checkbox = QtGui.QCheckBox('Dials')

        # spacegroup label
        sg_label = QtGui.QLabel('Space Group:')

        # spacegroup dropdown menu
        xce_object.reprocess_space_group_comboxbox = QtGui.QComboBox()
        xce_object.reprocess_space_group_comboxbox.addItem('ignore')
        for sg in XChemMain.space_group_list():
            xce_object.reprocess_space_group_comboxbox.addItem(sg)

        # mtz label
        mtz_label = QtGui.QLabel('Reference MTZ:')

        # file label
        xce_object.reprocess_reference_mtz_file_label = QtGui.QLabel(xce_object.diffraction_data_reference_mtz)

        # select button
        select_button = QtGui.QPushButton("Select")
        select_button.clicked.connect(xce_object.select_reprocess_reference_mtz)

        # define order of widgets to be added to options hbox
        hbox_options_widgets = [label, xce_object.xia2_3d_checkbox, xce_object.xia2_3dii_checkbox,
                                xce_object.xia2_dials_checkbox, sg_label, xce_object.reprocess_space_group_comboxbox,
                                mtz_label, xce_object.reprocess_reference_mtz_file_label, select_button]

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
        xce_object.reprocess_isigma_combobox = QtGui.QComboBox()
        misigma = ['default', '4', '3', '2.5', '2', '1.5', '1', '0.5']
        self.layout_funcs.populate_combobox(misigma, xce_object.reprocess_isigma_combobox)
        xce_object.reprocess_isigma_combobox.setCurrentIndex(0)
        xce_object.reprocess_isigma_combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")

        # create vertical box to add labels and dropdowns to, create box and put in frame
        vbox_isigma = QtGui.QVBoxLayout()
        vbox_isigma_widgets = [label, xce_object.reprocess_isigma_combobox]
        self.layout_funcs.add_to_box(vbox_isigma, vbox_isigma_widgets)
        frame_isigma = QtGui.QFrame()
        frame_isigma.setLayout(vbox_isigma)

        # res limit cc half label
        res_cc_label = QtGui.QLabel('Resolution\nLimit:\nCC 1/2')
        label.setAlignment(QtCore.Qt.AlignCenter)

        # res limit cc half dropdown
        xce_object.reprocess_cc_half_combobox = QtGui.QComboBox()
        cc_half = ['default', '0.9', '0.8', '0.7', '0.6', '0.5', '0.4', '0.3', '0.2', '0.1']
        self.layout_funcs.populate_combobox(cc_half, xce_object.reprocess_cc_half_combobox)
        xce_object.reprocess_cc_half_combobox.setCurrentIndex(0)
        xce_object.reprocess_cc_half_combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")

        # create a vbox for label and dropdown, and add items to it
        vbox_cc_half = QtGui.QVBoxLayout()
        vbox_cc_half_widgets = [res_cc_label, xce_object.reprocess_cc_half_combobox]
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
        xce_object.reprocess_vbox = QtGui.QVBoxLayout()  # box to hold reprocessing subtab content
        xce_object.reprocess_hbox_widgets = [frame_select, xce_object.datasets_reprocess_table, bottom_options_frame]
        self.layout_funcs.add_to_box(xce_object.reprocess_vbox, xce_object.reprocess_hbox_widgets)

        ################################################################################################################
        #                                                                                                              #
        #                                                    MAPS TAB                                                  #
        #                                                                                                              #
        ################################################################################################################
        # select box for dimple
        xce_object.select_sample_for_dimple_box = QtGui.QCheckBox('(de-)select all samples for DIMPLE')
        self.layout_funcs.add_checkbox(xce_object, xce_object.select_sample_for_dimple_box,
                                       'xce_object.set_run_dimple_flag')

        # set new reference button
        set_new_reference_button = QtGui.QPushButton("Set New Reference (if applicable)")
        set_new_reference_button.clicked.connect(xce_object.set_new_reference_if_applicable)

        # refresh button
        refresh_reference_file_list_button = QtGui.QPushButton("Refresh reference file list")
        refresh_reference_file_list_button.clicked.connect(xce_object.refresh_reference_file_list)

        # list and populate reference files
        xce_object.reference_file_list = xce_object.get_reference_file_list(' ')
        xce_object.reference_file_selection_combobox = QtGui.QComboBox()
        xce_object.populate_reference_combobox(xce_object.reference_file_selection_combobox)

        # setup hbox to hold everything and add widgets
        xce_object.maps_checkbutton_hbox = QtGui.QHBoxLayout()
        maps_checkbutton_widgets = [xce_object.select_sample_for_dimple_box, set_new_reference_button,
                                    refresh_reference_file_list_button, xce_object.reference_file_selection_combobox]

        self.layout_funcs.add_to_box(xce_object.maps_checkbutton_hbox, maps_checkbutton_widgets)

        # table setup
        xce_object.maps_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(xce_object.maps_table, xce_object.maps_table_columns)

        # box for table, add to box, add to tab
        xce_object.initial_model_vbox_for_table = QtGui.QVBoxLayout()
        xce_object.initial_model_vbox_for_table.addWidget(xce_object.maps_table)

        # create context menu... no idea where this lives again.
        xce_object.popMenu_for_maps_table = QtGui.QMenu()
        run_dimple = QtGui.QAction("mark selected for dimple run", xce_object.window)
        run_dimple.triggered.connect(xce_object.select_sample_for_dimple)
        xce_object.popMenu_for_maps_table.addAction(run_dimple)
        xce_object.maps_table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        xce_object.maps_table.customContextMenuRequested.connect(xce_object.on_context_menu_initial_model)

        ################################################################################################################
        #                                                                                                              #
        #                                                   PANDDA TAB                                                 #
        #                                                                                                              #
        ################################################################################################################
        # list of subtabs in PanDDA tab
        pandda_tab_list = ['pandda.analyse',
                           'Dataset Summary',
                           'Processing Output',
                           'pandda.inspect',
                           'Statistical Map Summaries']

        # setup tab widget, set up tab dict, and make tab dict
        xce_object.pandda_tab_widget = QtGui.QTabWidget()
        xce_object.pandda_tab_dict = {}
        self.layout_funcs.make_tab_dict(pandda_tab_list, xce_object.pandda_tab_widget, xce_object.pandda_tab_dict)

        # pandda analyse subtab
        # setup a grid to hold everything
        grid_pandda = QtGui.QGridLayout()
        grid_pandda.setColumnStretch(0, 20)
        grid_pandda.setRowStretch(0, 20)

        # table - left
        xce_object.pandda_analyse_data_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(xce_object.pandda_analyse_data_table, xce_object.pandda_table_columns)

        # add table to grid
        frame_pandda = QtGui.QFrame()
        grid_pandda.addWidget(xce_object.pandda_analyse_data_table, 0, 0)

        # status of pandda job - under table
        xce_object.pandda_status = 'UNKNOWN'
        xce_object.pandda_status_label = QtGui.QLabel()

        # status options [filename, test to output, colour of text]
        pandda_status_options = [['/pandda.done', 'Finished!', 'color: green'],
                                 ['/pandda.running', 'Running...', 'color: orange'],
                                 ['/pandda.errored', 'Error encountered... please check the log files for pandda!',
                                  'color: red']]

        # enumerate text options and set text under table
        for option in pandda_status_options:
            if os.path.exists(str(xce_object.panddas_directory + option[0])):
                xce_object.pandda_status = option[1]
                xce_object.pandda_status_label.setStyleSheet(option[2])

        xce_object.pandda_status_label.setText(str('STATUS: ' + xce_object.pandda_status))
        xce_object.pandda_status_label.setFont(QtGui.QFont("Arial", 25, QtGui.QFont.Bold))
        grid_pandda.addWidget(xce_object.pandda_status_label, 3, 0)

        # input parameters for PANDDAs run - right
        frame_right = QtGui.QFrame()

        xce_object.pandda_analyse_input_params_vbox = QtGui.QVBoxLayout()

        # data directory section
        pandda_input_dir_hbox = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Data directory:')
        xce_object.pandda_analyse_input_params_vbox.addWidget(label)
        xce_object.pandda_input_data_dir_entry = QtGui.QLineEdit()
        xce_object.pandda_input_data_dir_entry.setText(os.path.join(xce_object.initial_model_directory, '*'))
        pandda_input_dir_hbox.addWidget(xce_object.pandda_input_data_dir_entry)
        xce_object.select_pandda_input_dir_button = QtGui.QPushButton("Select Input Template")
        xce_object.select_pandda_input_dir_button.clicked.connect(xce_object.select_pandda_input_template)
        pandda_input_dir_hbox.addWidget(xce_object.select_pandda_input_dir_button)
        xce_object.pandda_analyse_input_params_vbox.addLayout(pandda_input_dir_hbox)

        # pdb style section
        pandda_pdb_style_hbox = QtGui.QHBoxLayout()
        label = QtGui.QLabel('pdb style')
        pandda_pdb_style_hbox.addWidget(label)
        xce_object.pandda_pdb_style_entry = QtGui.QLineEdit()
        xce_object.pandda_pdb_style_entry.setText('dimple.pdb')
        pandda_pdb_style_hbox.addWidget(xce_object.pandda_pdb_style_entry)
        xce_object.pandda_analyse_input_params_vbox.addLayout(pandda_pdb_style_hbox)

        # mtz style section
        pandda_mtz_style_hbox = QtGui.QHBoxLayout()
        label = QtGui.QLabel('mtz style')
        pandda_mtz_style_hbox.addWidget(label)
        xce_object.pandda_mtz_style_entry = QtGui.QLineEdit()
        xce_object.pandda_mtz_style_entry.setText('dimple.mtz')
        pandda_mtz_style_hbox.addWidget(xce_object.pandda_mtz_style_entry)
        xce_object.pandda_analyse_input_params_vbox.addLayout(pandda_mtz_style_hbox)

        # spacer to separate out sections
        spacer = QtGui.QLabel(' ')
        xce_object.pandda_analyse_input_params_vbox.addWidget(spacer)

        # output directory section
        pandda_output_dir_hbox = QtGui.QHBoxLayout()
        label = QtGui.QLabel('Output directory:')
        xce_object.pandda_analyse_input_params_vbox.addWidget(label)
        xce_object.pandda_output_data_dir_entry = QtGui.QLineEdit()
        xce_object.pandda_output_data_dir_entry.setText(xce_object.panddas_directory)
        pandda_output_dir_hbox.addWidget(xce_object.pandda_output_data_dir_entry)
        xce_object.select_pandda_output_dir_button = QtGui.QPushButton("Select PANNDAs Directory")
        xce_object.select_pandda_output_dir_button.clicked.connect(xce_object.settings_button_clicked)
        pandda_output_dir_hbox.addWidget(xce_object.select_pandda_output_dir_button)
        xce_object.pandda_analyse_input_params_vbox.addLayout(pandda_output_dir_hbox)

        xce_object.pandda_analyse_input_params_vbox.addWidget(spacer)

        # qstat or local machine
        label = QtGui.QLabel('Submit:')
        xce_object.pandda_analyse_input_params_vbox.addWidget(label)
        xce_object.pandda_submission_mode_selection_combobox = QtGui.QComboBox()
        if xce_object.external_software['qsub']:
            xce_object.pandda_submission_mode_selection_combobox.addItem('qsub')
        xce_object.pandda_submission_mode_selection_combobox.addItem('local machine')
        xce_object.pandda_analyse_input_params_vbox.addWidget(xce_object.pandda_submission_mode_selection_combobox)

        # number of processors section
        label = QtGui.QLabel('Number of processors:')
        xce_object.pandda_analyse_input_params_vbox.addWidget(label)
        xce_object.pandda_nproc = multiprocessing.cpu_count() - 1
        xce_object.pandda_nproc_entry = QtGui.QLineEdit()
        xce_object.pandda_nproc_entry.setText(str(xce_object.pandda_nproc).replace(' ', ''))
        xce_object.pandda_analyse_input_params_vbox.addWidget(xce_object.pandda_nproc_entry)

        xce_object.pandda_analyse_input_params_vbox.addWidget(spacer)

        # how to order events
        label = QtGui.QLabel('Order events by:')
        xce_object.pandda_analyse_input_params_vbox.addWidget(label)
        xce_object.pandda_sort_event_combobox = QtGui.QComboBox()
        pandda_events = ['cluster_size', 'z_peak']
        self.layout_funcs.populate_combobox(pandda_events, xce_object.pandda_sort_event_combobox)
        xce_object.pandda_analyse_input_params_vbox.addWidget(xce_object.pandda_sort_event_combobox)

        # crystal form option
        label = QtGui.QLabel('Use space group of reference file as filter:')
        xce_object.pandda_analyse_input_params_vbox.addWidget(label)
        # reference file combobox, label with spg display
        hbox = QtGui.QHBoxLayout()
        xce_object.reference_file_list = xce_object.get_reference_file_list(' ')
        xce_object.pandda_reference_file_selection_combobox = QtGui.QComboBox()
        xce_object.populate_reference_combobox(xce_object.pandda_reference_file_selection_combobox)
        xce_object.pandda_reference_file_selection_combobox.activated[str].connect(xce_object.change_pandda_spg_label)
        hbox.addWidget(xce_object.pandda_reference_file_selection_combobox)
        xce_object.pandda_reference_file_spg_label = QtGui.QLabel()
        hbox.addWidget(xce_object.pandda_reference_file_spg_label)
        xce_object.pandda_analyse_input_params_vbox.addLayout(hbox)

        xce_object.pandda_analyse_input_params_vbox.addWidget(spacer)

        # Expert parameters
        label = QtGui.QLabel('Expert Parameters (only change if you know what you are doing!):')
        xce_object.pandda_analyse_input_params_vbox.addWidget(label)

        xce_object.pandda_analyse_input_params_vbox.addWidget(spacer)

        xce_object.wilson_checkbox = QtGui.QCheckBox('Wilson B-factor Scaling')
        self.layout_funcs.add_checkbox(xce_object, xce_object.wilson_checkbox, 'xce_object.set_run_dimple_flag')

        # minimum number of datasets
        label = QtGui.QLabel('min_build_datasets')
        xce_object.pandda_analyse_input_params_vbox.addWidget(label)
        xce_object.pandda_min_build_dataset_entry = QtGui.QLineEdit()
        xce_object.pandda_min_build_dataset_entry.setText('40')
        xce_object.pandda_analyse_input_params_vbox.addWidget(xce_object.pandda_min_build_dataset_entry)

        # maximum number of datasets
        label = QtGui.QLabel('max_new_datasets')
        xce_object.pandda_analyse_input_params_vbox.addWidget(label)
        xce_object.pandda_max_new_datasets_entry = QtGui.QLineEdit()
        xce_object.pandda_max_new_datasets_entry.setText('500')
        xce_object.pandda_analyse_input_params_vbox.addWidget(xce_object.pandda_max_new_datasets_entry)

        # grid spacing
        label = QtGui.QLabel(
            'grid_spacing (default=0.5)\nNote: higher values speed up calculations, but maps might be less pretty)')
        xce_object.pandda_analyse_input_params_vbox.addWidget(label)
        xce_object.pandda_grid_spacing_entry = QtGui.QLineEdit()
        xce_object.pandda_grid_spacing_entry.setText('0.5')
        xce_object.pandda_analyse_input_params_vbox.addWidget(xce_object.pandda_grid_spacing_entry)
        frame_right.setLayout(xce_object.pandda_analyse_input_params_vbox)

        grid_pandda.addWidget(frame_right, 0, 1, 5, 5)
        frame_pandda.setLayout(grid_pandda)

        # these are still currently populated in XCE.py - change
        xce_object.pandda_map_list = QtGui.QComboBox()
        xce_object.pandda_maps_html = QtWebKit.QWebView()

        # statistical map summaries vbox, add to vbox and add to layout
        xce_object.pandda_map_layout = QtGui.QVBoxLayout()
        pandda_map_layout_widgets = [xce_object.pandda_map_list, xce_object.pandda_maps_html]
        self.layout_funcs.add_to_box(xce_object.pandda_map_layout, pandda_map_layout_widgets)
        xce_object.pandda_maps_html.show()

        xce_object.pandda_analyse_hbox = QtGui.QHBoxLayout()
        xce_object.pandda_analyse_hbox.addWidget(frame_pandda)

        # next three blocks display html documents created by pandda.analyse
        self.layout_funcs.pandda_html(xce_object)

        xce_object.pandda_initial_html = QtWebKit.QWebView()
        xce_object.pandda_initial_html.load(QtCore.QUrl(xce_object.pandda_initial_html_file))
        xce_object.pandda_initial_html.show()

        xce_object.pandda_analyse_html = QtWebKit.QWebView()
        xce_object.pandda_analyse_html.load(QtCore.QUrl(xce_object.pandda_analyse_html_file))
        xce_object.pandda_analyse_html.show()

        xce_object.pandda_inspect_html = QtWebKit.QWebView()
        xce_object.pandda_analyse_html.load(QtCore.QUrl(xce_object.pandda_inspect_html_file))
        xce_object.pandda_analyse_html.show()

        xce_object.panddas_results_vbox = QtGui.QVBoxLayout()
        xce_object.panddas_results_vbox.addWidget(xce_object.pandda_tab_widget)
        xce_object.show_pandda_html_summary()

        ################################################################################################################
        #                                                                                                              #
        #                                                REFINEMENT TAB                                                #
        #                                                                                                              #
        ################################################################################################################
        xce_object.summary_vbox_for_table = QtGui.QVBoxLayout()

        # table
        xce_object.refinement_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(xce_object.refinement_table, xce_object.refinement_table_columns)
        xce_object.summary_vbox_for_table.addWidget(xce_object.refinement_table)

        ################################################################################################################
        #                                                                                                              #
        #                                                DEPOSITION TAB                                                #
        #                                                                                                              #
        ################################################################################################################
        xce_object.deposition_vbox = QtGui.QVBoxLayout()

        scroll = QtGui.QScrollArea()
        xce_object.deposition_vbox.addWidget(scroll)
        scrollContent = QtGui.QWidget(scroll)
        scrollLayout = QtGui.QVBoxLayout(scrollContent)
        scrollContent.setLayout(scrollLayout)

        # deposition page heading
        deposition_page_heading = self.layout_funcs.add_depo_heading('HTML export & ZENODO upload')
        deposition_page_heading.setStyleSheet("font: 30pt Arial Bold")

        # introduction text
        introduction_text = self.layout_funcs.add_depo_text(XChemToolTips.html_summary_introduction())

        # zenodo example image
        zenodo_image = QtGui.QLabel()
        zenodo_pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'), 'image', 'html_summary_page.png'))
        zenodo_image.setPixmap(zenodo_pixmap)

        # export html directory heading
        export_html_heading = self.layout_funcs.add_depo_heading('1. Specify HTML export directory in the settings tab')

        # export html directory text
        export_html_text = self.layout_funcs.add_depo_text(XChemToolTips.html_export_directory_background())

        # default for dls
        dls_dir_text = self.layout_funcs.add_depo_text(
            'Note: default for labxchem project at DLS is <labxchem_directory>/processing/html.')

        # user specific directory
        usr_dir_text = self.layout_funcs.add_depo_text('In your case: ' + xce_object.html_export_directory)

        # html prep heading
        html_prep_heading = self.layout_funcs.add_depo_heading("2. Prepare files for HTML export")

        # html prep test
        html_prep_text = self.layout_funcs.add_depo_text(XChemToolTips.html_export_step())

        # html export button
        html_export_button = QtGui.QPushButton('Export to HTML')
        html_export_button.clicked.connect(xce_object.export_to_html)
        html_export_button.setMaximumWidth(200)

        # prepare ICM heading
        ICM_heading = self.layout_funcs.add_depo_heading("3. Prepare ICM files")

        # ICM background text
        ICM_bgr_text = self.layout_funcs.add_depo_text(XChemToolTips.icb_file_background())

        # ICM prep text
        ICM_prep_text = self.layout_funcs.add_depo_text(XChemToolTips.prepare_ICB_files())

        # Open ICM button
        ICM_button = QtGui.QPushButton('Open ICM-pro')
        ICM_button.clicked.connect(xce_object.open_icm)
        ICM_button.setMaximumWidth(200)

        # ICM loading image
        ICM_image = QtGui.QLabel()
        ICM_pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'), 'image', 'drag_and_drop_icb_file.png'))
        ICM_image.setPixmap(ICM_pixmap)

        # ICM run image
        ICM_run_image = QtGui.QLabel()
        ICM_run_pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'), 'image', 'run_icm_script.png'))
        ICM_run_image.setPixmap(ICM_run_pixmap)

        # zenodo upload heading 1
        zenodo_upload_heading = self.layout_funcs.add_depo_heading("4. Prepare files for ZENODO upload")

        # zenodo upload text 1
        zenodo_upload_text = self.layout_funcs.add_depo_text(
            XChemToolTips.zenodo_upload_start(xce_object.html_export_directory))

        # prepare files button
        prep_files_button = QtGui.QPushButton('prepare files')
        prep_files_button.clicked.connect(xce_object.prepare_files_for_zenodo_upload)
        prep_files_button.setMaximumWidth(200)

        # zenodo upload heading 2
        zenodo_upload_heading2 = self.layout_funcs.add_depo_heading("5. ZENODO")

        # zenodo upload text 2
        zenodo_upload_text2 = self.layout_funcs.add_depo_text(
            XChemToolTips.zenodo_upload_part_one(xce_object.html_export_directory))

        # zenodo upload image
        zenodo_upload_image = QtGui.QLabel()
        zenodo_upload_pixmap = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'), 'image',
                                                          'new_zenodo_upload.png'))
        zenodo_upload_image.setPixmap(zenodo_upload_pixmap)

        # zenodo upload ID heading
        zenodo_upload_ID_heading = self.layout_funcs.add_depo_heading("5. ZENODO upload ID")

        # zenodo upload ID text
        zenodo_upload_ID_text = self.layout_funcs.add_depo_text(XChemToolTips.zenodo_upload_part_two())

        # zenodo upload ID image
        zenodo_upload_image2 = QtGui.QLabel()
        zenodo_upload_pixmap2 = QtGui.QPixmap(os.path.join(os.getenv('XChemExplorer_DIR'), 'image',
                                                           'zenodo_upload_id.png'))
        zenodo_upload_image2.setPixmap(zenodo_upload_pixmap2)

        # zenodo upload ID text 2
        zenodo_upload_ID_text2 = self.layout_funcs.add_depo_text(XChemToolTips.zenodo_upload_part_three())

        # zenodo upload ID entry box
        hbox_zenodo_upload_id = QtGui.QHBoxLayout()
        hbox_zenodo_upload_id.addWidget(QtGui.QLabel('upload ID:'))
        xce_object.zenodo_upload_id_entry = QtGui.QLineEdit()
        xce_object.zenodo_upload_id_entry.setFixedWidth(200)
        hbox_zenodo_upload_id.addWidget(xce_object.zenodo_upload_id_entry)
        hbox_zenodo_upload_id.addStretch(1)

        # update html button
        update_html_button = QtGui.QPushButton('update html files with upload ID')
        update_html_button.clicked.connect(xce_object.update_html_for_zenodo_upload)
        update_html_button.setMaximumWidth(300)

        # zenodo upload html files heading
        upload_html_heading = self.layout_funcs.add_depo_heading("6. ZENODO upload HTML files")

        # zenodo upload html text
        upload_html_text = self.layout_funcs.add_depo_text(XChemToolTips.zenodo_upload_part_four(xce_object.
                                                                                                 html_export_directory))

        deposition_widget_list = [deposition_page_heading, QtGui.QLabel('  '), introduction_text, QtGui.QLabel('  '),
                                  zenodo_image, QtGui.QLabel('  '), export_html_heading, export_html_text, dls_dir_text,
                                  usr_dir_text, QtGui.QLabel('  '), html_prep_heading, html_prep_text,
                                  html_export_button, QtGui.QLabel('  '), ICM_heading, ICM_bgr_text, ICM_prep_text,
                                  ICM_button, QtGui.QLabel('  '), ICM_image, QtGui.QLabel('  '), ICM_run_image,
                                  QtGui.QLabel('  '), zenodo_upload_heading, zenodo_upload_text, prep_files_button,
                                  QtGui.QLabel('  '), zenodo_upload_heading2, zenodo_upload_text2, zenodo_upload_image,
                                  QtGui.QLabel('  '), zenodo_upload_ID_heading, zenodo_upload_ID_text,
                                  zenodo_upload_image2, zenodo_upload_ID_text2]

        deposition_widget_list2 = [update_html_button, QtGui.QLabel('  '),
                                   upload_html_heading, upload_html_text, QtGui.QLabel('  ')]

        self.layout_funcs.add_to_box(scrollLayout, deposition_widget_list)
        scrollLayout.addLayout(hbox_zenodo_upload_id)
        self.layout_funcs.add_to_box(scrollLayout, deposition_widget_list2)

        # container settings
        scrollLayout.addStretch(1)
        scroll.setWidget(scrollContent)

        ################################################################################################################
        #                                                                                                              #
        #                                                 SETTINGS TAB                                                 #
        #                                                                                                              #
        ################################################################################################################
        xce_object.settings_container = QtGui.QWidget()
        xce_object.buttons_etc = QtGui.QWidget()
        xce_object.settings_vbox = QtGui.QVBoxLayout()

        xce_object.scroll = QtGui.QScrollArea(xce_object.settings_container)
        xce_object.settings_vbox.addWidget(xce_object.scroll)
        scrollContent_settings = QtGui.QWidget(xce_object.scroll)

        scrollLayout_settings = QtGui.QVBoxLayout(scrollContent_settings)
        scrollContent_settings.setLayout(scrollLayout_settings)

        # Settings Tab
        xce_object.data_collection_vbox_for_settings = QtGui.QVBoxLayout()

        xce_object.buttons_etc.setLayout(xce_object.data_collection_vbox_for_settings)
        xce_object.scroll.setWidget(xce_object.buttons_etc)

        xce_object.initial_model_directory_label = self.layout_funcs.settings_section_setup\
            (xce_object.data_collection_vbox_for_settings,
                                                 '\n\nProject Directory: - REQUIRED -',
                                                 xce_object.initial_model_directory,
                                                 'Select Project Directory',
                                                 xce_object.settings_button_clicked)

        xce_object.reference_directory_label = self.layout_funcs.settings_section_setup\
            (xce_object.data_collection_vbox_for_settings,
                                                 '\n\nReference Structure Directory: - OPTIONAL -',
                                                 xce_object.reference_directory,
                                                 'Select Reference Structure Directory',
                                                 xce_object.settings_button_clicked)

        if xce_object.data_source_file != '':
            xce_object.data_source_file_label_text = os.path.join(xce_object.database_directory,
                                                                  xce_object.data_source_file)
        else:
            xce_object.data_source_file_label_text = ''

        self.layout_funcs.settings_section_setup(xce_object.data_collection_vbox_for_settings,
                                                 '\n\nData Source: - REQUIRED -',
                                                 xce_object.data_source_file_label_text,
                                                 'Select Data Source File',
                                                 xce_object.settings_button_clicked)

        xce_object.data_collection_vbox_for_settings.addWidget(
            QtGui.QLabel('\n\nData Collection Directory: - OPTIONAL -'))

        settings_beamline_frame = QtGui.QFrame()
        settings_beamline_frame.setFrameShape(QtGui.QFrame.StyledPanel)
        settings_beamline_vbox = QtGui.QVBoxLayout()

        settings_hbox_beamline_directory = QtGui.QHBoxLayout()
        xce_object.beamline_directory_label = QtGui.QLabel(xce_object.beamline_directory)
        settings_hbox_beamline_directory.addWidget(xce_object.beamline_directory_label)
        settings_button_beamline_directory = QtGui.QPushButton('Select Data Collection Directory')
        settings_button_beamline_directory.setMaximumWidth(500)
        settings_button_beamline_directory.clicked.connect(xce_object.settings_button_clicked)
        settings_hbox_beamline_directory.addWidget(settings_button_beamline_directory)
        settings_beamline_vbox.addLayout(settings_hbox_beamline_directory)

        settings_hbox_datasets_summary_file = QtGui.QHBoxLayout()
        xce_object.datasets_summary_file_label = QtGui.QLabel(xce_object.datasets_summary_file)
        settings_hbox_datasets_summary_file.addWidget(xce_object.datasets_summary_file_label)
        settings_button_datasets_summary_file = QtGui.QPushButton('Select Existing\nCollection Summary File')
        settings_button_datasets_summary_file.setMaximumWidth(247)
        settings_button_datasets_summary_file.clicked.connect(xce_object.settings_button_clicked)
        settings_hbox_datasets_summary_file.addWidget(settings_button_datasets_summary_file)

        settings_button_new_datasets_summary_file = QtGui.QPushButton('Assign New\nCollection Summary File')
        settings_button_new_datasets_summary_file.clicked.connect(xce_object.settings_button_clicked)
        settings_button_new_datasets_summary_file.setMaximumWidth(247)
        settings_hbox_datasets_summary_file.addWidget(settings_button_new_datasets_summary_file)

        settings_beamline_vbox.addLayout(settings_hbox_datasets_summary_file)

        settings_beamline_frame.setLayout(settings_beamline_vbox)
        xce_object.data_collection_vbox_for_settings.addWidget(settings_beamline_frame)

        xce_object.ccp4_scratch_directory_label = self.layout_funcs.settings_section_setup\
            (xce_object.data_collection_vbox_for_settings,
                                                 '\n\nCCP4_SCR Directory: - OPTIONAL -',
                                                 xce_object.ccp4_scratch_directory,
                                                 'Select CCP4_SCR Directory',
                                                 xce_object.settings_button_clicked)

        xce_object.panddas_directory_label = \
            self.layout_funcs.settings_section_setup(xce_object.data_collection_vbox_for_settings,
                                                 '\n\nPANDDAs directory: - OPTIONAL -',
                                                 xce_object.panddas_directory,
                                                 'Select PANNDAs Directory',
                                                 xce_object.settings_button_clicked)

        xce_object.html_export_directory_label = self.layout_funcs.settings_section_setup\
            (xce_object.data_collection_vbox_for_settings,
                                                 '\n\nHTML export directory: - OPTIONAL -',
                                                 xce_object.html_export_directory,
                                                 'Select HTML Export Directory',
                                                 xce_object.settings_button_clicked)

        xce_object.group_deposition_directory_label = self.layout_funcs.settings_section_setup\
            (xce_object.data_collection_vbox_for_settings,
                                                 '\n\nGroup deposition directory: - OPTIONAL -',
                                                 xce_object.group_deposit_directory,
                                                 'Select Group deposition Directory',
                                                 xce_object.settings_button_clicked)

        xce_object.data_collection_vbox_for_settings.setSpacing(0)
        xce_object.data_collection_vbox_for_settings.setContentsMargins(30, 0, 0, 0)

        xce_object.buttons_etc.resize(xce_object.screen.width() - 100, xce_object.buttons_etc.sizeHint().height())

        ################################################################################################################
        #                                                                                                              #
        #                                                  STATUS BAR                                                  #
        #                                                                                                              #
        ################################################################################################################
        xce_object.status_bar = QtGui.QStatusBar()
        xce_object.progress_bar = QtGui.QProgressBar()
        xce_object.progress_bar.setMaximum(100)
        xce_object.status_bar.setMaximumWidth(xce_object.screen.width())
        xce_object.progress_bar.setMaximumWidth(xce_object.screen.width())
        hbox_status = QtGui.QHBoxLayout()
        hbox_status.addWidget(xce_object.status_bar)
        hbox_status.addWidget(xce_object.progress_bar)

        vbox_main = QtGui.QVBoxLayout()
        menu_bar.setMaximumWidth(xce_object.screen.width())
        vbox_main.addWidget(menu_bar)
        xce_object.main_tab_widget.setMaximumSize(xce_object.screen.width(), xce_object.screen.height() - 245)
        vbox_main.addWidget(xce_object.main_tab_widget)

        hboxTaskFrames = QtGui.QHBoxLayout()

        hboxTaskFrames.addWidget(update_from_datasource_button)
        hboxTaskFrames.addWidget(frame_dataset_task)
        hboxTaskFrames.addWidget(frame_map_cif_file_task)
        hboxTaskFrames.addWidget(frame_panddas_file_task)
        hboxTaskFrames.addWidget(frame_refine_file_task)

        vbox_main.addLayout(hboxTaskFrames)

        vbox_main.addLayout(hbox_status)

        xce_object.window.setLayout(vbox_main)

        xce_object.status_bar.showMessage('Ready')
        xce_object.window.show()

        if xce_object.data_source_file != '':
            write_enabled = xce_object.check_write_permissions_of_data_source()
            if not write_enabled:
                xce_object.data_source_set = False

    def workflow(self, xce_object):
        ################################################################################################################
        #                                                                                                              #
        # ========================================== WORKFLOW TASK CONTAINER ========================================= #
        #                                                                                                              #
        ################################################################################################################
        xce_object.workflow_widget_dict = {}

        # workflow task container - order of tabs as they appear for the main window
        xce_object.workflow = ['Overview',  # 0
                               'Datasets',  # 1
                               'Maps',  # 2
                               'PANDDAs',  # 3
                               'Refinement',  # 4
                               'Deposition',  # 6
                               'Settings']  # 5

        # dictionary with keys corresponding to each stage in the workflow
        xce_object.workflow_dict = {xce_object.workflow[0]: 'Overview',
                                    xce_object.workflow[1]: 'Datasets',
                                    xce_object.workflow[2]: 'Maps',
                                    xce_object.workflow[3]: 'PANDDAs',
                                    xce_object.workflow[4]: 'Refinement',
                                    xce_object.workflow[6]: 'Settings',
                                    xce_object.workflow[5]: 'Deposition'}

        # tab widget
        xce_object.main_tab_widget = QtGui.QTabWidget()
        xce_object.tab_dict = {}
        self.layout_funcs.make_tab_dict(xce_object.workflow, xce_object.main_tab_widget, xce_object.tab_dict)


class LayoutFuncs():
    def __init__(self):
        pass

    def make_tab_dict(self, tab_list, tab_widget, tab_dict):
        for page in tab_list:
            tab = QtGui.QWidget()
            vbox = QtGui.QVBoxLayout(tab)
            tab_widget.addTab(tab, page)
            tab_dict[page] = [tab, vbox]

    def add_checkbox(self, xce_object, checkbox, function, checkopt=False):
        checkbox.toggle()
        checkbox.setChecked(checkopt)
        eval(str('checkbox.stateChanged.connect(' + function + ')'))

    def table_setup(self, table, table_columns, sortingopt=True):
        table.setColumnCount(len(table_columns))
        table.setSortingEnabled(sortingopt)
        table.setHorizontalHeaderLabels(table_columns)
        table.resizeRowsToContents()
        table.resizeColumnsToContents()

    def pandda_html(self, xce_object):
        if os.path.exists(str(xce_object.panddas_directory + '/interesting_datasets')):
            print('WARNING: USING RESULTS FROM OLD PANDDA ANALYSE! THIS IS NOT FULLY SUPPORTED IN XCE2')
            print('PLEASE CHANGE YOUR PANDDA DIRECTORY TO A NEW RUN, OR USE THE OLD VERSION OF XCE!')
            xce_object.pandda_initial_html_file = str(
                xce_object.panddas_directory + '/results_summareis/pandda_initial.html')
            xce_object.pandda_analyse_html_file = str(
                xce_object.panddas_directory + '/results_summaries/pandda_analyse.html')
        xce_object.pandda_initial_html_file = str(
            xce_object.panddas_directory + '/analyses/html_summaries/' + 'pandda_initial.html')
        xce_object.pandda_analyse_html_file = str(
            xce_object.panddas_directory + '/analyses/html_summaries/' + 'pandda_analyse.html')
        xce_object.pandda_inspect_html_file = str(
            xce_object.panddas_directory + '/analyses/html_summaries/' + 'pandda_inspect.html')

    # function for datasource, run and status button setup
    def setup_push_button(self, xce_object, button_dict):
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
    def bottom_box_setup(self, xce_object, label, dropdown_options, dropdown_tooltip, buttons, colour):

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)
        frame.setStyleSheet("QFrame { border: 1px solid black; border-radius: 1px; padding: 0px; margin: 0px }")

        vbox = QtGui.QVBoxLayout()
        label = QtGui.QLabel(label)
        label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        label.setFont(xce_object.headlineLabelfont)
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
        frame.setMaximumWidth((xce_object.screen.width() - 20) / 5)

        return frame, combobox

    # function to add items to top menu bar
    def setup_menubar(self, xce_object, menu_bar, menu_items_dict):
        # use iterkeys to determine order of key by letter
        for config in sorted(menu_items_dict.iterkeys()):
            # add current item to menu bar
            menu = eval('menu_bar.addMenu("' + str(menu_items_dict[config][0]) + '")')
            # for each configuration item
            for menu_item in menu_items_dict[config][1]:
                # add the drop down option
                action = eval(str('QtGui.QAction("' + str(menu_item[0]) + '", xce_object.window)'))
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

    def populate_combobox(self, combobox_list, combobox):
        for item in combobox_list:
            combobox.addItem(item)

    def add_depo_heading(self, heading_text):
        heading = QtGui.QLabel(str(heading_text))
        heading.setStyleSheet("font: bold 20pt Arial")

        return heading

    def add_depo_text(self, text):
        out_text = QtGui.QLabel(text)
        out_text.setStyleSheet("font: 17pt Arial")

        return out_text

    def settings_section_setup(self, vbox, label_text, directory, button_text, button_function):
        vbox.addWidget(QtGui.QLabel(label_text))
        hbox = QtGui.QHBoxLayout()
        directory_label = QtGui.QLabel(directory)
        hbox.addWidget(directory_label)
        button = QtGui.QPushButton(button_text)
        button.setMaximumWidth(500)
        button.clicked.connect(button_function)
        hbox.addWidget(button)
        vbox.addLayout(hbox)

        return directory_label

    def add_widgets_layouts(self, xce_object):
        tab_add_widget = [
            [xce_object.tab_dict[xce_object.workflow_dict['Overview']][1], xce_object.overview_tab_widget],
            [xce_object.overview_tab_dict['Data Source'][1], xce_object.overview_datasource_table],
            [xce_object.overview_tab_dict['Summary'][1], xce_object.overview_canvas],
            [xce_object.pandda_tab_dict['Dataset Summary'][1], xce_object.pandda_initial_html],
            [xce_object.pandda_tab_dict['Processing Output'][1], xce_object.pandda_analyse_html],
            [xce_object.pandda_tab_dict['pandda.inspect'][1], xce_object.pandda_inspect_html]]

        tab_add_layout = [
            [xce_object.tab_dict[xce_object.workflow_dict['Datasets']][1], xce_object.datasets_data_collection_vbox],
            [xce_object.datasets_tab_dict['Summary'][1], xce_object.datasets_summarys_vbox_for_table],
            [xce_object.datasets_tab_dict['Summary'][1], xce_object.datasets_summarys_vbox_for_details],
            [xce_object.datasets_tab_dict['Reprocess'][1], xce_object.reprocess_vbox],
            [xce_object.tab_dict[xce_object.workflow_dict['Maps']][1], xce_object.maps_checkbutton_hbox],
            [xce_object.tab_dict[xce_object.workflow_dict['Maps']][1], xce_object.initial_model_vbox_for_table],
            [xce_object.pandda_tab_dict['Statistical Map Summaries'][1], xce_object.pandda_map_layout],
            [xce_object.pandda_tab_dict['pandda.analyse'][1], xce_object.pandda_analyse_hbox],
            [xce_object.tab_dict[xce_object.workflow_dict['PANDDAs']][1], xce_object.panddas_results_vbox],
            [xce_object.tab_dict[xce_object.workflow_dict['Refinement']][1], xce_object.summary_vbox_for_table],
            [xce_object.tab_dict[xce_object.workflow_dict['Deposition']][1], xce_object.deposition_vbox],
            [xce_object.tab_dict[xce_object.workflow_dict['Settings']][1], xce_object.settings_vbox]]

        for item in tab_add_widget:
            item[0].addWidget(item[1])

        for item in tab_add_layout:
            item[0].addLayout(item[1])
