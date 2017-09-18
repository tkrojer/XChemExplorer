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

from overview_tab import OverviewTab
from datasets_tab import DatasetsTab
from maps_tab import MapsTab
from pandda_tab import PanddaTab


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

        # Setup tabs
        OverviewTab().setup(xce_object)
        DatasetsTab().setup(xce_object)
        MapsTab().setup(xce_object)
        PanddaTab().setup(xce_object)





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
        zenodo_upload_ID_heading = self.layout_funcs.add_depo_heading("6. ZENODO upload ID")

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
        upload_html_heading = self.layout_funcs.add_depo_heading("7. ZENODO upload HTML files")

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

        xce_object.initial_model_directory_label = self.layout_funcs.settings_section_setup \
            (xce_object.data_collection_vbox_for_settings,
             '\n\nProject Directory: - REQUIRED -',
             xce_object.initial_model_directory,
             'Select Project Directory',
             xce_object.settings_button_clicked)

        xce_object.reference_directory_label = self.layout_funcs.settings_section_setup \
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

            self.data_source_file_label = self.layout_funcs.settings_section_setup \
                (xce_object.data_collection_vbox_for_settings,
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

        xce_object.ccp4_scratch_directory_label = self.layout_funcs.settings_section_setup \
            (xce_object.data_collection_vbox_for_settings,
             '\n\nCCP4_SCR Directory: - OPTIONAL -',
             xce_object.ccp4_scratch_directory,
             'Select CCP4_SCR Directory',
             xce_object.settings_button_clicked)

        xce_object.panddas_directory_label = self.layout_funcs.settings_section_setup \
            (xce_object.data_collection_vbox_for_settings,
             '\n\nPANDDAs directory: - OPTIONAL -',
             xce_object.panddas_directory,
             'Select PANNDAs Directory',
             xce_object.settings_button_clicked)

        xce_object.html_export_directory_label = self.layout_funcs.settings_section_setup \
            (xce_object.data_collection_vbox_for_settings,
             '\n\nHTML export directory: - OPTIONAL -',
             xce_object.html_export_directory,
             'Select HTML Export Directory',
             xce_object.settings_button_clicked)

        xce_object.group_deposition_directory_label = self.layout_funcs.settings_section_setup \
            (xce_object.data_collection_vbox_for_settings,
             '\n\nGroup deposition directory: - OPTIONAL -',
             xce_object.group_deposit_directory,
             'Select Group deposition Directory',
             xce_object.settings_button_clicked)

        #xce_object.data_collection_vbox_for_settings.setSpacing(0)
        xce_object.data_collection_vbox_for_settings.setContentsMargins(30, 30, 30, 30)

        xce_object.buttons_etc.resize(xce_object.buttons_etc.sizeHint().width() + 100, xce_object.buttons_etc.sizeHint()
                                      .height())

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
