import os, sys
from PyQt4 import QtCore, QtGui

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'lib'))
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'web'))

import XChemUtils
import XChemDB
import XChemMain
import XChemLog


class setup():
    def __init__(self):
        pass

    def set_xce_logfile(self, object):
        XChemLog.startLog(object.xce_logfile).create_logfile(object.xce_version)
        object.update_log = XChemLog.updateLog(object.xce_logfile)

    def settings(self, object):
        # set XCE version
        object.xce_version = 'v1.0'

        # general settings
        object.allowed_unitcell_difference_percent = 12
        object.acceptable_low_resolution_limit_for_data = 3.5
        object.filename_root = '${samplename}'
        object.data_source_set = False
        object.max_queue_jobs = 100

        ## directory settings

        # set current directory and direct log to it
        object.current_directory = os.getcwd()
        object.xce_logfile = os.path.join(object.current_directory, 'xce.log')

        # if in the correct place, set the various directories
        if 'labxchem' in object.current_directory:
            object.labxchem_directory = '/' + os.path.join(
                *object.current_directory.split('/')[1:6])  # need splat operator: *
            object.beamline_directory = os.path.join(object.labxchem_directory, 'processing', 'beamline')
            object.initial_model_directory = os.path.join(object.labxchem_directory, 'processing', 'analysis',
                                                          'initial_model')
            object.reference_directory = os.path.join(object.labxchem_directory, 'processing', 'reference')
            object.database_directory = os.path.join(object.labxchem_directory, 'processing', 'database')
            object.panddas_directory = os.path.join(object.labxchem_directory, 'processing', 'analysis', 'panddas')
            object.datasets_summary_file = os.path.join(object.database_directory,
                                                        str(os.getcwd().split('/')[5]) + '_summary.pkl')
            object.data_source_file = ''
            object.html_export_directory = os.path.join(object.labxchem_directory, 'processing', 'html')
            object.group_deposit_directory = os.path.join(object.labxchem_directory, 'processing', 'group_deposition')
            if os.path.isfile(
                    os.path.join(object.labxchem_directory, 'processing', 'database', 'soakDBDataFile.sqlite')):
                object.data_source_file = 'soakDBDataFile.sqlite'
                object.database_directory = os.path.join(object.labxchem_directory, 'processing', 'database')
                object.data_source_set = True
                object.db = XChemDB.data_source(os.path.join(object.database_directory, object.data_source_file))
                object.db.create_missing_columns()

            object.ccp4_scratch_directory = os.path.join(object.labxchem_directory, 'processing', 'tmp')

            directory_list = [object.beamline_directory, os.path.join(object.labxchem_directory,
                                                                      'processing', 'analysis'),
                              object.initial_model_directory, object.panddas_directory, object.reference_directory,
                              object.database_directory, object.ccp4_scratch_directory, object.html_export_directory,
                              object.group_deposit_directory]

            for directory in directory_list:
                if not os.path.isdir(directory):
                    os.mkdir(directory)

        # otherwise, use the current working directory
        else:
            object.beamline_directory = object.current_directory
            object.initial_model_directory = object.current_directory
            object.reference_directory = object.current_directory
            object.database_directory = object.current_directory
            object.data_source_file = ''
            object.ccp4_scratch_directory = os.getenv('CCP4_SCR')
            object.panddas_directory = object.current_directory
            object.datasets_summary_file = ''
            object.group_deposit_directory = object.current_directory

        ## deposition

        object.deposit_dict = {}

        ## internal lists and dictionaries

        object.data_collection_list = []
        object.visit_list = []
        object.target = ''
        object.dataset_outcome_combobox_dict = {}
        object.data_collection_dict = {}
        object.xtal_db_dict = {}
        object.pandda_analyse_input_table_dict = {}
        object.dewar_configuration_dict = {}
        object.data_collection_statistics_dict = {}
        object.initial_model_dimple_dict = {}  # contains toggle button if dimple should be run
        object.reference_file_list = []
        object.all_columns_in_data_source = XChemDB.data_source(os.path.join(object.database_directory,
                                                                             object.data_source_file)).return_column_list()
        object.albula_button_dict = {}  # using dials.image_viewer instead of albula, but keep name for dictionary
        object.xtalform_dict = {}

        object.dataset_outcome_dict = {}  # contains the dataset outcome buttons
        object.data_collection_table_dict = {}  # contains the dataset table
        object.data_collection_image_dict = {}
        object.data_collection_column_three_dict = {}
        object.datasets_summary_dict = {}
        object.diffraction_data_table_dict = {}
        object.refinement_table_dict = {}
        object.main_data_collection_table_exists = False
        object.timer_to_check_for_new_data_collection = QtCore.QTimer()

        object.target_list, object.visit_list = XChemMain.get_target_and_visit_list(object.beamline_directory)

        object.diffraction_data_dict = {}

        ## internal switches and flags

        object.explorer_active = 0
        object.coot_running = 0
        object.progress_bar_start = 0
        object.progress_bar_step = 0
        object.albula = None
        object.albula_subframes = []
        object.show_diffraction_image = None
        # can be any widget to be displayed in data collection summary tab
        object.data_collection_details_currently_on_display = None

        object.dataset_outcome = ["success",
                                  "Failed - centring failed",
                                  "Failed - no diffraction",
                                  "Failed - processing",
                                  "Failed - loop empty",
                                  "Failed - loop broken",
                                  "Failed - low resolution",
                                  "Failed - no X-rays",
                                  "Failed - unknown"]

        object.refinement_stage = ['0 - All Datasets',
                                   '1 - Analysis Pending',
                                   '2 - PANDDA model',
                                   '3 - In Refinement',
                                   '4 - CompChem ready',
                                   '5 - Deposition ready',
                                   '6 - Deposited']

        self.set_xce_logfile(object)

        ## external software packages
        object.using_remote_qsub_submission = False
        object.remote_qsub_submission = "ssh <dls fed ID>@nx.diamond.ac.uk 'module load global/cluster; qsub'"

        object.update_log = XChemLog.updateLog(object.xce_logfile)
        object.update_log.insert('new session started')
        object.diffraction_data_directory = object.current_directory
        object.diffraction_data_search_info = 'n/a'
        object.diffraction_data_reference_mtz = 'ignore'
        object.html_export_directory = os.getcwd()
        object.external_software = XChemUtils.external_software(object.xce_logfile).check()

        if object.external_software['acedrg']:
            object.restraints_program = 'acedrg'
            object.update_log.insert('will use ACEDRG for generation of ligand coordinates and restraints')
        elif object.external_software['phenix.elbow']:
            object.restraints_program = 'phenix.elbow'
            object.update_log.insert('will use PHENIX.ELBOW for generation of ligand coordinates and restraints')
        elif object.external_software['grade']:
            object.restraints_program = 'grade'
            object.update_log.insert('will use GRADE for generation of ligand coordinates and restraints')
        else:
            object.restraints_program = ''
            object.update_log.insert('No program for generation of ligand coordinates and restraints available!')

    def preferences(self, object):
        ## preferences

        object.preferences_data_to_copy = [['aimless logiles and merged mtz only', 'mtz_log_only'], ]
        # ['All Files in the respective auto-processing directory','everything'],

        object.preferences_selection_mechanism = ['IsigI*Comp*UniqueRefl',
                                                  'highest_resolution',
                                                  'lowest_Rfree']

        object.preferences = {'processed_data_to_copy': 'mtz_log_only',
                              'dataset_selection_mechanism': 'IsigI*Comp*UniqueRefl'}

        ## settings

        object.settings = {'current_directory': object.current_directory,
                           'beamline_directory': object.beamline_directory,
                           'datasets_summary': object.datasets_summary_file,
                           'initial_model_directory': object.initial_model_directory,
                           'panddas_directory': object.panddas_directory,
                           'reference_directory': object.reference_directory,
                           'database_directory': object.database_directory,
                           'data_source': os.path.join(object.database_directory, object.data_source_file),
                           'ccp4_scratch': object.ccp4_scratch_directory,
                           'unitcell_difference': object.allowed_unitcell_difference_percent,
                           'too_low_resolution_data': object.acceptable_low_resolution_limit_for_data,
                           'filename_root': object.filename_root,
                           'preferences': object.preferences,
                           'xce_logfile': object.xce_logfile,
                           'max_queue_jobs': object.max_queue_jobs,
                           'diffraction_data_directory': object.diffraction_data_directory,
                           'html_export_directory': object.html_export_directory,
                           'group_deposit_directory': object.group_deposit_directory,
                           'remote_qsub': ''}

    def tables(self, object):
        # Table column settings

        # functions that use tables.overview_datasource_table_columns:
        #
        # - select_datasource_columns_to_display() - dropdown in datasource top menu (select columns to show)
        # - populate_data_source_table()           - appears to be completely unused, so commented out
        # - populate_and_update_datasource_table() - used within select_datasource_columns_to_display and
        #                                            update_all_tables()

        object.overview_datasource_table_columns = ['Sample ID',
                                                    'Compound ID',
                                                    'Smiles',
                                                    'Visit',
                                                    'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                                    'Refinement\nRfree',
                                                    'Data Collection\nDate',
                                                    'Puck',
                                                    'PuckPosition',
                                                    'Ligand\nConfidence']

        # functions that use tables.datasets_summary_table_columns:
        #
        # - populate_datasets_summary_table() - appears in create_widgets_for_autoprocessing_results_only()
        # - user_update_selected_autoproc_datasets_summary_table()
        #                                            - appears in create_widgets_for_autoprocessing_results_only()

        object.datasets_summary_table_columns = ['Sample ID',
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
                                                 'Show Diffraction\nImage'
                                                 ]

        # functions that use tables.datasets_reprocess_columns:
        #
        # - update_datasets_reprocess_table() - appears in search_for_datasets()

        object.datasets_reprocess_columns = ['Dataset ID',
                                             'Sample ID',
                                             'Run\nxia2',
                                             'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                             'Rmerge\nLow',
                                             'Dimple\nRfree',
                                             'DataProcessing\nSpaceGroup',
                                             'DataProcessing\nUnitCell',
                                             'DataProcessing\nStatus']

        # functions that use tables.maps_table_columns:
        #
        # - update_datasets_reprocess_table() - appears in create_maps_table()

        object.maps_table_columns = ['Sample ID',
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

        # functions that use tables.pandda_table_columns:
        #
        # - populate_pandda_analyse_input_table() - appears in update_all_tables()

        object.pandda_table_columns = ['Sample ID',
                                       'Refinement\nSpace Group',
                                       'Resolution\n[Mn<I/sig(I)> = 1.5]',
                                       'Dimple\nRcryst',
                                       'Dimple\nRfree',
                                       'Crystal Form\nName', ]

        # functions that use tables.refinement_table_columns:
        #
        # - populate_and_update_refinement_table() - appears in update_all_tables

        object.refinement_table_columns = ['Sample ID',
                                           'Compound ID',
                                           'Refinement\nSpace Group',
                                           'Refinement\nResolution',
                                           'Refinement\nRcryst',
                                           'Refinement\nRfree',
                                           'Refinement\nOutcome',
                                           'PanDDA site details',
                                           'Refinement\nStatus']

    def top_menu_dict(self, object):
        # dictionary containing config for top menu setup
        # menu dict = { 'order letter: menu_item_name': ["name_in_menu",
        #                                       [
        #                                           ['text_in_menu', 'shortcut', 'trigger function']],
        #                                           [...],
        #                                       ]]
        #             }
        object.menu_dict = {'A: file': ["&File",
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
                                           ['Select columns to show', '',
                                            'object.select_datasource_columns_to_display'],
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
                                         'lambda: self.layout_funcs.openFile("/dls/science/groups/i04-1/software/docs/'
                                         'XChemExplorer.pdf")')],
                                     ['Troubleshooting', '', str(
                                         'lambda: self.layout_funcs.openFile("/dls/science/groups/i04-1/software/'
                                         'xce_troubleshooting.pdf")')]
                                 ]]

                     }

    def bottom_box_buttons(self, object):
        object.datasource_button_dict = {'datasource_button': [r"Update Tables\nFrom Datasource",
                                                        [
                                                            ['XChemToolTips.update_from_datasource_button_tip()',
                                                             # tooltip
                                                             'QPushButton { padding: 1px; margin: 1px; '
                                                             'background: rgb(140,140,140) }',
                                                             # stylesheet
                                                             'object.headlineLabelfont',  # font
                                                             'object.datasource_menu_reload_samples']  # action
                                                        ]]}

        object.dataset_task_run_button_dict = {'dataset_run_button':
                                            [r"Run",
                                             [
                                                 ['XChemToolTips.dataset_task_run_button_tip()',  # tooltip
                                                  'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                  '',  # font
                                                  'object.button_clicked']  # action
                                             ]]}

        object.dataset_task_status_button_dict = {'dataset_status_button':
                                               [r"Status",
                                                [
                                                    ['XChemToolTips.dataset_task_status_button_tip()',  # tooltip
                                                     'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                     '',  # font
                                                     'object.button_clicked']  # action
                                                ]]}

        object.map_cif_file_task_run_button_dict = {'dataset_run_button':
                                                 [r"Run",
                                                  [
                                                      ['XChemToolTips.map_cif_file_task_run_button_tip()',  # tooltip
                                                       'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                       '',  # font
                                                       'object.button_clicked']  # action
                                                  ]]}

        object.map_cif_file_task_status_button_dict = {'dataset_status_button':
                                                    [r"Status",
                                                     [
                                                         ['XChemToolTips.map_cif_file_task_status_button_tip()',
                                                          # tooltip
                                                          'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                          '',  # font
                                                          'object.button_clicked']  # action
                                                     ]]}

        object.panddas_file_task_run_button_dict = {'dataset_run_button':
                                                 [r"Run",
                                                  [
                                                      ['XChemToolTips.panddas_file_task_run_button_tip()',  # tooltip
                                                       'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                       '',  # font
                                                       'object.button_clicked']  # action
                                                  ]]}

        object.panddas_file_task_status_button_dict = {'dataset_status_button':
                                                    [r"Status",
                                                     [
                                                         ['XChemToolTips.panddas_file_task_status_button_tip()',
                                                          # tooltip
                                                          'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                          '',  # font
                                                          'object.button_clicked']  # action
                                                     ]]}

        object.refine_file_task_run_button_dict = {'dataset_run_button':
                                                [r"Run",
                                                 [
                                                     ['XChemToolTips.refine_file_task_run_button_tip()',  # tooltip
                                                      'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                      '',  # font
                                                      'object.button_clicked']  # action
                                                 ]]}

        object.refine_file_task_status_button_dict = {'dataset_status_button':
                                                   [r"Status",
                                                    [
                                                        ['XChemToolTips.refine_file_task_status_button_tip()',
                                                         # tooltip
                                                         'QPushButton { padding: 1px; margin: 1px }',  # stylesheet
                                                         '',  # font
                                                         'object.button_clicked']  # action
                                                    ]]}


