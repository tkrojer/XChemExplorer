from import_modules import *

def settings_directories(self):
    self.xce_version = 'v1.0-beta.3.4'

    ## Misc. Settings
    self.allowed_unitcell_difference_percent = 12
    self.acceptable_low_resolution_limit_for_data = 3.5
    self.filename_root = '${samplename}'
    self.data_source_set = False
    self.max_queue_jobs = 100

    ## Directories
    # get current working directory and create log file
    self.current_directory = os.getcwd()
    self.xce_logfile = os.path.join(self.current_directory, 'xce.log')
    try:
        XChemLog.startLog(self.xce_logfile).create_logfile(self.xce_version)
    except IOError:
        self.xce_logfile = os.path.join(self.current_directory, 'xce_' + getpass.getuser() + '.log')
        XChemLog.startLog(self.xce_logfile).create_logfile(self.xce_version)
    self.update_log = XChemLog.updateLog(self.xce_logfile)
    self.update_log.insert('new session started')
    self.diffraction_data_directory = self.current_directory
    self.diffraction_data_search_info = 'n/a'
    self.diffraction_data_reference_mtz = 'ignore'
    self.html_export_directory = os.getcwd()

    # automatically define directory structure for labxchem locations
    if 'labxchem' in self.current_directory:
        self.labxchem_directory = '/' + os.path.join(*self.current_directory.split('/')[1:6])  # need splat operator: *
        self.beamline_directory = os.path.join(self.labxchem_directory, 'processing', 'beamline')
        self.initial_model_directory = os.path.join(self.labxchem_directory, 'processing', 'analysis', 'initial_model')
        self.reference_directory = os.path.join(self.labxchem_directory, 'processing', 'reference')
        self.database_directory = os.path.join(self.labxchem_directory, 'processing', 'database')
        self.panddas_directory = os.path.join(self.labxchem_directory, 'processing', 'analysis', 'panddas')
        self.data_collection_summary_file = os.path.join(self.database_directory,
                                                         str(os.getcwd().split('/')[5]) + '_summary.pkl')
        self.data_source_file = ''
        self.html_export_directory = os.path.join(self.labxchem_directory, 'processing', 'html')
        self.group_deposit_directory = os.path.join(self.labxchem_directory, 'processing', 'group_deposition')

        if os.path.isfile(os.path.join(self.labxchem_directory, 'processing', 'database', 'soakDBDataFile.sqlite')):
            self.data_source_file = 'soakDBDataFile.sqlite'
            self.database_directory = os.path.join(self.labxchem_directory, 'processing', 'database')
            self.data_source_set = True
            self.db = XChemDB.data_source(os.path.join(self.database_directory, self.data_source_file))
            self.db.create_missing_columns()

        self.ccp4_scratch_directory = os.path.join(self.labxchem_directory, 'processing', 'tmp')

        # set up missing directories in labxchem locations
        if not os.path.isdir(self.beamline_directory):
            os.mkdir(self.beamline_directory)
        if not os.path.isdir(os.path.join(self.labxchem_directory, 'processing', 'analysis')):
            os.mkdir(os.path.join(self.labxchem_directory, 'processing', 'analysis'))
        if not os.path.isdir(self.initial_model_directory):
            os.mkdir(self.initial_model_directory)
        if not os.path.isdir(self.panddas_directory):
            os.mkdir(self.panddas_directory)
        if not os.path.isdir(self.reference_directory):
            os.mkdir(self.reference_directory)
        if not os.path.isdir(self.database_directory):
            os.mkdir(self.database_directory)
        if not os.path.isdir(self.ccp4_scratch_directory):
            os.mkdir(self.ccp4_scratch_directory)
        if not os.path.isdir(self.html_export_directory):
            os.mkdir(self.html_export_directory)
        if not os.path.isdir(self.group_deposit_directory):
            os.mkdir(self.group_deposit_directory)

    # if not in labxchem, do nothing
    else:
        self.beamline_directory = self.current_directory
        self.initial_model_directory = self.current_directory
        self.reference_directory = self.current_directory
        self.database_directory = self.current_directory
        self.data_source_file = ''
        self.ccp4_scratch_directory = os.getenv('CCP4_SCR')
        self.panddas_directory = self.current_directory
        self.data_collection_summary_file = ''
        self.group_deposit_directory = self.current_directory

    ## Preferences
    # ['All Files in the respective auto-processing directory',           'everything'],
    self.preferences_data_to_copy = [['aimless logiles and merged mtz only', 'mtz_log_only'],]

    self.preferences_selection_mechanism = ['IsigI*Comp*UniqueRefl',
                                            'highest_resolution',
                                            'lowest_Rfree']

    self.preferences = {'processed_data_to_copy': 'mtz_log_only',
                        'dataset_selection_mechanism': 'IsigI*Comp*UniqueRefl'}

    ## Settings
    self.settings = {'current_directory': self.current_directory,
                     'beamline_directory': self.beamline_directory,
                     'data_collection_summary': self.data_collection_summary_file,
                     'initial_model_directory': self.initial_model_directory,
                     'panddas_directory': self.panddas_directory,
                     'reference_directory': self.reference_directory,
                     'database_directory': self.database_directory,
                     'data_source': os.path.join(self.database_directory, self.data_source_file),
                     'ccp4_scratch': self.ccp4_scratch_directory,
                     'unitcell_difference': self.allowed_unitcell_difference_percent,
                     'too_low_resolution_data': self.acceptable_low_resolution_limit_for_data,
                     'filename_root': self.filename_root,
                     'preferences': self.preferences,
                     'xce_logfile': self.xce_logfile,
                     'max_queue_jobs': self.max_queue_jobs,
                     'diffraction_data_directory': self.diffraction_data_directory,
                     'html_export_directory': self.html_export_directory,
                     'group_deposit_directory': self.group_deposit_directory}

    ## Deposition
    self.deposit_dict = {}

    ## internal lists and directories
    self.data_collection_list = []
    self.visit_list = []
    self.target = ''
    self.dataset_outcome_combobox_dict = {}
    self.data_collection_dict = {}
    self.xtal_db_dict = {}
    self.pandda_analyse_input_table_dict = {}
    self.dewar_configuration_dict = {}
    self.data_collection_statistics_dict = {}
    self.initial_model_dimple_dict = {}  # contains toggle button if dimple should be run
    self.reference_file_list = []
    self.all_columns_in_data_source = XChemDB.data_source(os.path.join(self.database_directory,
                                                                       self.data_source_file)).return_column_list()
    self.albula_button_dict = {}  # using dials.image_viewer instead of albula, but keep name for dictionary
    self.xtalform_dict = {}
    self.dataset_outcome_dict = {}  # contains the dataset outcome buttons
    self.data_collection_table_dict = {}  # contains the dataset table
    self.data_collection_image_dict = {}
    self.data_collection_column_three_dict = {}
    self.data_collection_summary_dict = {}
    self.diffraction_data_table_dict = {}
    self.summary_table_dict = {}
    self.main_data_collection_table_exists = False
    self.target_list, self.visit_list = XChemMain.get_target_and_visit_list(self.beamline_directory)
    self.diffraction_data_dict = {}

    ## internal switches and flags
    self.explorer_active = 0
    self.coot_running = 0
    self.progress_bar_start = 0
    self.progress_bar_step = 0
    self.albula = None
    self.albula_subframes = []
    self.show_diffraction_image = None
    self.data_collection_details_currently_on_display = None  # can be any widget to be displayed in data collection summary tab

    self.dataset_outcome = ["success",
                            "Failed - centring failed",
                            "Failed - no diffraction",
                            "Failed - processing",
                            "Failed - loop empty",
                            "Failed - loop broken",
                            "Failed - low resolution",
                            "Failed - no X-rays",
                            "Failed - unknown"]

    self.refinement_stage = ['0 - All Datasets',
                             '1 - Analysis Pending',
                             '2 - PANDDA model',
                             '3 - In Refinement',
                             '4 - CompChem ready',
                             '5 - Deposition ready',
                             '6 - Deposited']

    ## checking for external software packages
    self.external_software = external_software(self.xce_logfile).check()
    if self.external_software['acedrg']:
        self.restraints_program = 'acedrg'
        self.update_log.insert('will use ACEDRG for generation of ligand coordinates and restraints')
    elif self.external_software['phenix.elbow']:
        self.restraints_program = 'phenix.elbow'
        self.update_log.insert('will use PHENIX.ELBOW for generation of ligand coordinates and restraints')
    elif self.external_software['grade']:
        self.restraints_program = 'grade'
        self.update_log.insert('will use GRADE for generation of ligand coordinates and restraints')
    else:
        self.restraints_program = ''
        self.update_log.insert('No program for generation of ligand coordinates and restraints available!')