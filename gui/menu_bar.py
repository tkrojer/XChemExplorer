from import_modules import *

## operating system handling for opening external files
def openFile(file):
    if sys.platform == 'linux2':
        subprocess.call(["xdg-open", file])
    else:
        os.startfile(file)

def quit_xce(self):
    # save pkl file
    if self.data_collection_dict != {}:
        if os.path.isfile(self.data_collection_summary_file):
            self.update_log.insert('saving results to PKL file')
            pickle.dump(self.data_collection_dict,open(self.data_collection_summary_file,'wb'))
    self.update_log.insert('quitting XCE... bye,bye!')
    QtGui.qApp.quit()

def datasource_menu_save_samples(self):
    print 'hallo'

def datasource_menu_import_csv_file(self):
    if self.data_source_set:
        file_name = QtGui.QFileDialog.getOpenFileName(self.window,'Open file', self.database_directory)
        self.db.import_csv_file(file_name)
    else:
        self.update_status_bar('Please load a data source file first')

def datasource_menu_export_csv_file(self):
    file_name = str(QtGui.QFileDialog.getSaveFileName(self.window,'Save file', self.database_directory))
    if file_name.rfind('.') != -1:
        file_name=file_name[:file_name.rfind('.')]+'.csv'
    else:
        file_name=file_name+'.csv'
    self.db.export_to_csv_file(file_name)

def datasource_menu_update_datasource(self):
    self.work_thread=XChemThread.synchronise_db_and_filesystem(self.initial_model_directory,os.path.join(self.database_directory,self.data_source_file),self.panddas_directory,self.xce_logfile,'project_directory')
    self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
    self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
    self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
    self.connect(self.work_thread, QtCore.SIGNAL("datasource_menu_reload_samples"),self.datasource_menu_reload_samples)
    self.work_thread.start()

def select_datasource_columns_to_display(self):
    columns_to_show = QtGui.QMessageBox()
    columns_to_showLayout = columns_to_show.layout()
    columns_in_data_source=self.db.return_column_list()
    try:
        columns_in_data_source=self.db.return_column_list()
    except AttributeError:
        print '==> XCE: please select a datasource file'
        self.status_bar.showMessage('please select a datasource file')
        return

    column_dict={}
    vbox = QtGui.QVBoxLayout()
    number_of_entries=len(columns_in_data_source)
    columns_shown_in_dialog_column=15
    grid = QtGui.QGridLayout()
    x=0
    y=0

    columns_to_ignore=self.db.columns_not_to_display()
    for entries_added in range(number_of_entries):
        if not columns_in_data_source[entries_added][1] in columns_to_ignore:
            data_source_column = QtGui.QCheckBox(columns_in_data_source[entries_added][1])
            column_dict[entries_added]=data_source_column
            if columns_in_data_source[entries_added][1] in self.data_source_columns_to_display:
                data_source_column.setChecked(True)
            grid.addWidget(data_source_column, y,x)
            y+=1
        if y==columns_shown_in_dialog_column:
            y=0
            x+=1
    vbox.addLayout(grid)
    columns_to_showLayout.addLayout(vbox,0,0)

    columns_to_show.addButton(QtGui.QPushButton('OK'), QtGui.QMessageBox.YesRole)
    columns_to_show.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)
    reply=columns_to_show.exec_();
    if reply == 0:
        columns_to_show_list=['Sample ID']
        for key in column_dict:
            if column_dict[key].isChecked():
                columns_to_show_list.append(columns_in_data_source[key][1])
        self.data_source_columns_to_display=columns_to_show_list
        self.populate_and_update_data_source_table()

def create_new_data_source_func(self):
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
    self.database_directory=file_name[:file_name.rfind('/')]
    self.data_source_file=file_name[file_name.rfind('/')+1:]
    self.data_source_file_label.setText(os.path.join(self.database_directory,self.data_source_file))
    self.settings['database_directory']=self.database_directory
    self.settings['data_source']=self.data_source_file
    self.data_source_set=True
    self.datasource_menu_reload_samples()

def export_data_for_WONKA(self):
    self.update_log.insert('exporting CSV file for input into WONKA')
    self.db.export_csv_for_WONKA()

## menu bar options and processes
def menu_bar(self):
    # initiate a menubar instance
    self.menu_bar = QtGui.QMenuBar()

    ## File menu
    file = self.menu_bar.addMenu("&File")

    # open configuration file
    load = QtGui.QAction("Open Config File", self.window)
    load.setShortcut('Ctrl+O')

    # save configuration file
    load.triggered.connect(self.open_config_file)
    save = QtGui.QAction("Save Config File", self.window)
    save.setShortcut('Ctrl+S')
    save.triggered.connect(self.save_config_file)

    # quit xce
    quit = QtGui.QAction("Quit", self.window)
    quit.setShortcut('Ctrl+Q')
    quit.triggered.connect(quit_xce)

    # add the actions for the file menu
    file.addAction(load)
    file.addAction(save)
    file.addAction(quit)

    ## Datasource menu
    datasource_menu = self.menu_bar.addMenu("&Data Source")

    # reload samples from datasource
    reload_samples_from_datasource = QtGui.QAction('Reload Samples from Datasource', self.window)
    reload_samples_from_datasource.triggered.connect(self.datasource_menu_reload_samples)

    # save samples to datasource
    save_samples_to_datasource = QtGui.QAction('Save Samples to Datasource', self.window)
    save_samples_to_datasource.triggered.connect(datasource_menu_save_samples)

    # import csv
    import_csv_file_into_datasource = QtGui.QAction('Import CSV file into Datasource', self.window)
    import_csv_file_into_datasource.triggered.connect(datasource_menu_import_csv_file)

    # export csv
    export_csv_file_into_datasource = QtGui.QAction('Export CSV file from Datasource', self.window)
    export_csv_file_into_datasource.triggered.connect(datasource_menu_export_csv_file)

    # update datasource
    update_datasource = QtGui.QAction('Update Datasource from file system', self.window)
    update_datasource.triggered.connect(datasource_menu_update_datasource)

    # select columns to show
    select_columns_to_show = QtGui.QAction('Select columns to show', self.window)
    select_columns_to_show.triggered.connect(select_datasource_columns_to_display)

    # create new datasource
    create_new_data_source = QtGui.QAction('Create New Data Source (SQLite)', self.window)
    create_new_data_source.triggered.connect(create_new_data_source_func)

    # export csv (wonka)
    export_csv_for_WONKA = QtGui.QAction('export CSV for WONKA', self.window)
    export_csv_for_WONKA.triggered.connect(export_data_for_WONKA)

    # add the actions for the datasource menu
    datasource_menu.addAction(reload_samples_from_datasource)
    datasource_menu.addAction(save_samples_to_datasource)
    datasource_menu.addAction(import_csv_file_into_datasource)
    datasource_menu.addAction(export_csv_file_into_datasource)
    datasource_menu.addAction(update_datasource)
    datasource_menu.addAction(select_columns_to_show)
    datasource_menu.addAction(create_new_data_source)
    datasource_menu.addAction(export_csv_for_WONKA)

    ## Preferences menu
    preferences_menu = self.menu_bar.addMenu("&Preferences")

    # edit preferences
    show_preferences = QtGui.QAction('Edit Preferences', self.window)
    show_preferences.triggered.connect(self.show_preferences)

    # add the actions for the preferences menu
    preferences_menu.addAction(show_preferences)

    ## Deposition menu
    deposition_menu = self.menu_bar.addMenu("&Deposition")

    # edit information
    edit_deposition_info = QtGui.QAction('Edit Information', self.window)
    edit_deposition_info.triggered.connect(self.deposition_data)

    # export results to html
    export_results_to_html = QtGui.QAction('Export to HTML', self.window)
    export_results_to_html.triggered.connect(self.export_to_html)

    # find apo structures (PanDDA)
    find_apo_structures = QtGui.QAction('find PanDDA apo structures', self.window)
    find_apo_structures.triggered.connect(self.create_missing_apo_records_in_depositTable)

    # update apo records
    update_file_information_of_apo_records = QtGui.QAction('update file info of apo structures', self.window)
    update_file_information_of_apo_records.triggered.connect(self.update_file_information_of_apo_records)

    ## mmcif files
    # initiate a  dictionary object for mmcif preparation
    self.prepare_mmcif_files_dict = {}

    # prepare mmcif for apo structures
    prepare_mmcif_files_for_apo_structures = QtGui.QAction('prepare mmcif for apo structures', self.window)
    prepare_mmcif_files_for_apo_structures.triggered.connect(self.prepare_models_for_deposition)
    # add value to mmcif dict
    self.prepare_mmcif_files_dict['apo'] = prepare_mmcif_files_for_apo_structures

    # prepare mmcif for ligand bound structures
    prepare_mmcif_files_for_ligand_bound_structures = QtGui.QAction('prepare mmcif for ligand bound structures',
                                                                    self.window)
    prepare_mmcif_files_for_ligand_bound_structures.triggered.connect(self.prepare_models_for_deposition)
    # add value to mmcif dict
    self.prepare_mmcif_files_dict['ligand_bound'] = prepare_mmcif_files_for_ligand_bound_structures

    # prepare for group depo upload
    prepare_for_group_deposition_upload = QtGui.QAction('copy files to group deposition directory', self.window)
    prepare_for_group_deposition_upload.triggered.connect(self.prepare_for_group_deposition_upload)

    # enter pdb codes
    enter_pdb_codes = QtGui.QAction('Update DB with PDB codes', self.window)
    enter_pdb_codes.triggered.connect(self.enter_pdb_codes)

    # check smiles
    check_smiles = QtGui.QAction('Check SMILES', self.window)
    check_smiles.triggered.connect(self.check_smiles_in_db_and_pdb)


    # add actions to the deposition menu
    deposition_menu.addAction(edit_deposition_info)
    deposition_menu.addAction(export_results_to_html)
    deposition_menu.addAction(find_apo_structures)
    deposition_menu.addAction(update_file_information_of_apo_records)
    deposition_menu.addAction(prepare_mmcif_files_for_apo_structures)
    deposition_menu.addAction(prepare_mmcif_files_for_ligand_bound_structures)
    deposition_menu.addAction(prepare_for_group_deposition_upload)
    deposition_menu.addAction(enter_pdb_codes)
    deposition_menu.addAction(check_smiles)

    ## Help menu
    help_menu = self.menu_bar.addMenu("&Help")

    # load xce tutorial
    load_xce_tutorial = QtGui.QAction('Open XCE tutorial', self.window)
    file = '/dls/science/groups/i04-1/software/docs/XChemExplorer.pdf'
    load_xce_tutorial.triggered.connect(lambda: openFile(file))

    # load troubleshooting guide
    load_xce_troubleshoot = QtGui.QAction('Troubleshooting', self.window)
    file2 = '/dls/science/groups/i04-1/software/xce_troubleshooting.pdf'
    load_xce_troubleshoot.triggered.connect(lambda: openFile(file2))

    # add actions to the help menu
    help_menu.addAction(load_xce_tutorial)
    help_menu.addAction(load_xce_troubleshoot)

    # size constraints
    self.menu_bar.setMaximumWidth(self.screen.width())