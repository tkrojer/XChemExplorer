from import_modules import *

def settings_button_clicked(self):
    if self.sender().text()=='Select Project Directory':
        self.initial_model_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
        self.initial_model_directory_label.setText(self.initial_model_directory)
        self.settings['initial_model_directory']=self.initial_model_directory
    if self.sender().text()=='Select Reference Structure Directory':
        reference_directory_temp = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
        if reference_directory_temp != self.reference_directory:
            self.reference_directory=reference_directory_temp
            self.update_reference_files(' ')
        self.reference_directory_label.setText(self.reference_directory)
        self.settings['reference_directory']=self.reference_directory
    if self.sender().text()=='Select Data Source File':
        filepath_temp=QtGui.QFileDialog.getOpenFileNameAndFilter(self.window,'Select File', self.database_directory,'*.sqlite')
        filepath=str(tuple(filepath_temp)[0])
        self.data_source_file =   filepath.split('/')[-1]
        self.database_directory = filepath[:filepath.rfind('/')]
        self.settings['database_directory']=self.database_directory
        self.settings['data_source']=os.path.join(self.database_directory,self.data_source_file)
        write_enabled=self.check_write_permissions_of_data_source()
        if not write_enabled:
            self.data_source_set=False
        else:
            self.data_source_set=True
            self.data_source_file_label.setText(os.path.join(self.database_directory,self.data_source_file))
            self.db=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
            self.db.create_missing_columns()
            self.datasource_menu_reload_samples()
#                self.update_header_and_data_from_datasource()
#                self.populate_and_update_data_source_table()
#                self.create_initial_model_table()
    if self.sender().text()=='Select Data Collection Directory':
        dir_name = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
        if dir_name != self.beamline_directory:
            self.beamline_directory=dir_name
            self.target_list,self.visit_list=XChemMain.get_target_and_visit_list(self.beamline_directory)
#                self.target_list,self.visit_list=XChemMain.get_target_and_visit_list_for_Pietro(self.beamline_directory)
            self.populate_target_selection_combobox(self.target_selection_combobox)
        self.beamline_directory_label.setText(self.beamline_directory)
        self.settings['beamline_directory']=self.beamline_directory

    if self.sender().text()=='Select Existing\nCollection Summary File':
        if self.data_collection_summary_file != '':
            filepath_temp=QtGui.QFileDialog.getOpenFileNameAndFilter(self.window,'Select File', self.data_collection_summary_file[:self.data_collection_summary_file.rfind('/')],'*.pkl')
        else:
            filepath_temp=QtGui.QFileDialog.getOpenFileNameAndFilter(self.window,'Select File', os.getcwd(),'*.pkl')
        filepath=str(tuple(filepath_temp)[0])
        self.data_collection_summary_file=filepath
        self.data_collection_summary_file_label.setText(self.data_collection_summary_file)
        self.settings['data_collection_summary']=self.data_collection_summary_file

    if self.sender().text()=='Assign New\nCollection Summary File':
        if self.data_collection_summary_file != '':
            file_name = str(QtGui.QFileDialog.getSaveFileName(self.window,'New file', self.data_collection_summary_file[:self.data_collection_summary_file.rfind('/')]))
        else:
            file_name = str(QtGui.QFileDialog.getSaveFileName(self.window,'New file', self.current_directory))
        #make sure that the file always has .pkl extension
        if str(file_name).rfind('.') != -1:
            file_name=file_name[:file_name.rfind('.')]+'.pkl'
        else:
            file_name=file_name+'.pkl'
        self.data_collection_summary_file=file_name
        self.data_collection_summary_file_label.setText(self.data_collection_summary_file)
        self.settings['data_collection_summary']=self.data_collection_summary_file


    if self.sender().text()=='Select CCP4_SCR Directory':
        self.ccp4_scratch_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
        self.ccp4_scratch_directory_label.setText(self.ccp4_scratch_directory)
        self.settings['ccp4_scratch']=self.ccp4_scratch_directory
    if self.sender().text()=='Select PANNDAs Directory':
        self.panddas_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
        self.panddas_directory_label.setText(self.panddas_directory)
        self.pandda_output_data_dir_entry.setText(self.panddas_directory)
        print 'PANDDA',self.panddas_directory
        self.settings['panddas_directory']=self.panddas_directory
        self.pandda_initial_html_file=os.path.join(self.panddas_directory,'results_summaries','pandda_initial.html')
        self.pandda_analyse_html_file=os.path.join(self.panddas_directory,'results_summaries','pandda_analyse.html')
        self.pandda_inspect_html_file=os.path.join(self.panddas_directory,'results_summaries','pandda_inspect.html')

    if self.sender().text()=='Select HTML Export Directory':
        self.html_export_directory=str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
        self.html_export_directory_label.setText(self.html_export_directory)
        self.settings['html_export_directory']=self.html_export_directory

    if self.sender().text()=='Select Group deposition Directory':
        self.group_deposit_directory=str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
        self.group_deposition_directory_label.setText(self.group_deposit_directory)
        self.settings['group_deposit_directory']=self.group_deposit_directory
