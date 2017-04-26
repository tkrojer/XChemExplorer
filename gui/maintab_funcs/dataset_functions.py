from import_modules import *

def autoprocessing_or_rescore(self, rescore_only):
    self.update_log.insert('checking for new data collection')
    start_thread = False
    if rescore_only:
        # first pop up a warning message as this will overwrite all user selections
        msgBox = QtGui.QMessageBox()
        msgBox.setText("*** WARNING ***\nThis will overwrite all your manual selections!\nDo you want to continue?")
        msgBox.addButton(QtGui.QPushButton('Yes'), QtGui.QMessageBox.YesRole)
        msgBox.addButton(QtGui.QPushButton('No'), QtGui.QMessageBox.RejectRole)
        reply = msgBox.exec_()
        if reply == 0:
            start_thread = True
        else:
            start_thread = False
    else:
        start_thread = True

def collection_autoupdate(self, state):
    self.timer = QtCore.QTimer()
    self.timer.timeout.connect(
        lambda: self.autoprocessing_or_rescore(False))
    if state == QtCore.Qt.Checked:
        print '==> XCE: checking automatically every 120s for new data collection'
        self.timer.start(120000)
    else:
        print '==> XCE: stopped checking for new data collections'
        self.timer.stop()

def target_selection_combobox_activated(self,text):
    self.target=str(text)

## WARNING: CURRENTLY UNUSED
def flag_sample_for_recollection(self):
    self.dewar_configuration_dict[self.dewar_label_active].setStyleSheet("background-color: yellow")

def undo_flag_sample_for_recollection(self):
    self.dewar_configuration_dict[self.dewar_label_active].setStyleSheet("background-color: gray")

def on_context_menu(self, point):
    # show context menu
    for key in self.dewar_configuration_dict:
        if self.dewar_configuration_dict[key]==self.sender():
            self.dewar_label_active=key
    self.popMenu.exec_(self.sender().mapToGlobal(point))

## END OF WARNING!

def select_diffraction_data_directory(self):
    self.diffraction_data_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
    self.diffraction_data_dir_label.setText(self.diffraction_data_directory)
    self.settings['diffraction_data_directory']=self.diffraction_data_directory
    self.update_log.insert('setting diffraction data directory to '+self.diffraction_data_directory)

def update_reprocess_datasets_table(self,data_dict):
    self.update_log.insert('updating reprocess datasets table')
    self.diffraction_data_table_dict={}
    self.diffraction_data_dict=data_dict

    self.diffraction_data_search_info='found '+str(len(self.diffraction_data_dict))+' datasets'
    self.diffraction_data_search_label.setText(self.diffraction_data_search_info)
    self.update_log.insert(self.diffraction_data_search_info)
    self.datasource_menu_reload_samples()

    # update table
    column_name=self.db.translate_xce_column_list_to_sqlite(self.reprocess_datasets_column_list)
    # set rows to 0
    self.reprocess_datasets_table.setRowCount(0)
    for entry in sorted(self.diffraction_data_dict):
        self.update_log.insert(str(self.diffraction_data_dict[entry]))
        if entry in self.xtal_db_dict:
            db_dict=self.xtal_db_dict[entry]
        else:
            db_dict={}
        row=self.reprocess_datasets_table.rowCount()
        self.reprocess_datasets_table.insertRow(row)
        for column,header in enumerate(column_name):
            if header[0]=='Dataset ID' or header[0]=='Sample ID':
                cell_text=QtGui.QTableWidgetItem()
                cell_text.setText(str(entry))
                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                self.reprocess_datasets_table.setItem(row, column, cell_text)
            elif header[0]=='Run\nxia2':
                run_xia2 = QtGui.QCheckBox()
                run_xia2.toggle()
                self.reprocess_datasets_table.setCellWidget(row, column, run_xia2)
                run_xia2.setChecked(False)
                self.diffraction_data_table_dict[entry]=[run_xia2]
            else:
                cell_text=QtGui.QTableWidgetItem()
                if db_dict != {}:
                    if header[0]=='DataProcessing\nStatus':
                        if str( db_dict[ header[1] ]  ) == 'running':
                            cell_text.setBackground(QtGui.QColor(100,230,150))
                        elif str( db_dict[ header[1] ]  ) == 'pending':
                            cell_text.setBackground(QtGui.QColor(20,100,230))
                        elif str( db_dict[ header[1] ]  ) == 'started':
                            cell_text.setBackground(QtGui.QColor(230,240,110))
                        elif str( db_dict[ header[1] ]  ) == 'finished':
                            cell_text.setBackground(QtGui.QColor(255,255,255))
                    cell_text.setText(str( db_dict[ header[1] ]  ))
                else:
                    cell_text.setText('')
                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                self.reprocess_datasets_table.setItem(row, column, cell_text)

def search_for_datasets(self):
    self.update_log.insert('search diffraction data directory for datasets...')
    self.work_thread=XChemMain.find_diffraction_image_directory_fast(self.diffraction_data_directory)
    self.explorer_active=1
    self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
    self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
    self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
    self.connect(self.work_thread, QtCore.SIGNAL("update_reprocess_datasets_table"), update_reprocess_datasets_table)
    self.work_thread.start()

def open_csv_file_translate_datasetID_to_sampleID(self):
    file_name_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(self.window,'Open file', self.current_directory,'*.csv')
    file_name=tuple(file_name_temp)[0]
    self.translate_datasetID_to_sampleID_csv_label.setText(file_name)
    self.translate_datasetID_to_sampleID_file=file_name

def translate_datasetID_to_sampleID(self):
    translate = QtGui.QMessageBox()
    translateLayout = translate.layout()
    self.translate_datasetID_to_sampleID_file='-'
    vbox = QtGui.QVBoxLayout()
    button=QtGui.QPushButton('Open CSV')
    button.clicked.connect(self.open_csv_file_translate_datasetID_to_sampleID)
    vbox.addWidget(button)
    self.translate_datasetID_to_sampleID_csv_label=QtGui.QLabel(self.translate_datasetID_to_sampleID_file)
    vbox.addWidget(self.translate_datasetID_to_sampleID_csv_label)
    translateLayout.addLayout(vbox,0,0)
    translate.addButton(QtGui.QPushButton('OK'), QtGui.QMessageBox.YesRole)
    translate.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)
    reply=translate.exec_()
    if reply == 0:
        if os.path.isfile(self.translate_datasetID_to_sampleID_file):
            trans_dict={}
            for line in open(self.translate_datasetID_to_sampleID_file):
                if len(line.split(','))==2:
                    dataset=line.split(',')[0]
                    new_sample_id=line.split(',')[1]
                    trans_dict[dataset]=new_sample_id
            if len(trans_dict) >= 1:
                allRows = self.reprocess_datasets_table.rowCount()
                for row in xrange(0,allRows):
                    dataset_id=str(self.reprocess_datasets_table.item(row,0).text())
                    sample_id=str(self.reprocess_datasets_table.item(row,1).text())
                    if dataset_id in trans_dict:
                        cell_text=QtGui.QTableWidgetItem()
                        cell_text.setText(trans_dict[dataset_id])
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.reprocess_datasets_table.setItem(row, 1, cell_text)
                        self.update_log.insert('dataset: %s -> changing sampleID to: %s' %(dataset_id,trans_dict[dataset_id]))

def on_context_menu_reprocess_data(self, point):
    # show context menu
    self.popMenu_for_reprocess_datasets_table.exec_(self.sender().mapToGlobal(point))

def select_reprocess_reference_mtz(self):
    self.update_log.insert('trying to set new reference mtz file for reprocessing with xia2')
    file_name = str(QtGui.QFileDialog.getOpenFileName(self.window,'Select file', self.database_directory))
    if os.path.isfile(file_name):
        if file_name.endswith('.mtz'):
            self.diffraction_data_reference_mtz=file_name
            self.update_log.insert('new reference file for data processing with xia2: '+self.diffraction_data_reference_mtz)
            self.reprocess_reference_mtz_file_label.setText(self.diffraction_data_reference_mtz)
        else:
            self.update_log.insert('this does not seem to be a mtz file: '+file_name)
