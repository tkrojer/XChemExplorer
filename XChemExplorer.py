########################################################################################################################
# DEVELOPER README:                                                                                                    #
# This is the main script, where the GUI is initialised from. All of the main layout objects live in their own scripts #
# under ./gui_scripts (i.e. the tab content). The settings and preferences script sets up all of the directory paths   #
# and contains dictionaries defining the top menu, push buttons and the tables held in the main tabs. The layout       #
# script contains functions for performing simple layout tasks, such as adding a combobox, and contains init.          #
# functions for all of the main layout functions.                                                                      #
#                                                                                                                      #
# In the future, the functions associated with buttons and frames etc. should be moved into the relevant script, but   #
# this is a bit more complicated. For now, they are separated out into sections within this script. The only GUI stuff #
# going on in here is calling the initialisation functions. To change the layout of a tab, edit it in it's own script, #
# and add any new functions in this script, in the relevant section. (If there is one yet)                             #
#                                                                                                                      #
# There's still a lot of cleaning up to be done in the future...                                                       #
########################################################################################################################

# solve gtk startup error
import gtk

gtk.set_interactive(False)

import base64
import getpass
import glob
import math
import multiprocessing
import pickle
import subprocess
import sys, os
import webbrowser
from datetime import datetime
from PyQt4 import QtGui, QtCore, QtWebKit

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'lib'))
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'web'))
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'gui_scripts'))

from settings_preferences import *
from layout import *
from stylesheet import set_stylesheet
from proasis_functions import Proasis

from messagebox import DepositionDetails
from messagebox import CheckAutoProcessing


from XChemUtils import parse
import XChemThread
import XChemDB
import XChemPANDDA
import XChemToolTips
import XChemMain
import XChemPlots
import XChemLog
import XChemProcess
import XChemDeposit
import XChemWeb

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

class XChemExplorer(QtGui.QApplication):
    def __init__(self, args):

        # init a QApplication object to hold XCE
        QtGui.QApplication.__init__(self, args)

        # start GUI
        self.start_GUI()

        # set stylesheet - how the gui looks
        set_stylesheet(self)

        self.exec_()

    def start_GUI(self):

        # check http://doc.qt.io/qt-4.8/stylesheet-customizing.html#the-box-model
        # This needs moving somewhere more appropriate...
        self.headlineLabelfont = QtGui.QFont("Arial", 20, QtGui.QFont.Bold)

        setup().settings(self)
        setup().preferences(self)
        setup().tables(self)
        Proasis().proasis_menu(self)

        self.layout_funcs = LayoutFuncs()

        # GUI setup
        self.window = QtGui.QWidget()
        self.window.setWindowTitle("XChemExplorer")
        self.screen = QtGui.QDesktopWidget().screenGeometry()

        LayoutObjects(self).workflow(self)
        LayoutObjects(self).generic_message_box(self)
        LayoutObjects(self).main_layout(self)
        LayoutFuncs().add_widgets_layouts(self)



    # function to update datasource
    def datasource_menu_reload_samples(self):
        self.update_log.insert(
            'reading samples from data source: ' + os.path.join(self.database_directory, self.data_source_file))
        self.update_status_bar(
            'reading samples from data source: ' + os.path.join(self.database_directory, self.data_source_file))
        self.update_header_and_data_from_datasource()
        self.update_all_tables()
        self.overview_datasource_table.resizeColumnsToContents()

    # function to create new datasource
    def create_new_data_source(self):
        file_name = str(QtGui.QFileDialog.getSaveFileName(self.window, 'Save file', self.database_directory))
        # make sure that the file always has .sqlite extension
        if file_name.rfind('.') != -1:
            file_name = file_name[:file_name.rfind('.')] + '.sqlite'
        else:
            file_name = file_name + '.sqlite'
        self.db = XChemDB.data_source(file_name)
        print('==> XCE: creating new data source')
        self.db.create_empty_data_source_file()
        self.db.create_missing_columns()
        self.database_directory = file_name[:file_name.rfind('/')]
        self.data_source_file = file_name[file_name.rfind('/') + 1:]
        self.data_source_file_label.setText(os.path.join(self.database_directory, self.data_source_file))
        self.settings['database_directory'] = self.database_directory
        self.settings['data_source'] = self.data_source_file
        self.data_source_set = True
        self.datasource_menu_reload_samples()


    ####################################################################################################################
    #                                                                                                                  #
    #                                                 DATASETS TAB                                                     #
    #                                                                                                                  #
    ####################################################################################################################
    def continously_check_for_new_data_collection(self, state):
        self.timer_to_check_for_new_data_collection.timeout.connect(
            lambda: self.check_for_new_autoprocessing_or_rescore(False))
        if state == QtCore.Qt.Checked:
            print('==> XCE: checking automatically every 120s for new data collection')
            self.timer_to_check_for_new_data_collection.start(120000)
        else:
            print('==> XCE: stopped checking for new data collections')
            self.timer_to_check_for_new_data_collection.stop()

    def target_selection_combobox_activated(self, text):
        self.target = str(text)

    def select_diffraction_data_directory(self):
        self.diffraction_data_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
        self.diffraction_data_dir_label.setText(self.diffraction_data_directory)
        self.settings['diffraction_data_directory'] = self.diffraction_data_directory
        self.update_log.insert('setting diffraction data directory to ' + self.diffraction_data_directory)

    def search_for_datasets(self):
        self.update_log.insert('search diffraction data directory for datasets...')
        print('will search ' + str(self.diffraction_data_directory))
        self.work_thread = XChemMain.find_diffraction_image_directory_fast(self.diffraction_data_directory)
        self.explorer_active = 1

        self.connect(self.work_thread, QtCore.SIGNAL("update_datasets_reprocess_table"),
                     self.update_datasets_reprocess_table)
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)

        self.work_thread.start()

        #self.work_thread = self.update_datasets_reprocess_table(self.diffraction_data_directory)

    def translate_datasetID_to_sampleID(self):
        translate = QtGui.QMessageBox()
        translateLayout = translate.layout()
        self.translate_datasetID_to_sampleID_file = '-'
        vbox = QtGui.QVBoxLayout()
        button = QtGui.QPushButton('Open CSV')
        button.clicked.connect(self.open_csv_file_translate_datasetID_to_sampleID)
        vbox.addWidget(button)
        self.translate_datasetID_to_sampleID_csv_label = QtGui.QLabel(self.translate_datasetID_to_sampleID_file)
        vbox.addWidget(self.translate_datasetID_to_sampleID_csv_label)
        translateLayout.addLayout(vbox, 0, 0)
        translate.addButton(QtGui.QPushButton('OK'), QtGui.QMessageBox.YesRole)
        translate.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)
        reply = translate.exec_();
        if reply == 0:
            if os.path.isfile(self.translate_datasetID_to_sampleID_file):
                trans_dict = {}
                for line in open(self.translate_datasetID_to_sampleID_file):
                    if len(line.split(',')) == 2:
                        dataset = line.split(',')[0]
                        new_sample_id = line.split(',')[1]
                        trans_dict[dataset] = new_sample_id
                if len(trans_dict) >= 1:
                    allRows = self.datasets_reprocess_table.rowCount()
                    for row in xrange(0, allRows):
                        dataset_id = str(self.datasets_reprocess_table.item(row, 0).text())
                        sample_id = str(self.datasets_reprocess_table.item(row, 1).text())
                        if dataset_id in trans_dict:
                            cell_text = QtGui.QTableWidgetItem()
                            cell_text.setText(trans_dict[dataset_id])
                            cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                            self.datasets_reprocess_table.setItem(row, 1, cell_text)
                            self.update_log.insert(
                                'dataset: {0!s} -> changing sampleID to: {1!s}'.format(dataset_id,
                                                                                       trans_dict[dataset_id]))

    def select_sample_for_xia2(self):
        indexes = self.datasets_reprocess_table.selectionModel().selectedRows()
        for index in sorted(indexes):
            xtal = str(self.datasets_reprocess_table.item(index.row(), 1).text())
            print(xtal, self.diffraction_data_table_dict[xtal][0])
            self.update_log.insert('{0!s} marked for reprocessing'.format(index.row()))
            self.diffraction_data_table_dict[xtal][0].setChecked(True)

    def select_reprocess_reference_mtz(self):
        self.update_log.insert('trying to set new reference mtz file for reprocessing with xia2')
        file_name = str(QtGui.QFileDialog.getOpenFileName(self.window, 'Select file', self.database_directory))
        if os.path.isfile(file_name):
            if file_name.endswith('.mtz'):
                self.diffraction_data_reference_mtz = file_name
                self.update_log.insert(
                    'new reference file for data processing with xia2: ' + self.diffraction_data_reference_mtz)
                self.reprocess_reference_mtz_file_label.setText(self.diffraction_data_reference_mtz)
            else:
                self.update_log.insert('this does not seem to be a mtz file: ' + file_name)


    def select_best_autoprocessing_result(self):
        start_thread = True
        if self.rescore:
            # first pop up a warning message as this will overwrite all user selections
            msgBox = QtGui.QMessageBox()
            msgBox.setText("*** WARNING ***\nThis will overwrite all your manual selections!\nDo you want to continue?")
            msgBox.addButton(QtGui.QPushButton('Yes'), QtGui.QMessageBox.YesRole)
            msgBox.addButton(QtGui.QPushButton('No'), QtGui.QMessageBox.RejectRole)
            reply = msgBox.exec_();
            if reply != 0:
                start_thread = False

        if start_thread:
            self.update_log.insert('selecting best autoprocessing result')
            self.update_log.insert('samples where user made manual changes will be ignored!')
            rescore = False
            processedDir=os.path.join(self.beamline_directory,'processed',self.target)
            visit,beamline = XChemMain.getVisitAndBeamline(processedDir)
            self.work_thread = XChemThread.choose_autoprocessing_outcome(os.path.join(self.database_directory,
                                                                                      self.data_source_file),
                                                                                      visit,
                                                                                      self.reference_file_list,
                                                                                      self.preferences,
                                                                                      self.initial_model_directory,
                                                                                      self.rescore,
                                                                                      self.xce_logfile)

            self.explorer_active = 1
            self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.connect(self.work_thread, QtCore.SIGNAL("populate_datasets_summary_table"),
                     self.populate_datasets_summary_table)
            self.work_thread.start()

    def update_gdaLog_parsing_instructions_and_score(self, gdaLogInstructions):
        self.gdaLogInstructions = gdaLogInstructions
        self.select_best_autoprocessing_result()

    def read_pinIDs_from_gda_logs(self):
        self.update_log.insert('reading pinIDs from gda logfiles...')
        visit, beamline = XChemMain.getVisitAndBeamline(self.beamline_directory)
        self.work_thread = XChemThread.read_pinIDs_from_gda_logs(beamline,
                                                                 visit,
                                                                 os.path.join(
                                                                     self.database_directory,
                                                                     self.data_source_file),
                                                                 self.gdaLogInstructions,
                                                                 self.xce_logfile)

        self.explorer_active = 1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.connect(self.work_thread, QtCore.SIGNAL("update_gdaLog_parsing_instructions_and_score"),
                     self.update_gdaLog_parsing_instructions_and_score)
        self.work_thread.start()

    def check_for_new_autoprocessing_results(self):
        self.update_log.insert('checking for new data collection')
        start_thread = True
        if start_thread:
            if self.target == '=== SELECT TARGET ===':
                start_thread = CheckAutoProcessing().query(self)

        if start_thread:
            processedDir=os.path.join(self.beamline_directory,'processed',self.target)
            self.work_thread = XChemThread.read_write_autoprocessing_results_from_to_disc(processedDir,
                                                                                          os.path.join(
                                                                                              self.database_directory,
                                                                                              self.data_source_file),
                                                                                          self.initial_model_directory,
                                                                                          self.xce_logfile)

            self.explorer_active = 1
            self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.connect(self.work_thread, QtCore.SIGNAL("read_pinIDs_from_gda_logs"),
                         self.read_pinIDs_from_gda_logs)
            self.work_thread.start()

    ####################################################################################################################
    #                                                                                                                  #
    #                                                 MAPS TAB                                                         #
    #                                                                                                                  #
    ####################################################################################################################
    def set_new_reference_if_applicable(self):
        print('hallo')
        reference_root = str(self.reference_file_selection_combobox.currentText())
        pg_ref = ''
        ucVol_ref = 0.0
        for reference in self.reference_file_list:
            print(reference[0], reference_root)
            if reference[0] == reference_root:
                pg_ref = reference[5]
                ucVol_ref = reference[4]
                break
        if ucVol_ref == 0.0:
            self.update_log.insert('cannot set reference file since unit cell volume of reference pdb is 0!')
            return

        for xtal in self.initial_model_dimple_dict:
            reference_file_selection_combobox = self.initial_model_dimple_dict[xtal][1]
            self.populate_reference_combobox(reference_file_selection_combobox)
            db_dict = self.xtal_db_dict[xtal]
            pg_xtal = db_dict['DataProcessingPointGroup']
            ucVol_xtal = db_dict['DataProcessingUnitCellVolume']

            try:
                difference = math.fabs(1 - (float(ucVol_xtal) / float(ucVol_ref))) * 100
            except ValueError:
                self.update_log.insert(xtal + ' -> cannot calculate unit cell volume difference')
                continue

            if pg_xtal == pg_ref and difference < self.allowed_unitcell_difference_percent:
                print(xtal, pg_xtal, ucVol_xtal)
                index = reference_file_selection_combobox.findText(reference_root, QtCore.Qt.MatchFixedString)
                reference_file_selection_combobox.setCurrentIndex(index)
                self.update_log.insert(xtal + ' -> setting ' + reference_root + ' as input PDB file for DIMPLE')

    def refresh_reference_file_list(self):
        self.reference_file_list = self.get_reference_file_list(' ')
        self.populate_reference_combobox(self.reference_file_selection_combobox)

    def on_context_menu_initial_model(self, point):
        # show context menu
        self.popMenu_for_maps_table.exec_(self.sender().mapToGlobal(point))

    ####################################################################################################################
    #                                                                                                                  #
    #                                                 PANDDA TAB                                                       #
    #                                                                                                                  #
    ####################################################################################################################
    def select_pandda_input_template(self):
        mtzin = ''
        filepath_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(self.window, 'Select Example PDB or MTZ File',
                                                                   self.initial_model_directory, '*.pdb;;*.mtz')
        filepath = str(tuple(filepath_temp)[0])
        pdbin = filepath.split('/')[-1]
        if filepath.endswith('.pdb'):
            pdbin = filepath.split('/')[-1]
            mtzin_temp = pdbin.replace('.pdb', '.mtz')
            if os.path.isfile(filepath.replace(pdbin, mtzin_temp)):
                mtzin = mtzin_temp
            else:
                mtzin = ''
        if filepath.endswith('.mtz'):
            mtzin = filepath.split('/')[-1]
            pdbin_temp = pdbin.replace('.mtz', '.pdb')
            if os.path.isfile(filepath.replace(mtzin, pdbin_temp)):
                pdbin = pdbin_temp
            else:
                pdbin = ''
        if len(filepath.split('/')) - len(self.initial_model_directory.split('/')) == 2:
            self.pandda_input_data_dir_entry.setText(os.path.join(self.initial_model_directory, '*'))
        elif len(filepath.split('/')) - len(self.initial_model_directory.split('/')) > 2:
            subdir = os.path.join(
                *filepath.split('/')[len(self.initial_model_directory.split('/')) + 1:len(filepath.split('/')) - 1])
            self.pandda_input_data_dir_entry.setText(os.path.join(self.initial_model_directory, '*', subdir))
        else:
            pass
        self.pandda_pdb_style_entry.setText(pdbin)
        self.pandda_mtz_style_entry.setText(mtzin)

    def change_pandda_spg_label(self):
        combo_text = str(self.pandda_reference_file_selection_combobox.currentText())
        for file in self.reference_file_list:
            if file[0] == combo_text:
                self.pandda_reference_file_spg_label.setText(file[1])
                break

    ####################################################################################################################
    #                                                                                                                  #
    #                                                 DEPO TAB                                                         #
    #                                                                                                                  #
    ####################################################################################################################
    def export_to_html(self):
        self.update_log.insert('exporting contents of SQLite database into ' + self.html_export_directory)
        os.system(
            'ccp4-python ' + os.getenv('XChemExplorer_DIR') + '/web/process_sqlite.py -t Summary -s ' + os.path.join(
                self.database_directory, self.data_source_file) + ' -d ' + self.html_export_directory)
        XChemWeb.create_ICM_input_file(self.html_export_directory,
                                       os.path.join(self.database_directory, self.data_source_file))
        self.update_log.insert('open ICMpro:')
        self.update_log.insert('/dls/science/groups/i04-1/software/icm-3.8-5/icm64 -g')
        self.update_log.insert('open file browser and navigate to ' + self.html_export_directory)
        self.update_log.insert('drag and drop dsEvent_sqlite.icm into the main window')
        self.update_log.insert('the script will appear in the Workspace Panel')
        self.update_log.insert('right click on the script and select RUN')
        self.update_log.insert('be patient, this may take a while, depending on the number of events')
        self.status_bar.showMessage('please check terminal window for further information')

    def open_icm(self):
        self.update_log.insert('starting ICM...')
        self.work_thread = XChemThread.start_ICM(self.html_export_directory)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def prepare_files_for_zenodo_upload(self):
        self.update_log.insert('preparing files for ZENODO upload...')
        os.system('ccp4-python ' + os.getenv(
            'XChemExplorer_DIR') + '/helpers/prepare_for_zenodo_upload.py ' + self.html_export_directory)

    def update_html_for_zenodo_upload(self):
        try:
            uploadID = int(self.zenodo_upload_id_entry.text())
            self.update_log.insert('updating html files for ZENODO upload,...')
            self.update_log.insert('ZENODO upload = ' + str(uploadID))
            os.system('ccp4-python ' + os.getenv(
                'XChemExplorer_DIR') + '/helpers/prepare_for_zenodo_upload.py {0!s} {1!s}'.format(
                self.html_export_directory, uploadID))
        except ValueError:
            self.update_log.insert('zenodo upload ID must be an integer!')

    ####################################################################################################################
    #                                                                                                                  #
    #                                                 SETTINGS TAB                                                     #
    #                                                                                                                  #
    ####################################################################################################################
    def settings_button_clicked(self):
        if self.sender().text() == 'Select Project Directory':
            self.initial_model_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.initial_model_directory_label.setText(self.initial_model_directory)
            self.pandda_input_data_dir_entry.setText(self.initial_model_directory)
            self.settings['initial_model_directory'] = self.initial_model_directory
        if self.sender().text() == 'Select Reference Structure Directory':
            reference_directory_temp = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            if reference_directory_temp != self.reference_directory:
                self.reference_directory = reference_directory_temp
                self.update_reference_files(' ')
            self.reference_directory_label.setText(self.reference_directory)
            self.settings['reference_directory'] = self.reference_directory
        if self.sender().text() == 'Select Data Source File':
            filepath_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(self.window, 'Select File',
                                                                       self.database_directory, '*.sqlite')
            filepath = str(tuple(filepath_temp)[0])
            self.data_source_file = filepath.split('/')[-1]
            self.database_directory = filepath[:filepath.rfind('/')]
            self.settings['database_directory'] = self.database_directory
            self.settings['data_source'] = os.path.join(self.database_directory, self.data_source_file)
            write_enabled = self.check_write_permissions_of_data_source()
            if not write_enabled:
                self.data_source_set = False
            else:
                self.data_source_set = True
                self.data_source_file_label.setText(os.path.join(self.database_directory, self.data_source_file))
                self.db = XChemDB.data_source(os.path.join(self.database_directory, self.data_source_file))
                self.db.create_missing_columns()
                self.datasource_menu_reload_samples()
        if self.sender().text() == 'Select Data Collection Directory':
            dir_name = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            if dir_name != self.beamline_directory:
                self.beamline_directory = dir_name
                self.target_list, self.visit_list = XChemMain.get_target_and_visit_list(self.beamline_directory)
                self.populate_target_selection_combobox(self.target_selection_combobox)
            self.beamline_directory_label.setText(self.beamline_directory)
            self.settings['beamline_directory'] = self.beamline_directory

        if self.sender().text() == 'Select Existing\nCollection Summary File':
            if self.datasets_summary_file != '':
                filepath_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(self.window, 'Select File',
                                                                           self.datasets_summary_file[
                                                                           :self.datasets_summary_file.rfind(
                                                                               '/')], '*.pkl')
            else:
                filepath_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(self.window, 'Select File', os.getcwd(),
                                                                           '*.pkl')
            filepath = str(tuple(filepath_temp)[0])
            self.datasets_summary_file = filepath
            self.datasets_summary_file_label.setText(self.datasets_summary_file)
            self.settings['datasets_summary'] = self.datasets_summary_file

        if self.sender().text() == 'Assign New\nCollection Summary File':
            if self.datasets_summary_file != '':
                file_name = str(QtGui.QFileDialog.getSaveFileName(self.window, 'New file',
                                                                  self.datasets_summary_file[
                                                                  :self.datasets_summary_file.rfind('/')]))
            else:
                file_name = str(QtGui.QFileDialog.getSaveFileName(self.window, 'New file', self.current_directory))
            # make sure that the file always has .pkl extension
            if str(file_name).rfind('.') != -1:
                file_name = file_name[:file_name.rfind('.')] + '.pkl'
            else:
                file_name = file_name + '.pkl'
            self.datasets_summary_file = file_name
            self.datasets_summary_file_label.setText(self.datasets_summary_file)
            self.settings['datasets_summary'] = self.datasets_summary_file

        if self.sender().text() == 'Select CCP4_SCR Directory':
            self.ccp4_scratch_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.ccp4_scratch_directory_label.setText(self.ccp4_scratch_directory)
            self.settings['ccp4_scratch'] = self.ccp4_scratch_directory
        if self.sender().text() == 'Select PANNDAs Directory':
            self.panddas_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.panddas_directory_label.setText(self.panddas_directory)
            self.pandda_output_data_dir_entry.setText(self.panddas_directory)
            print('PANDDA', self.panddas_directory)
            self.settings['panddas_directory'] = self.panddas_directory

            self.layout_funcs.pandda_html(self)

        if self.sender().text() == 'Select HTML Export Directory':
            self.html_export_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.html_export_directory_label.setText(self.html_export_directory)
            self.settings['html_export_directory'] = self.html_export_directory

        if self.sender().text() == 'Select Group deposition Directory':
            self.group_deposit_directory = str(QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory"))
            self.group_deposition_directory_label.setText(self.group_deposit_directory)
            self.settings['group_deposit_directory'] = self.group_deposit_directory

        #self.datasource_menu_reload_samples()



    ######################################### sort stuff below here ####################################################
    def select_sample_for_dimple(self):
        indexes = self.maps_table.selectionModel().selectedRows()
        for index in sorted(indexes):
            xtal = str(self.maps_table.item(index.row(), 0).text())
            self.update_log.insert('{0!s} is marked for DIMPLE'.format(index.row()))
            self.initial_model_dimple_dict[xtal][0].setChecked(True)

    def update_summary_plot(self):
        if self.data_source_set:
            XChemPlots.summary_plot(os.path.join(self.database_directory, self.data_source_file),
                                    self.overview_axes).update_overview()
            self.overview_canvas.draw()

    def show_preferences(self):
        preferences = QtGui.QMessageBox()
        preferencesLayout = preferences.layout()

        vbox = QtGui.QVBoxLayout()
        settings_hbox_filename_root = QtGui.QHBoxLayout()
        filename_root_label = QtGui.QLabel('filename root:')
        settings_hbox_filename_root.addWidget(filename_root_label)
        filename_root_input = QtGui.QLineEdit()
        filename_root_input.setFixedWidth(400)
        filename_root_input.setText(str(self.filename_root))
        filename_root_input.textChanged[str].connect(self.change_filename_root)
        settings_hbox_filename_root.addWidget(filename_root_input)
        vbox.addLayout(settings_hbox_filename_root)

        settings_hbox_adjust_allowed_unit_cell_difference = QtGui.QHBoxLayout()
        adjust_allowed_unit_cell_difference_label = QtGui.QLabel(
            'Max. Allowed Unit Cell Difference between Reference and Target (%):')
        settings_hbox_adjust_allowed_unit_cell_difference.addWidget(adjust_allowed_unit_cell_difference_label)
        adjust_allowed_unit_cell_difference = QtGui.QLineEdit()
        adjust_allowed_unit_cell_difference.setFixedWidth(200)
        adjust_allowed_unit_cell_difference.setText(str(self.allowed_unitcell_difference_percent))
        adjust_allowed_unit_cell_difference.textChanged[str].connect(self.change_allowed_unitcell_difference_percent)
        settings_hbox_adjust_allowed_unit_cell_difference.addWidget(adjust_allowed_unit_cell_difference)
        vbox.addLayout(settings_hbox_adjust_allowed_unit_cell_difference)

        settings_hbox_acceptable_low_resolution_limit = QtGui.QHBoxLayout()
        adjust_acceptable_low_resolution_limit_label = QtGui.QLabel(
            'Acceptable low resolution limit for datasets (in Angstrom):')
        settings_hbox_acceptable_low_resolution_limit.addWidget(adjust_acceptable_low_resolution_limit_label)
        adjust_acceptable_low_resolution_limit = QtGui.QLineEdit()
        adjust_acceptable_low_resolution_limit.setFixedWidth(200)
        adjust_acceptable_low_resolution_limit.setText(str(self.acceptable_low_resolution_limit_for_data))
        adjust_acceptable_low_resolution_limit.textChanged[str].connect(self.change_acceptable_low_resolution_limit)
        settings_hbox_acceptable_low_resolution_limit.addWidget(adjust_acceptable_low_resolution_limit)
        vbox.addLayout(settings_hbox_acceptable_low_resolution_limit)

        vbox_data = QtGui.QVBoxLayout()
        vbox_data.addWidget(
            QtGui.QLabel('Select amount of processed data you wish to copy to initial_model directory:'))
        self.preferences_data_to_copy_combobox = QtGui.QComboBox()
        for item in self.preferences_data_to_copy:
            self.preferences_data_to_copy_combobox.addItem(item[0])
        self.preferences_data_to_copy_combobox.currentIndexChanged.connect(
            self.preferences_data_to_copy_combobox_changed)
        vbox_data.addWidget(self.preferences_data_to_copy_combobox)
        vbox.addLayout(vbox_data)

        vbox_select = QtGui.QVBoxLayout()
        vbox_select.addWidget(QtGui.QLabel('Dataset Selection Mechanism:'))
        self.preferences_selection_mechanism_combobox = QtGui.QComboBox()
        for item in self.preferences_selection_mechanism:
            self.preferences_selection_mechanism_combobox.addItem(item)
        self.preferences_selection_mechanism_combobox.currentIndexChanged.connect(
            self.preferences_selection_mechanism_combobox_changed)
        vbox_select.addWidget(self.preferences_selection_mechanism_combobox)
        vbox.addLayout(vbox_select)

        vbox_restraints = QtGui.QVBoxLayout()
        vbox_restraints.addWidget(QtGui.QLabel('Restraints generation program:'))
        self.preferences_restraints_generation_combobox = QtGui.QComboBox()
        program_list = []

        if self.external_software['acedrg']:
            program_list.append('acedrg')
            self.restraints_program = 'acedrg'
        if self.external_software['phenix.elbow']: program_list.append('phenix.elbow')
        if self.external_software['grade']:        program_list.append('grade')
        for item in program_list:
            self.preferences_restraints_generation_combobox.addItem(item)
        self.preferences_restraints_generation_combobox.currentIndexChanged.connect(
            self.preferences_restraints_generation_combobox_changed)
        index = self.preferences_restraints_generation_combobox.findText(self.restraints_program,
                                                                         QtCore.Qt.MatchFixedString)
        self.preferences_restraints_generation_combobox.setCurrentIndex(index)
        vbox_restraints.addWidget(self.preferences_restraints_generation_combobox)
        vbox.addLayout(vbox_restraints)

        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(QtGui.QLabel('XCE logfile:'))
        self.xce_logfile_label = QtGui.QLabel(self.xce_logfile)
        hbox.addWidget(self.xce_logfile_label)
        button = QtGui.QPushButton("Change")
        button.clicked.connect(self.set_xce_logfile)
        hbox.addWidget(button)
        vbox.addLayout(hbox)

        settings_hbox_max_queue_jobs = QtGui.QHBoxLayout()
        adjust_max_queue_jobs_label = QtGui.QLabel('Max. number of jobs running at once on DLS cluster:')
        settings_hbox_max_queue_jobs.addWidget(adjust_max_queue_jobs_label)
        adjust_max_queue_jobs = QtGui.QLineEdit()
        adjust_max_queue_jobs.setFixedWidth(200)
        adjust_max_queue_jobs.setText(str(self.max_queue_jobs))
        adjust_max_queue_jobs.textChanged[str].connect(self.change_max_queue_jobs)
        settings_hbox_max_queue_jobs.addWidget(adjust_max_queue_jobs)
        vbox.addLayout(settings_hbox_max_queue_jobs)

        settings_hbox_remote_qsub = QtGui.QHBoxLayout()
        remote_qsub_label = QtGui.QLabel('remote qsub:')
        settings_hbox_remote_qsub.addWidget(remote_qsub_label)
        self.remote_qsub_checkbox = QtGui.QCheckBox('use')
        #self.remote_qsub_checkbox.toggled.connect(self.run_qsub_remotely)

        if self.using_remote_qsub_submission:
            self.remote_qsub_checkbox.setChecked(True)
        settings_hbox_remote_qsub.addWidget(self.remote_qsub_checkbox)
        self.remote_qsub_command = QtGui.QLineEdit()
        self.remote_qsub_command.setFixedWidth(550)
        self.remote_qsub_command.setText(self.remote_qsub_submission)
        settings_hbox_remote_qsub.addWidget(self.remote_qsub_command)
        vbox.addLayout(settings_hbox_remote_qsub)

        apply_button = QtGui.QPushButton('Apply')
        apply_button.clicked.connect(self.run_qsub_remotely)
        settings_hbox_remote_qsub.addWidget(apply_button)


        preferencesLayout.addLayout(vbox, 0, 0)

        preferences.exec_();

    def run_qsub_remotely(self):
        self.remote_qsub_submission = str(self.remote_qsub_command.text())
        print(str(self.remote_qsub_submission))
        if self.remote_qsub_checkbox.isChecked():
            self.update_log.insert('submitting jobs to remote machine with: %s' % self.remote_qsub_submission)
            self.external_software['qsub_remote'] = self.remote_qsub_submission
            self.using_remote_qsub_submission = True
            self.settings['remote_qsub'] = self.remote_qsub_submission
        else:
            self.update_log.insert('switching off remote job submission')
            self.external_software['qsub_remote'] = ''
            self.settings['remote_qsub'] = ''
            self.using_remote_qsub_submission = False

    def enter_pdb_codes(self):
        pdbID_entry = QtGui.QMessageBox()
        pdbID_entryLayout = pdbID_entry.layout()

        vbox = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel('Text from PDB email'), 0, 0)
        self.pdb_code_entry = QtGui.QTextEdit()
        self.pdb_code_entry.setText('')
        self.pdb_code_entry.setFixedWidth(500)
        grid.addWidget(self.pdb_code_entry, 1, 0, 20, 1)

        frame.setLayout(grid)
        vbox.addWidget(frame)

        hbox = QtGui.QHBoxLayout()
        button = QtGui.QPushButton('Update Database')
        button.clicked.connect(self.update_database_with_pdb_codes)
        hbox.addWidget(button)

        vbox.addLayout(hbox)
        pdbID_entryLayout.addLayout(vbox, 0, 0)
        pdbID_entry.exec_();

    def create_missing_apo_records_in_depositTable(self):
        self.db.create_missing_apo_records_for_all_structures_in_depositTable(self.initial_model_directory,
                                                                              self.xce_logfile)

    def update_file_information_of_apo_records(self):
        XChemDeposit.update_file_locations_of_apo_structuresin_DB(
            os.path.join(self.database_directory, self.data_source_file), self.initial_model_directory,
            self.xce_logfile)

    def prepare_models_for_deposition_ligand_bound(self):

        structureType = "ligand_bound"

        overwrite_existing_mmcif = True
        self.work_thread = XChemDeposit.prepare_mmcif_files_for_deposition(
            os.path.join(self.database_directory, self.data_source_file),
            self.xce_logfile,
            overwrite_existing_mmcif,
            self.initial_model_directory,
            structureType)
        self.explorer_active = 1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def prepare_models_for_deposition_apo(self):

        structureType = "apo"

        overwrite_existing_mmcif = True
        self.work_thread = XChemDeposit.prepare_mmcif_files_for_deposition(
            os.path.join(self.database_directory, self.data_source_file),
            self.xce_logfile,
            overwrite_existing_mmcif,
            self.initial_model_directory,
            structureType)
        self.explorer_active = 1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()


    def prepare_for_group_deposition_upload(self):

        self.work_thread = XChemDeposit.prepare_for_group_deposition_upload(
            os.path.join(self.database_directory, self.data_source_file),
            self.xce_logfile,
            self.group_deposit_directory)
        self.explorer_active = 1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def check_smiles_in_db_and_pdb(self):

        self.work_thread = XChemDeposit.compare_smiles_in_db_with_ligand_in_pdb(self.initial_model_directory,
                                                                                os.path.join(self.database_directory,
                                                                                             self.data_source_file),
                                                                                self.xce_logfile)
        self.explorer_active = 1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.connect(self.work_thread, QtCore.SIGNAL("show_error_dict"), self.show_error_dict)
        self.work_thread.start()

    def edit_deposition_data(self):
        DepositionDetails().setup(self)
        self.messageBox.exec_()

    def save_deposit_config_file(self):
        self.update_deposit_dict()
        file_name = str(QtGui.QFileDialog.getSaveFileName(self.window, 'Save file', self.current_directory))
        # make sure that the file always has .deposit extension
        if str(file_name).rfind('.') != -1:
            file_name = file_name[:file_name.rfind('.')] + '.deposit'
        else:
            file_name = file_name + '.deposit'
        pickle.dump(self.deposit_dict, open(file_name, 'wb'))

    def update_database_with_pdb_codes(self):
        self.work_thread = XChemDeposit.import_PDB_IDs(str(self.pdb_code_entry.toPlainText()),
                                                       os.path.join(self.database_directory, self.data_source_file),
                                                       self.xce_logfile)
        self.explorer_active = 1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def load_deposit_config_file(self):
        file_name_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(self.window, 'Open file', self.current_directory,
                                                                    '*.deposit')
        file_name = tuple(file_name_temp)[0]
        self.deposit_dict = pickle.load(open(file_name, "rb"))
        self.update_deposit_input()

    def load_deposit_from_database(self):
        print('hallo')

    def save_deposit_to_database(self):
        self.update_deposit_dict()
        msgBox = QtGui.QMessageBox()
        msgBox.setText(
            "*** WARNING ***\nAre you sure you want to update the database?\nThis will overwrite previous entries!")
        msgBox.addButton(QtGui.QPushButton('Yes'), QtGui.QMessageBox.YesRole)
        msgBox.addButton(QtGui.QPushButton('No'), QtGui.QMessageBox.RejectRole)
        reply = msgBox.exec_();
        if reply == 0:
            self.work_thread = XChemDeposit.update_depositTable(self.deposit_dict,
                                                                os.path.join(self.database_directory,
                                                                             self.data_source_file),
                                                                self.xce_logfile)
            self.explorer_active = 1
            self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.work_thread.start()

    def update_deposit_input(self):
        try:
            self.contact_author_PI_salutation.setText(self.deposit_dict['contact_author_PI_salutation'])
            self.contact_author_PI_first_name.setText(self.deposit_dict['contact_author_PI_first_name'])
            self.contact_author_PI_last_name.setText(self.deposit_dict['contact_author_PI_last_name'])
            self.contact_author_PI_middle_name.setText(self.deposit_dict['contact_author_PI_middle_name'])
            index = self.contact_author_PI_role.findText(self.deposit_dict['contact_author_PI_role'],
                                                         QtCore.Qt.MatchFixedString)
            self.contact_author_PI_role.setCurrentIndex(index)
            index = self.contact_author_PI_organization_type.findText(
                self.deposit_dict['contact_author_PI_organization_type'], QtCore.Qt.MatchFixedString)
            self.contact_author_PI_organization_type.setCurrentIndex(index)
            self.contact_author_PI_organization_name.setText(self.deposit_dict['contact_author_PI_organization_name'])
            self.contact_author_PI_email.setText(self.deposit_dict['contact_author_PI_email'])
            self.contact_author_PI_address.setText(self.deposit_dict['contact_author_PI_address'])
            self.contact_author_PI_city.setText(self.deposit_dict['contact_author_PI_city'])
            self.contact_author_PI_State_or_Province.setText(self.deposit_dict['contact_author_PI_State_or_Province'])
            self.contact_author_PI_Zip_Code.setText(self.deposit_dict['contact_author_PI_Zip_Code'])
            self.contact_author_PI_Country.setText(self.deposit_dict['contact_author_PI_Country'])
            self.contact_author_PI_phone_number.setText(self.deposit_dict['contact_author_PI_phone_number'])

            self.contact_author_salutation.setText(self.deposit_dict['contact_author_salutation'])
            self.contact_author_first_name.setText(self.deposit_dict['contact_author_first_name'])
            self.contact_author_last_name.setText(self.deposit_dict['contact_author_last_name'])
            self.contact_author_middle_name.setText(self.deposit_dict['contact_author_middle_name'])
            index = self.contact_author_role.findText(self.deposit_dict['contact_author_role'],
                                                      QtCore.Qt.MatchFixedString)
            self.contact_author_role.setCurrentIndex(index)
            index = self.contact_author_organization_type.findText(
                self.deposit_dict['contact_author_organization_type'], QtCore.Qt.MatchFixedString)
            self.contact_author_organization_type.setCurrentIndex(index)
            self.contact_author_organization_name.setText(self.deposit_dict['contact_author_organization_name'])
            self.contact_author_email.setText(self.deposit_dict['contact_author_email'])
            self.contact_author_address.setText(self.deposit_dict['contact_author_address'])
            self.contact_author_city.setText(self.deposit_dict['contact_author_city'])
            self.contact_author_State_or_Province.setText(self.deposit_dict['contact_author_State_or_Province'])
            self.contact_author_Zip_Code.setText(self.deposit_dict['contact_author_Zip_Code'])
            self.contact_author_Country.setText(self.deposit_dict['contact_author_Country'])
            self.contact_author_phone_number.setText(self.deposit_dict['contact_author_phone_number'])
            index = self.Release_status_for_coordinates.findText(self.deposit_dict['Release_status_for_coordinates'],
                                                                 QtCore.Qt.MatchFixedString)
            self.Release_status_for_coordinates.setCurrentIndex(index)
            index = self.Release_status_for_sequence.findText(self.deposit_dict['Release_status_for_sequence'],
                                                              QtCore.Qt.MatchFixedString)
            self.Release_status_for_sequence.setCurrentIndex(index)

            self.group_deposition_title.setText(self.deposit_dict['group_deposition_title'])
            self.group_description.setText(self.deposit_dict['group_description'])

            self.structure_title.setText(self.deposit_dict['structure_title'])
            self.structure_title_apo.setText(self.deposit_dict['structure_title_apo'])

            for n, name in enumerate(self.deposit_dict['structure_author_name'].split(';')):
                self.structure_author_name_List[n].setText(name)

            self.primary_citation_id.setText(self.deposit_dict['primary_citation_id'])
            self.primary_citation_journal_abbrev.setText(self.deposit_dict['primary_citation_journal_abbrev'])
            self.primary_citation_title.setText(self.deposit_dict['primary_citation_title'])
            self.primary_citation_year.setText(self.deposit_dict['primary_citation_year'])
            self.primary_citation_journal_volume.setText(self.deposit_dict['primary_citation_journal_volume'])
            self.primary_citation_page_first.setText(self.deposit_dict['primary_citation_page_first'])
            self.primary_citation_page_last.setText(self.deposit_dict['primary_citation_page_last'])

            for n, name in enumerate(self.deposit_dict['primary_citation_author_name'].split(';')):
                self.primary_citation_author_name_List[n].setText(name)

            self.molecule_name.setText(self.deposit_dict['molecule_name'])
            self.fragment_name_one_specific_mutation.setText(self.deposit_dict['fragment_name_one_specific_mutation'])
            index = self.Source_organism_scientific_name.findText(self.deposit_dict['Source_organism_scientific_name'],
                                                                  QtCore.Qt.MatchFixedString)
            self.Source_organism_scientific_name.setCurrentIndex(index)

            self.Source_organism_gene.setText(self.deposit_dict['Source_organism_gene'])
            self.Source_organism_strain.setText(self.deposit_dict['Source_organism_strain'])
            index = self.Expression_system_scientific_name.findText(
                self.deposit_dict['Expression_system_scientific_name'], QtCore.Qt.MatchFixedString)
            self.Expression_system_scientific_name.setCurrentIndex(index)

            self.Expression_system_strain.setText(self.deposit_dict['Expression_system_strain'])
            self.Expression_system_vector_type.setText(self.deposit_dict['Expression_system_vector_type'])
            self.Expression_system_plasmid_name.setText(self.deposit_dict['Expression_system_plasmid_name'])
            self.Manipulated_source_details.setText(self.deposit_dict['Manipulated_source_details'])

            self.structure_keywords.setText(self.deposit_dict['structure_keywords'])
            self.biological_assembly_chain_number.setText(self.deposit_dict['biological_assembly_chain_number'])
            self.molecule_one_letter_sequence_uniprot_id.setText(
                self.deposit_dict['molecule_one_letter_sequence_uniprot_id'])
            self.molecule_one_letter_sequence.setText(self.deposit_dict['molecule_one_letter_sequence'])
            self.SG_project_name.setText(self.deposit_dict['SG_project_name'])
            self.full_name_of_SG_center.setText(self.deposit_dict['full_name_of_SG_center'])

            index = self.crystallization_method.findText(self.deposit_dict['crystallization_method'],
                                                         QtCore.Qt.MatchFixedString)
            self.crystallization_method.setCurrentIndex(index)

            self.crystallization_pH.setText(self.deposit_dict['crystallization_pH'])
            self.crystallization_temperature.setText(self.deposit_dict['crystallization_temperature'])
            self.crystallization_details.setText(self.deposit_dict['crystallization_details'])
            index = self.radiation_source.findText(self.deposit_dict['radiation_source'], QtCore.Qt.MatchFixedString)
            self.radiation_source.setCurrentIndex(index)

            index = self.radiation_source_type.findText(self.deposit_dict['radiation_source_type'],
                                                        QtCore.Qt.MatchFixedString)
            self.radiation_source_type.setCurrentIndex(index)

            self.radiation_wavelengths.setText(self.deposit_dict['radiation_wavelengths'])
            index = self.radiation_detector.findText(self.deposit_dict['radiation_detector'],
                                                     QtCore.Qt.MatchFixedString)
            self.radiation_detector.setCurrentIndex(index)

            index = self.radiation_detector_type.findText(self.deposit_dict['radiation_detector_type'],
                                                          QtCore.Qt.MatchFixedString)
            self.radiation_detector_type.setCurrentIndex(index)

            self.data_collection_date.setText(self.deposit_dict['data_collection_date'])
            self.data_collection_temperature.setText(self.deposit_dict['data_collection_temperature'])
            self.data_collection_protocol.setText(self.deposit_dict['data_collection_protocol'])

            self.pdbx_starting_model.setText(self.deposit_dict['pdbx_starting_model'])
            index = self.data_integration_software.findText(self.deposit_dict['data_integration_software'],
                                                            QtCore.Qt.MatchFixedString)
            self.data_integration_software.setCurrentIndex(index)
            index = self.phasing_software.findText(self.deposit_dict['phasing_software'], QtCore.Qt.MatchFixedString)
            self.phasing_software.setCurrentIndex(index)

        except ValueError:
            self.update_status_bar('Sorry, this is not a XChemExplorer deposit file!')
            self.update_log.insert('Sorry, this is not a XChemExplorer deposit file!')

    def update_deposit_dict(self):
        self.deposit_dict = {
            'contact_author_PI_salutation': str(self.contact_author_PI_salutation.text()),
            'contact_author_PI_first_name': str(self.contact_author_PI_first_name.text()),
            'contact_author_PI_last_name': str(self.contact_author_PI_last_name.text()),
            'contact_author_PI_middle_name': str(self.contact_author_PI_middle_name.text()),
            'contact_author_PI_role': str(self.contact_author_PI_role.currentText()),
            'contact_author_PI_organization_type': str(self.contact_author_PI_organization_type.currentText()),
            'contact_author_PI_organization_name': str(self.contact_author_PI_organization_name.text()),
            'contact_author_PI_email': str(self.contact_author_PI_email.text()),
            'contact_author_PI_address': str(self.contact_author_PI_address.text()),
            'contact_author_PI_city': str(self.contact_author_PI_city.text()),
            'contact_author_PI_State_or_Province': str(self.contact_author_PI_State_or_Province.text()),
            'contact_author_PI_Zip_Code': str(self.contact_author_PI_Zip_Code.text()),
            'contact_author_PI_Country': str(self.contact_author_PI_Country.text()),
            'contact_author_PI_phone_number': str(self.contact_author_PI_phone_number.text()),

            'contact_author_salutation': str(self.contact_author_salutation.text()),
            'contact_author_first_name': str(self.contact_author_first_name.text()),
            'contact_author_last_name': str(self.contact_author_last_name.text()),
            'contact_author_middle_name': str(self.contact_author_middle_name.text()),
            'contact_author_role': str(self.contact_author_role.currentText()),
            'contact_author_organization_type': str(self.contact_author_organization_type.currentText()),
            'contact_author_organization_name': str(self.contact_author_organization_name.text()),
            'contact_author_email': str(self.contact_author_email.text()),
            'contact_author_address': str(self.contact_author_address.text()),
            'contact_author_city': str(self.contact_author_city.text()),
            'contact_author_State_or_Province': str(self.contact_author_State_or_Province.text()),
            'contact_author_Zip_Code': str(self.contact_author_Zip_Code.text()),
            'contact_author_Country': str(self.contact_author_Country.text()),
            'contact_author_phone_number': str(self.contact_author_phone_number.text()),

            'Release_status_for_coordinates': str(self.Release_status_for_coordinates.currentText()),
            'Release_status_for_sequence': str(self.Release_status_for_sequence.currentText()),

            'group_deposition_title': str(self.group_deposition_title.text()),
            'group_description': str(self.group_description.text()),

            'structure_title': str(self.structure_title.text()),
            'structure_title_apo': str(self.structure_title_apo.text()),

            'primary_citation_id': str(self.primary_citation_id.text()),
            'primary_citation_journal_abbrev': str(self.primary_citation_journal_abbrev.text()),
            'primary_citation_title': str(self.primary_citation_title.text()),
            'primary_citation_year': str(self.primary_citation_year.text()),
            'primary_citation_journal_volume': str(self.primary_citation_journal_volume.text()),
            'primary_citation_page_first': str(self.primary_citation_page_first.text()),
            'primary_citation_page_last': str(self.primary_citation_page_last.text()),

            'molecule_name': str(self.molecule_name.text()),
            'Source_organism_scientific_name': str(self.Source_organism_scientific_name.currentText()),
            'Source_organism_gene': str(self.Source_organism_gene.text()),
            'Source_organism_strain': str(self.Source_organism_strain.text()),
            'Expression_system_scientific_name': str(self.Expression_system_scientific_name.currentText()),
            'Expression_system_strain': str(self.Expression_system_strain.text()),
            'Expression_system_plasmid_name': str(self.Expression_system_plasmid_name.text()),
            'Expression_system_vector_type': str(self.Expression_system_vector_type.text()),
            'Manipulated_source_details': str(self.Manipulated_source_details.text()),
            'fragment_name_one_specific_mutation': str(self.fragment_name_one_specific_mutation.text()),

            'structure_keywords': str(self.structure_keywords.text()),
            'biological_assembly_chain_number': str(self.biological_assembly_chain_number.text()),
            'molecule_one_letter_sequence_uniprot_id': str(self.molecule_one_letter_sequence_uniprot_id.text()),
            'SG_project_name': str(self.SG_project_name.text()),
            'full_name_of_SG_center': str(self.full_name_of_SG_center.text()),
            'molecule_one_letter_sequence': str(self.molecule_one_letter_sequence.toPlainText()).replace(' ',
                                                                                                         '').replace(
                '\n', '').replace('\r', ''),

            'crystallization_method': str(self.crystallization_method.currentText()),
            'crystallization_pH': str(self.crystallization_pH.text()),
            'crystallization_temperature': str(self.crystallization_temperature.text()),
            'crystallization_details': str(self.crystallization_details.text()),

            'radiation_source': str(self.radiation_source.currentText()),
            'radiation_source_type': str(self.radiation_source_type.currentText()),
            'radiation_wavelengths': str(self.radiation_wavelengths.text()),
            'radiation_detector': str(self.radiation_detector.currentText()),
            'radiation_detector_type': str(self.radiation_detector_type.currentText()),
            'data_collection_date': str(self.data_collection_date.text()),
            'data_collection_temperature': str(self.data_collection_temperature.text()),
            'data_collection_protocol': str(self.data_collection_protocol.text()),
            'pdbx_starting_model': str(self.pdbx_starting_model.text()),
            'data_integration_software': str(self.data_integration_software.currentText()),
            'phasing_software': str(self.phasing_software.currentText())
        }

        structure_author_name = ''
        for widget in self.structure_author_name_List:
            structure_author_name += str(widget.text()) + ';'
        self.deposit_dict['structure_author_name'] = structure_author_name[:-1]

        primary_citation_author_name = ''
        for widget in self.primary_citation_author_name_List:
            primary_citation_author_name += str(widget.text()) + ';'
        self.deposit_dict['primary_citation_author_name'] = primary_citation_author_name[:-1]

    def set_primary_citation_as_structure_authors(self, state):
        if state == QtCore.Qt.Checked:
            for n, entry in enumerate(self.structure_author_name_List):
                self.primary_citation_author_name_List[n].setText(str(entry.text()))
        else:
            for n, entry in enumerate(self.primary_citation_author_name_List):
                entry.setText('')

    def set_xce_logfile(self):
        file_name = str(QtGui.QFileDialog.getSaveFileName(self.window, 'Save file', self.current_directory))
        self.xce_logfile = str(file_name)
        self.xce_logfile_label.setText(str(self.xce_logfile))
        if self.xce_logfile == '' or self.xce_logfile[self.xce_logfile.rfind('/') + 1:] == '':
            print('==> XCE: invalid file format')
        else:
            XChemLog.startLog(self.xce_logfile).create_logfile(self.xce_version)
            self.update_log = XChemLog.updateLog(self.xce_logfile)

    def select_datasource_columns_to_display(self):
        columns_to_show = QtGui.QMessageBox()
        columns_to_showLayout = columns_to_show.layout()
        columns_in_data_source = self.db.return_column_list()
        try:
            columns_in_data_source = self.db.return_column_list()
        except AttributeError:
            print('==> XCE: please select a datasource file')
            self.status_bar.showMessage('please select a datasource file')
            return

        column_dict = {}
        vbox = QtGui.QVBoxLayout()
        number_of_entries = len(columns_in_data_source)
        columns_shown_in_dialog_column = 15
        grid = QtGui.QGridLayout()
        x = 0
        y = 0
        columns_to_ignore = self.db.columns_not_to_display()
        for entries_added in range(number_of_entries):
            if not columns_in_data_source[entries_added][1] in columns_to_ignore:
                data_source_column = QtGui.QCheckBox(columns_in_data_source[entries_added][1])
                column_dict[entries_added] = data_source_column
                if columns_in_data_source[entries_added][1] in self.overview_datasource_table_columns:
                    data_source_column.setChecked(True)
                grid.addWidget(data_source_column, y, x)
                y += 1
            if y == columns_shown_in_dialog_column:
                y = 0
                x += 1
        vbox.addLayout(grid)
        columns_to_showLayout.addLayout(vbox, 0, 0)

        columns_to_show.addButton(QtGui.QPushButton('OK'), QtGui.QMessageBox.YesRole)
        columns_to_show.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)
        reply = columns_to_show.exec_();
        if reply == 0:
            columns_to_show_list = ['Sample ID']
            for key in column_dict:
                if column_dict[key].isChecked():
                    columns_to_show_list.append(columns_in_data_source[key][1])
            self.overview_datasource_table_columns = columns_to_show_list
            self.populate_and_update_datasource_table()

    def update_header_and_data_from_datasource(self):
        self.update_log.insert('getting information for all samples from data source...')
        self.db = XChemDB.data_source(os.path.join(self.database_directory, self.data_source_file))
        self.update_log.insert('creating missing columns in data source')
        self.db.create_missing_columns()
        self.update_log.insert('load header and data from data source')
        self.header, self.data = self.db.load_samples_from_data_source()
        self.update_log.insert('get all samples in data source')
        all_samples_in_db = self.db.execute_statement("select CrystalName from mainTable where CrystalName is not '';")

        self.xtal_db_dict = {}
        sampleID_column = 0
        for n, entry in enumerate(self.header):
            if entry == 'CrystalName':
                sampleID_column = n
                break
        for line in self.data:
            if str(line[sampleID_column]) != '':
                db_dict = {}
                for n, entry in enumerate(line):
                    if n != sampleID_column:
                        db_dict[str(self.header[n])] = str(entry)
                self.xtal_db_dict[str(line[sampleID_column])] = db_dict

        print('==> XCE: found ' + str(len(self.xtal_db_dict)) + ' samples')

    def datasource_menu_save_samples(self):
        print('hallo')

    def datasource_menu_export_csv_file(self):
        file_name = str(QtGui.QFileDialog.getSaveFileName(self.window, 'Save file', self.database_directory))
        if file_name.rfind('.') != -1:
            file_name = file_name[:file_name.rfind('.')] + '.csv'
        else:
            file_name = file_name + '.csv'
        self.db.export_to_csv_file(file_name)

    def datasource_menu_import_csv_file(self):
        if self.data_source_set:
            file_name = QtGui.QFileDialog.getOpenFileName(self.window, 'Open file', self.database_directory)
            self.db.import_csv_file(file_name)
        else:
            self.update_status_bar('Please load a data source file first')

    def datasource_menu_update_datasource(self):
        self.work_thread = XChemThread.synchronise_db_and_filesystem(self.initial_model_directory,
                                                                     os.path.join(self.database_directory,
                                                                                  self.data_source_file),
                                                                     self.panddas_directory, self.xce_logfile,
                                                                     'project_directory')
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.connect(self.work_thread, QtCore.SIGNAL("datasource_menu_reload_samples"),
                     self.datasource_menu_reload_samples)
        self.work_thread.start()

    def export_data_for_WONKA(self):
        self.update_log.insert('exporting CSV file for input into WONKA')
        self.db.export_csv_for_WONKA()

    def on_context_menu(self, point):
        # show context menu
        for key in self.dewar_configuration_dict:
            if self.dewar_configuration_dict[key] == self.sender():
                self.dewar_label_active = key
        self.popMenu.exec_(self.sender().mapToGlobal(point))



    def on_context_menu_reprocess_data(self, point):
        # show context menu
        self.popMenu_for_datasets_reprocess_table.exec_(self.sender().mapToGlobal(point))

    def flag_sample_for_recollection(self):
        self.dewar_configuration_dict[self.dewar_label_active].setStyleSheet("background-color: yellow")

    def undo_flag_sample_for_recollection(self):
        self.dewar_configuration_dict[self.dewar_label_active].setStyleSheet("background-color: gray")

    def show_html_summary_in_firefox(self, xtal):
        html_summary = self.albula_button_dict[xtal][2]
        print('html_summary', html_summary)
        new = 2
        webbrowser.open(html_summary, new=new)

    def update_pandda_crystal_from_combobox(self):
        self.pandda_analyse_crystal_from_selection_combobox.clear()
        self.pandda_analyse_crystal_from_selection_combobox.addItem('use all datasets')
        if os.path.isfile(os.path.join(self.database_directory, self.data_source_file)):
            self.load_crystal_form_from_datasource()
            if self.xtalform_dict != {}:
                print(self.xtalform_dict)
                for key in self.xtalform_dict:
                    self.pandda_analyse_crystal_from_selection_combobox.addItem(key)

    def populate_reference_combobox(self, combobox):
        combobox.clear()
        for reference_file in self.reference_file_list:
            combobox.addItem(reference_file[0])



    def populate_refinement_outcome_combobox(self, combobox):
        combobox.clear()
        for stage in self.refinement_stage:
            combobox.addItem(stage)



    def populate_target_selection_combobox(self, combobox):
        combobox.clear()
        for target in self.target_list:
            combobox.addItem(target)

    def combo_selected(self, text):
        self.map_url = str(self.panddas_directory + '/analyses/html_summaries/pandda_map_' + text + '.html')
        self.pandda_maps_html.load(QtCore.QUrl(self.map_url))
        self.pandda_maps_html.show()

    def add_map_html(self):
        self.map_list = glob.glob(str(self.panddas_directory + '/analyses/html_summaries/pandda_map_*.html'))
        self.list_options = []
        for i in range(0, len(self.map_list)):
            string = self.map_list[i]
            string = string.replace('/analyses/html_summaries/pandda_map_', '')
            string = string.replace('.html', '')
            string = string.replace(self.panddas_directory, '')
            self.list_options.append(string)
        self.pandda_map_list.clear()
        for i in range(0, len(self.list_options)):
            self.pandda_map_list.addItem(self.list_options[i])
        self.connect(self.pandda_map_list, QtCore.SIGNAL('activated(QString)'), self.combo_selected)

    def open_config_file(self):
        file_name_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(self.window, 'Open file', self.current_directory,
                                                                    '*.conf')
        file_name = tuple(file_name_temp)[0]

        try:
            pickled_settings = pickle.load(open(file_name, 'rb'))

        except:
             print('==> XCE: failed to open config file...')

        key_list = {#'beamline_directory': 'beamline_directory',
                    'initial_model_directory': 'initial_model_directory',
                   'panddas_directory': 'panddas_directory',
                   'html_export_directory': 'html_export_directory',
                    'group_deposit_directory': 'group_deposit_directory',
                   'database_directory': 'database_directory',
                   'datasets_summary_file': 'datasets_summary',
                   #"'data_source_file': 'data_source',
                   'ccp4_scratch_directory': 'ccp4_scratch',
                    'allowed_unitcell_difference_percent': 'unitcell_difference',
                   'acceptable_low_resolution_limit_for_data': 'too_low_resolution_data',
                   #'reference_directory_temp': 'reference_directory'
                     }
        self.pandda_input_data_dir_entry.setText(os.path.join(self.initial_model_directory, '*'))

        for current_key in key_list:
            try:
                command = str('self.' + current_key + " = pickled_settings['" + key_list[current_key] +"']")
                exec(command)
                command = str('self.settings["' + key_list[current_key]+ '"]= self.' + current_key)
                exec(command)
                print('==> XCE: found ' + key_list[current_key])
            except:
                print('==> XCE: WARNING: Failed to find settings for: ' + key_list[current_key] + ' Error type: '
                      + str(sys.exc_info()[0]))
                exec(str(current_key + " = ''"))
                continue


        try:
            pickled_settings = pickle.load(open(file_name, "rb"))
            if pickled_settings['beamline_directory'] != self.beamline_directory:
                self.beamline_directory = pickled_settings['beamline_directory']
                self.target_list, self.visit_list = XChemMain.get_target_and_visit_list(self.beamline_directory)
                self.settings['beamline_directory'] = self.beamline_directory
                self.populate_target_selection_combobox(self.target_selection_combobox)


            self.layout_funcs.pandda_html(self)
            self.show_pandda_html_summary()

            self.html_export_directory_label.setText(self.html_export_directory)

            self.group_deposition_directory_label.setText(self.group_deposit_directory)

            self.datasets_summary_file_label.setText(self.datasets_summary_file)

            self.data_source_file = pickled_settings['data_source']
            if self.data_source_file != '':
                self.settings['data_source'] = os.path.join(self.database_directory, self.data_source_file)
                # this is probably not necessary
                if os.path.isfile(self.settings['data_source']):
                    write_enabled = self.check_write_permissions_of_data_source()
                    if not write_enabled:
                        self.data_source_file_label.setText('')
                        self.data_source_set = False
                    else:
                        self.data_source_file_label.setText(
                            os.path.join(self.database_directory, self.data_source_file))
                        self.data_source_set = True
                        self.db = XChemDB.data_source(os.path.join(self.database_directory, self.data_source_file))
                        self.datasource_menu_reload_samples()

            reference_directory_temp = pickled_settings['reference_directory']
            if reference_directory_temp != self.reference_directory:
                self.reference_directory = reference_directory_temp
                self.settings['reference_directory'] = self.reference_directory
                self.update_reference_files(' ')
                for xtal in self.initial_model_dimple_dict:
                    reference_file_selection_combobox = self.initial_model_dimple_dict[xtal][1]
                    self.populate_reference_combobox(reference_file_selection_combobox)

            self.initial_model_directory_label.setText(self.initial_model_directory)
            self.panddas_directory_label.setText(self.panddas_directory)
            self.pandda_output_data_dir_entry.setText(self.panddas_directory)
            self.reference_directory_label.setText(self.reference_directory)
            self.beamline_directory_label.setText(self.beamline_directory)
            self.ccp4_scratch_directory_label.setText(self.ccp4_scratch_directory)
            self.reference_file_list = self.get_reference_file_list(' ')

            self.update_all_tables()

        except KeyError:
            self.update_status_bar('Sorry, this is not a XChemExplorer config file!')
            self.update_log.insert('Sorry, this is not a XChemExplorer config file!')

        except:
            print "Unexpected error:", sys.exc_info()[0]
            raise

    def save_config_file(self):
        file_name = str(QtGui.QFileDialog.getSaveFileName(self.window, 'Save file', self.current_directory))
        # make sure that the file always has .conf extension
        if str(file_name).rfind('.') != -1:
            file_name = file_name[:file_name.rfind('.')] + '.conf'
        else:
            file_name = file_name + '.conf'
        pickle.dump(self.settings, open(file_name, 'wb'))

    def update_reference_files(self, reference_root):
        self.reference_file_list = self.get_reference_file_list(reference_root)
        self.populate_reference_combobox(self.reference_file_selection_combobox)
        self.populate_reference_combobox(self.pandda_reference_file_selection_combobox)



    def check_status_rerun_dimple_on_all_autoprocessing_files(self):
        print('hallo')

    def rerun_dimple_on_all_autoprocessing_files(self):
        job_list = []
        self.update_log.insert('preparing to run DIMPLE on all autoprocessing files')
        for xtal in self.data_collection_dict:
            for entry in self.data_collection_dict[xtal]:
                if entry[0] == 'logfile':
                    db_dict = entry[6]
                    try:
                        if os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'],
                                                       db_dict['DataProcessingMTZfileName'])) or \
                                os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'])):
                            job_list = self.get_job_list_for_dimple_rerun(xtal, job_list, db_dict, entry)
                    except KeyError:
                        try:
                            if os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'])):
                                job_list = self.get_job_list_for_dimple_rerun(xtal, job_list, db_dict, entry)
                        except KeyError:
                            continue
        if job_list != []:
            self.update_log.insert('trying to run DIMPLE on ALL auto-processing files')
            self.check_before_running_dimple(job_list)

    def run_dimple_on_selected_autoprocessing_file(self):
        job_list = []
        for xtal in sorted(self.initial_model_dimple_dict):
            if self.initial_model_dimple_dict[xtal][0].isChecked():
                db_dict = self.xtal_db_dict[xtal]

                # the if statement below is so convoluted, so that it is compatible with older data source files

                if os.path.isfile(
                        os.path.join(db_dict['ProjectDirectory'], xtal, db_dict['DataProcessingPathToMTZfile'],
                                     db_dict['DataProcessingMTZfileName'])) or \
                        os.path.isfile(
                            os.path.join(db_dict['ProjectDirectory'], xtal, db_dict['DataProcessingPathToMTZfile'])) or \
                        os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'],
                                                    db_dict['DataProcessingMTZfileName'])) or \
                        os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'])):

                    if os.path.isfile(
                            os.path.join(db_dict['DataProcessingPathToMTZfile'], db_dict['DataProcessingMTZfileName'])):
                        mtzin = os.path.join(db_dict['DataProcessingPathToMTZfile'],
                                             db_dict['DataProcessingMTZfileName'])
                    elif os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'])):
                        mtzin = os.path.join(db_dict['DataProcessingPathToMTZfile'])
                    elif os.path.isfile(
                            os.path.join(db_dict['ProjectDirectory'], xtal, db_dict['DataProcessingPathToMTZfile'],
                                         db_dict['DataProcessingMTZfileName'])):
                        mtzin = os.path.join(db_dict['ProjectDirectory'], xtal, db_dict['DataProcessingPathToMTZfile'],
                                             db_dict['DataProcessingMTZfileName'])
                    elif os.path.isfile(
                            os.path.join(db_dict['ProjectDirectory'], xtal, db_dict['DataProcessingPathToMTZfile'])):
                        mtzin = os.path.join(db_dict['ProjectDirectory'], xtal, db_dict['DataProcessingPathToMTZfile'])

                    reference_file = str(self.initial_model_dimple_dict[xtal][1].currentText())

                    reference_file_pdb = os.path.join(self.reference_directory, reference_file + '.pdb')

                    if not os.path.isfile(reference_file_pdb):
                        continue

                    if os.path.isfile(os.path.join(self.reference_directory, reference_file + '.mtz')):
                        reference_file_mtz = ' -R ' + os.path.join(self.reference_directory, reference_file + '.mtz')
                    else:
                        reference_file_mtz = ''

                    if os.path.isfile(os.path.join(self.reference_directory, reference_file + '.cif')):
                        reference_file_cif = ' --libin ' + os.path.join(self.reference_directory,
                                                                        reference_file + '.cif')
                    else:
                        reference_file_cif = ''

                    job_list.append([xtal,
                                     'dimple_rerun_on_selected_file',
                                     mtzin,
                                     reference_file_pdb,
                                     reference_file_mtz,
                                     reference_file_cif])

        if job_list != []:
            self.update_log.insert('trying to run DIMPLE on SELECTED auto-processing files')
            self.check_before_running_dimple(job_list)

    def remove_selected_dimple_files(self):
        job_list = []
        for xtal in sorted(self.initial_model_dimple_dict):
            if self.initial_model_dimple_dict[xtal][0].isChecked():
                job_list.append(xtal)

        if job_list != []:
            msgBox = QtGui.QMessageBox()
            msgBox.setText("Do you really want to delete {0!s} Dimple files?".format(len(job_list)))
            msgBox.addButton(QtGui.QPushButton('Go'), QtGui.QMessageBox.YesRole)
            msgBox.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)
            reply = msgBox.exec_();

            if reply == 0:
                self.status_bar.showMessage('preparing to remove DIMPLE files')
                self.update_log.insert('preparing to remove DIMPLE files')
                self.work_thread = XChemThread.remove_selected_dimple_files(job_list,
                                                                            self.initial_model_directory,
                                                                            self.xce_logfile,
                                                                            self.database_directory,
                                                                            self.data_source_file)
                self.explorer_active = 1
                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                self.connect(self.work_thread, QtCore.SIGNAL("datasource_menu_reload_samples"),
                             self.datasource_menu_reload_samples)
                self.work_thread.start()

    def run_xia2_on_selected_datasets(self, overwrite):

        # check which programs should be run
        protocol = []
        if self.xia2_3d_checkbox.isChecked():
            protocol.append('3d')
        if self.xia2_3dii_checkbox.isChecked():
            protocol.append('3dii')
        if self.xia2_dials_checkbox.isChecked():
            protocol.append('dials')

        # space group
        spg = []
        if str(self.reprocess_space_group_comboxbox.currentText()) != 'ignore':
            spg.append(str(self.reprocess_space_group_comboxbox.currentText()))

        # reference file
        ref = []
        if os.path.isfile(self.diffraction_data_reference_mtz):
            ref.append(self.diffraction_data_reference_mtz)

        # resolution limit
        reso_limit = []
        if str(self.reprocess_isigma_combobox.currentText()) != 'default':
            reso_limit.append(str(self.reprocess_isigma_combobox.currentText()))

        # cc 1/2
        cc_half = []
        if str(self.reprocess_cc_half_combobox.currentText()) != 'default':
            cc_half.append(str(self.reprocess_cc_half_combobox.currentText()))

        run_dict = {}
        allRows = self.datasets_reprocess_table.rowCount()
        for row in xrange(0, allRows):
            dataset_id = str(self.datasets_reprocess_table.item(row, 0).text())
            sample_id = str(self.datasets_reprocess_table.item(row, 1).text())
            if self.diffraction_data_table_dict[dataset_id][0].isChecked():
                run_dict[sample_id] = self.diffraction_data_dict[dataset_id]

        if protocol != [] and run_dict != {}:
            self.work_thread = XChemProcess.run_xia2(self.initial_model_directory,
                                                     run_dict,
                                                     protocol,
                                                     spg,
                                                     ref,
                                                     reso_limit,
                                                     cc_half,
                                                     self.xce_logfile,
                                                     self.external_software,
                                                     self.ccp4_scratch_directory,
                                                     self.max_queue_jobs,
                                                     os.path.join(self.database_directory, self.data_source_file),
                                                     overwrite)
            self.explorer_active = 1
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.work_thread.start()
        else:
            self.update_log.insert('please select datasets and/ or data processing protocol')
            self.update_status_bar('please select datasets and/ or data processing protocol')

    def update_reprocessing_table(self):
        allRows = self.datasets_reprocess_table.rowCount()
        for row in xrange(0, allRows):
            sample_id = str(self.datasets_reprocess_table.item(row, 1).text())
            if sample_id in self.xtal_db_dict:
                db_dict = self.xtal_db_dict[sample_id]
                cell_text = QtGui.QTableWidgetItem()
                cell_text.setText(db_dict['DataProcessingStatus'])
                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                if db_dict['DataProcessingStatus'] == 'running':
                    cell_text.setBackground(QtGui.QColor(100, 230, 150))
                elif db_dict['DataProcessingStatus'] == 'pending':
                    cell_text.setBackground(QtGui.QColor(20, 100, 230))
                elif db_dict['DataProcessingStatus'] == 'started':
                    cell_text.setBackground(QtGui.QColor(230, 240, 110))
                elif db_dict['DataProcessingStatus'] == 'finished':
                    cell_text.setBackground(QtGui.QColor(255, 255, 255))
                self.datasets_reprocess_table.setItem(row, 7, cell_text)

    def get_job_list_for_dimple_rerun(self, xtal, job_list, db_dict, entry):
        self.status_bar.showMessage('checking: ' + str(
            os.path.join(db_dict['DataProcessingPathToMTZfile'], db_dict['DataProcessingMTZfileName'])))
        suitable_reference = []
        for reference in self.reference_file_list:
            # first we need one in the same pointgroup
            if reference[5] == db_dict['DataProcessingPointGroup']:
                try:
                    difference = math.fabs(1 - (float(db_dict['DataProcessingUnitCellVolume']) / float(reference[4])))
                    suitable_reference.append([reference[0], difference])
                except ValueError:
                    continue
        if suitable_reference != []:
            reference_file = min(suitable_reference, key=lambda x: x[1])[0]
            visit = entry[1]
            run = entry[2]
            autoproc = entry[4]

            reference_file_pdb = os.path.join(self.reference_directory, reference_file + '.pdb')

            if os.path.isfile(os.path.join(self.reference_directory, reference_file + '.mtz')):
                reference_file_mtz = ' -R ' + os.path.join(self.reference_directory, reference_file + '.mtz')
            else:
                reference_file_mtz = ''

            if os.path.isfile(os.path.join(self.reference_directory, reference_file + '.cif')):
                reference_file_cif = ' --libin ' + os.path.join(self.reference_directory, reference_file + '.cif')
            else:
                reference_file_cif = ''

            if os.path.isfile(
                    os.path.join(db_dict['DataProcessingPathToMTZfile'], db_dict['DataProcessingMTZfileName'])):
                mtzin = os.path.join(db_dict['DataProcessingPathToMTZfile'], db_dict['DataProcessingMTZfileName'])
            elif os.path.isfile(os.path.join(db_dict['DataProcessingPathToMTZfile'])):
                mtzin = os.path.join(db_dict['DataProcessingPathToMTZfile'])

            self.update_log.insert('adding ' + xtal + visit + '-' + run + autoproc + ' to list')
            job_list.append([xtal,
                             visit + '-' + run + autoproc,
                             mtzin,
                             reference_file_pdb,
                             reference_file_mtz,
                             reference_file_cif])
        self.status_bar.showMessage('idle')
        return job_list

    def check_before_running_dimple(self, job_list):

        msgBox = QtGui.QMessageBox()
        msgBox.setText(
            "Do you really want to run {0!s} Dimple jobs?\nNote: we will not run more than 100 at once on the cluster!".format(
                len(job_list)))
        msgBox.addButton(QtGui.QPushButton('Go'), QtGui.QMessageBox.YesRole)
        msgBox.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)
        reply = msgBox.exec_();

        if reply == 0:
            self.status_bar.showMessage('preparing {0!s} DIMPLE jobs'.format(len(job_list)))
            self.update_log.insert('preparing to run {0!s} DIMPLE jobs'.format(len(job_list)))
            if self.external_software['qsub_array']:
                self.update_log.insert('we will be running an ARRAY job on the DLS computer cluster')
                self.update_log.insert(
                    'please note that the maximum number of jobs that will be running at once is {0!s}'.format(
                        self.max_queue_jobs))
                self.update_log.insert(
                    'you can change this in the PREFERENCES menu, but be warned that to high a number might break the cluster!')
            self.update_log.insert('preparing input files for DIMPLE...')
            self.work_thread = XChemThread.run_dimple_on_all_autoprocessing_files(job_list,
                                                                                  self.initial_model_directory,
                                                                                  self.external_software,
                                                                                  self.ccp4_scratch_directory,
                                                                                  self.database_directory,
                                                                                  self.data_source_file,
                                                                                  self.max_queue_jobs,
                                                                                  self.xce_logfile,
                                                                                  self.using_remote_qsub_submission,
                                                                                  self.remote_qsub_submission)
            self.explorer_active = 1
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.connect(self.work_thread, QtCore.SIGNAL("datasource_menu_reload_samples"),
                         self.datasource_menu_reload_samples)
            self.work_thread.start()









    def open_csv_file_translate_datasetID_to_sampleID(self):
        file_name_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(self.window, 'Open file', self.current_directory,
                                                                    '*.csv')
        file_name = tuple(file_name_temp)[0]
        self.translate_datasetID_to_sampleID_csv_label.setText(file_name)
        self.translate_datasetID_to_sampleID_file = file_name



    def update_datasets_reprocess_table(self, data_dict):
        self.update_log.insert('updating reprocess datasets table')
        print('updating reprocess datasets table')
        self.diffraction_data_table_dict = {}
        self.diffraction_data_dict = data_dict

        self.diffraction_data_search_info = 'found ' + str(len(self.diffraction_data_dict)) + ' datasets'
        self.diffraction_data_search_label.setText(self.diffraction_data_search_info)
        self.update_log.insert(self.diffraction_data_search_info)
        self.datasource_menu_reload_samples()
        # update table
        column_name = self.db.translate_xce_column_list_to_sqlite(self.datasets_reprocess_columns)
        # set rows to 0
        self.datasets_reprocess_table.setRowCount(0)
        for entry in sorted(self.diffraction_data_dict):
            self.update_log.insert(str(self.diffraction_data_dict[entry]))
            if entry in self.xtal_db_dict:
                db_dict = self.xtal_db_dict[entry]
            else:
                db_dict = {}
            row = self.datasets_reprocess_table.rowCount()
            self.datasets_reprocess_table.insertRow(row)
            for column, header in enumerate(column_name):
                if header[0] == 'Dataset ID' or header[0] == 'Sample ID':
                    cell_text = QtGui.QTableWidgetItem()
                    cell_text.setText(str(entry))
                    cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                    self.datasets_reprocess_table.setItem(row, column, cell_text)
                elif header[0] == 'Run\nxia2':
                    run_xia2 = QtGui.QCheckBox()
                    run_xia2.toggle()
                    self.datasets_reprocess_table.setCellWidget(row, column, run_xia2)
                    run_xia2.setChecked(False)
                    self.diffraction_data_table_dict[entry] = [run_xia2]
                else:
                    cell_text = QtGui.QTableWidgetItem()
                    if db_dict != {}:
                        if header[0] == 'DataProcessing\nStatus':
                            if str(db_dict[header[1]]) == 'running':
                                cell_text.setBackground(QtGui.QColor(100, 230, 150))
                            elif str(db_dict[header[1]]) == 'pending':
                                cell_text.setBackground(QtGui.QColor(20, 100, 230))
                            elif str(db_dict[header[1]]) == 'started':
                                cell_text.setBackground(QtGui.QColor(230, 240, 110))
                            elif str(db_dict[header[1]]) == 'finished':
                                cell_text.setBackground(QtGui.QColor(255, 255, 255))
                        cell_text.setText(str(db_dict[header[1]]))
                    else:
                        cell_text.setText('')
                    cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                    self.datasets_reprocess_table.setItem(row, column, cell_text)

    def update_all_tables(self):
        self.update_log.insert('checking for new reference files')
        self.update_status_bar('checking for new reference files')
        self.reference_file_list = self.get_reference_file_list(' ')
        self.update_log.insert('updating Overview table')
        self.update_status_bar('updating Overview table')
        self.populate_and_update_datasource_table()
        self.update_log.insert('updating Maps table')
        self.update_status_bar('updating Maps table')
        self.create_maps_table()
        self.update_log.insert('updating PANDDA table')
        self.update_status_bar('updating PANDDA table')
        self.populate_pandda_analyse_input_table()
        self.update_log.insert('updating REFINEMENT table')
        self.update_status_bar('updating REFINEMENT table')
        self.populate_and_update_refinement_table()
        self.update_log.insert('updating REPROCESSING table')
        self.update_status_bar('updating REPROCESSING table')
        self.update_reprocessing_table()
        self.update_status_bar('idle')
        self.update_summary_plot()



    def change_allowed_unitcell_difference_percent(self, text):
        try:
            self.allowed_unitcell_difference_percent = int(text)
            self.settings['unitcell_difference'] = self.allowed_unitcell_difference_percent
            self.update_log.insert(
                'changing max allowed unit cell difference between reference and xtal to {0!s} percent'.format(
                    self.allowed_unitcell_difference_percent))
        except ValueError:
            if str(text).find('.') != -1:
                self.allowed_unitcell_difference_percent = int(str(text)[:str(text).find('.')])
                self.settings['unitcell_difference'] = self.allowed_unitcell_difference_percent
                self.update_log.insert(
                    'changing max allowed unit cell difference between reference and xtal to {0!s} percent'.format(
                        self.allowed_unitcell_difference_percent))
            else:
                pass

    def change_max_queue_jobs(self, text):
        try:
            self.max_queue_jobs = int(text)
            self.settings['max_queue_jobs'] = self.max_queue_jobs
            self.update_log.insert('changing max number of jobs running simultaneously on DLS cluster to {0!s}'.format(
                self.max_queue_jobs))
        except ValueError:
            if str(text).find('.') != -1:
                self.max_queue_jobs = int(str(text)[:str(text).find('.')])
                self.settings['max_queue_jobs'] = self.max_queue_jobs
                self.update_log.insert(
                    'changing max number of jobs running simultaneously on DLS cluster to {0!s}'.format(
                        self.max_queue_jobs))
            else:
                pass

    def change_acceptable_low_resolution_limit(self, text):
        try:
            self.acceptable_low_resolution_limit_for_data = float(text)
            self.settings['too_low_resolution_data'] = self.acceptable_low_resolution_limit_for_data
        except ValueError:
            pass

    def change_filename_root(self, text):
        self.filename_root = str(text)
        self.settings['filename_root'] = self.filename_root

    def button_clicked(self):
        if self.data_source_set == False:
            print('sender text bit')
            if self.sender().text() == "Create New Data\nSource (SQLite)":
                file_name = str(QtGui.QFileDialog.getSaveFileName(self.window, 'Save file', self.database_directory))
                # make sure that the file always has .sqlite extension
                if file_name.rfind('.') != -1:
                    file_name = file_name[:file_name.rfind('.')] + '.sqlite'
                else:
                    file_name = file_name + '.sqlite'
                self.db = XChemDB.data_source(file_name)
                print('==> XCE: creating new data source')
                self.db.create_empty_data_source_file()
                self.db.create_missing_columns()
                if self.data_source_file == '':
                    self.database_directory = file_name[:file_name.rfind('/')]
                    self.data_source_file = file_name[file_name.rfind('/') + 1:]
                    self.data_source_file_label.setText(os.path.join(self.database_directory, self.data_source_file))
                    self.settings['database_directory'] = self.database_directory
                    self.settings['data_source'] = self.data_source_file
                    self.data_source_set = True
            else:
                self.no_data_source_selected()
                print('No datasource selected')
                pass

        # first find out which of the 'Run' or 'Status' buttons is sending
        for item in self.workflow_widget_dict:
            for widget in self.workflow_widget_dict[item]:
                if widget == self.sender():
                    # get index of item in self.workflow; Note this index should be the same as the index
                    # of the self.main_tab_widget which belongs to this task
                    task_index = self.workflow.index(item)
                    instruction = str(self.workflow_widget_dict[item][0].currentText())
                    action = str(self.sender().text())
                    if self.main_tab_widget.currentIndex() == task_index:
                        if self.explorer_active == 0 and self.data_source_set == True:
                            if action == 'Run':
                                print('==> XCE: Remote submission status = ' + str(self.using_remote_qsub_submission))
                                self.prepare_and_run_task(instruction)
                            elif action == 'Status':
                                self.get_status_of_workflow_milestone(instruction)
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
                    else:
                        self.need_to_switch_main_tab(task_index)

    def get_status_of_workflow_milestone(self, instruction):
        # first update all tables
        self.datasource_menu_reload_samples()

        cluster_dict = XChemMain.get_jobs_running_on_cluster()

        self.update_log.insert('getting status updates...')

        self.status_bar.showMessage('please check terminal window for further information')

        self.update_log.insert('{0!s} samples are currently in database'.format(str(len(self.xtal_db_dict))))

        if 'DIMPLE' in instruction:
            XChemMain.print_cluster_status_message('dimple', cluster_dict, self.xce_logfile)

        elif 'Create CIF/PDB/PNG file' in instruction:
            XChemMain.print_acedrg_status(self.xce_logfile, self.xtal_db_dict)
            XChemMain.print_cluster_status_message('acedrg', cluster_dict, self.xce_logfile)

        elif instruction.startswith('Run xia2 on selected datasets'):
            XChemMain.print_cluster_status_message('xia2', cluster_dict, self.xce_logfile)

        elif 'pandda' in instruction.lower():
            XChemMain.print_cluster_status_message('pandda', cluster_dict, self.xce_logfile)

        elif 'coot' in instruction.lower():
            XChemMain.print_cluster_status_message('refmac', cluster_dict, self.xce_logfile)

    def prepare_and_run_task(self, instruction):

        if instruction == 'Get New Results from Autoprocessing':
            self.rescore = False
            self.check_for_new_autoprocessing_results()

        elif instruction == 'Rescore Datasets':
            self.rescore = True
            self.select_best_autoprocessing_result()

        elif instruction == 'Run xia2 on selected datasets':
            self.run_xia2_on_selected_datasets(False)

        elif instruction == 'Run xia2 on selected datasets - overwrite':
            self.run_xia2_on_selected_datasets(True)

        elif instruction == 'Run DIMPLE on All Autoprocessing MTZ files':
            self.rerun_dimple_on_all_autoprocessing_files()

        elif instruction == 'Run DIMPLE on selected MTZ files':
            self.run_dimple_on_selected_autoprocessing_file()

        elif instruction == 'Remove selected DIMPLE PDB/MTZ files':
            self.remove_selected_dimple_files()

        elif instruction == 'Create CIF/PDB/PNG file of ALL compounds':
            self.create_cif_pdb_png_files('ALL')

        elif instruction == 'Create CIF/PDB/PNG file of NEW compounds':
            self.create_cif_pdb_png_files('NEW')

        elif instruction == 'Create CIF/PDB/PNG file of SELECTED compounds':
            self.create_cif_pdb_png_files('SELECTED')

        elif instruction == 'pandda.analyse':
            self.run_pandda_analyse('production_run')

        elif instruction == 'pre-run for ground state model':
            self.run_pandda_analyse('pre_run')

        elif instruction == 'pandda.inspect':
            self.run_pandda_inspect()

        elif instruction == 'run pandda.inspect at home':
            self.run_pandda_inspect_at_home()

        elif instruction == 'Export NEW PANDDA models':
            update_datasource_only = False
            which_models = 'new'
            self.run_pandda_export(update_datasource_only, which_models)

        elif instruction == 'Export ALL PANDDA models':
            update_datasource_only = False
            which_models = 'all'
            self.run_pandda_export(update_datasource_only, which_models)

        elif instruction == 'Export SELECTED PANDDA models':
            update_datasource_only = False
            which_models = 'selected'
            self.run_pandda_export(update_datasource_only, which_models)

        elif instruction == 'cluster datasets':
            self.cluster_datasets_for_pandda()

        elif instruction == 'Update datasource with results from pandda.inspect':
            update_datasource_only = True
            which_models = 'all'
            self.run_pandda_export(update_datasource_only, which_models)

        elif instruction == 'Show HTML summary':
            self.show_pandda_html_summary()

        elif instruction == 'Event Map -> SF':
            self.convert_event_maps_to_SF()

        elif instruction == 'check modelled ligands':
            self.compare_modelled_ligands_and_panddaTable()

        elif instruction.startswith("Open COOT") or instruction == 'Build ground state model':
            if not self.coot_running:
                self.update_log.insert('starting coot...')
                if instruction == "Open COOT":
                    interface = 'new'
                elif instruction == "Open COOT for old PanDDA":
                    interface = 'panddaV1'
                elif instruction == 'Build ground state model':
                    interface = 'reference'
                else:
                    interface = 'old'
                self.work_thread = XChemThread.start_COOT(self.settings, interface)
                self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
                self.work_thread.start()


        elif instruction == 'Update Deposition Table':
            self.update_deposition_table()




    def run_pandda_analyse(self, run):
        pandda_params = {
            'data_dir': str(self.pandda_input_data_dir_entry.text()),
            'out_dir': str(self.pandda_output_data_dir_entry.text()),
            'submit_mode': str(self.pandda_submission_mode_selection_combobox.currentText()),
            'nproc': str(self.pandda_nproc_entry.text()),
            'min_build_datasets': str(self.pandda_min_build_dataset_entry.text()),
            'pdb_style': str(self.pandda_pdb_style_entry.text()),
            'mtz_style': str(self.pandda_mtz_style_entry.text()),
            'sort_event': str(self.pandda_sort_event_combobox.currentText()),
            'max_new_datasets': str(self.pandda_max_new_datasets_entry.text()),
            'grid_spacing': str(self.pandda_grid_spacing_entry.text()),
            'pandda_dir_structure': str(self.pandda_input_data_dir_entry.text()),
            'perform_diffraction_data_scaling': str(self.wilson_checkbox.isChecked()),
            'filter_pdb': str(self.pandda_reference_file_selection_combobox.currentText()),
            'reference_dir': self.reference_directory,
            'appendix': '',
            'N_datasets': len(glob.glob(os.path.join(self.initial_model_directory, '*', 'dimple.pdb'))),
            'write_mean_map': 'interesting',
            'pandda_table': self.pandda_analyse_data_table
        }

        if run == 'pre_run':
            msgBox = QtGui.QMessageBox()
            msgBoxLayout = msgBox.layout()
            vbox = QtGui.QVBoxLayout()
            vbox.addWidget(QtGui.QLabel(XChemToolTips.pandda_pre_run(self.reference_directory)))
            hbox = QtGui.QHBoxLayout()
            hbox.addWidget(QtGui.QLabel('appendix:'))
            appendix = QtGui.QLineEdit()
            appendix.setText('pre')
            appendix.setFixedWidth(200)
            hbox.addWidget(appendix)
            vbox.addLayout(hbox)

            msgBoxLayout.addLayout(vbox, 0, 0)
            msgBox.addButton(QtGui.QPushButton('Go'), QtGui.QMessageBox.YesRole)
            msgBox.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)
            reply = msgBox.exec_();
            if reply == 0:
                pandda_params['appendix'] = str(appendix.text())
                pandda_params['max_new_datasets'] = '100'
                pandda_params['N_datasets'] = 100
                pandda_params['write_mean_map'] = 'all'
            else:
                return None

        self.update_log.insert('preparing pandda.analyse input script')
        self.work_thread = XChemPANDDA.run_pandda_analyse(pandda_params, self.xce_logfile,
                                                          os.path.join(self.database_directory, self.data_source_file))
        self.work_thread = XChemPANDDA.run_pandda_analyse(pandda_params, self.xce_logfile,
                                                          os.path.join(self.database_directory, self.data_source_file))
        #self.connect(self.work_thread, QtCore.SIGNAL("datasource_menu_reload_samples"),
                     #self.datasource_menu_reload_samples)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def cluster_datasets_for_pandda(self):

        pandda_params = {
            'out_dir': str(self.pandda_output_data_dir_entry.text()),
            'pdb_style': str(self.pandda_pdb_style_entry.text()),
            'mtz_style': str(self.pandda_mtz_style_entry.text())
        }
        self.update_log.insert('starting giant.cluster_mtzs_and_pdbs')
        self.work_thread = XChemPANDDA.giant_cluster_datasets(self.initial_model_directory, pandda_params,
                                                              self.xce_logfile, os.path.join(self.database_directory,
                                                                                             self.data_source_file),
                                                              run_pandda_analyse)
        self.explorer_active = 1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("datasource_menu_reload_samples"),
                     self.datasource_menu_reload_samples)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def run_pandda_inspect(self):
        self.settings['panddas_directory'] = str(self.pandda_output_data_dir_entry.text())
        print('==> XCE: starting pandda.inspect')
        self.work_thread = XChemThread.start_pandda_inspect(self.settings, self.xce_logfile)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def run_pandda_inspect_at_home(self):
        self.work_thread = XChemPANDDA.run_pandda_inspect_at_home(self.panddas_directory, self.xce_logfile)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)

    def convert_event_maps_to_SF(self):
        self.update_log.insert('converting all event maps in {0!s} to mtz files'.format(self.initial_model_directory))
        self.work_thread = XChemPANDDA.convert_all_event_maps_in_database(self.initial_model_directory,
                                                                          self.xce_logfile,
                                                                          os.path.join(self.database_directory,
                                                                                       self.data_source_file))
        self.explorer_active = 1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.work_thread.start()

    def compare_modelled_ligands_and_panddaTable(self):
        self.update_log.insert('checking agreement of ligands in refine.pdb and entries in panddaTable')
        self.work_thread = XChemPANDDA.check_number_of_modelled_ligands(self.initial_model_directory,
                                                                        self.xce_logfile,
                                                                        os.path.join(self.database_directory,
                                                                                     self.data_source_file))
        self.explorer_active = 1
        self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
        self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
        self.connect(self.work_thread, QtCore.SIGNAL("show_error_dict"), self.show_error_dict)
        self.work_thread.start()

    def run_pandda_export(self, update_datasource_only, which_models):
        modelList=[]
        self.settings['panddas_directory'] = str(self.pandda_output_data_dir_entry.text())
        if update_datasource_only:
            self.update_log.insert('updating data source with results from pandda.inspect')
        else:
            self.update_log.insert(
                'exporting PANDDA models, updating data source and launching inital refinement for new models')

        start_thread = False
        if which_models == 'all':
            self.update_log.insert('exporting ALL models! *** WARNING *** This may overwrite previous refinements!!!')
            msgBox = QtGui.QMessageBox()
            msgBox.setText("*** WARNING ***\nThis will overwrite all your manual selections!\nDo you want to continue?")
            msgBox.addButton(QtGui.QPushButton('Yes'), QtGui.QMessageBox.YesRole)
            msgBox.addButton(QtGui.QPushButton('No'), QtGui.QMessageBox.RejectRole)
            reply = msgBox.exec_();
            if reply == 0:
                if update_datasource_only:
                    self.update_log.insert('will update panddaTable in database only')
                else:
                    self.update_log.insert('will export ALL models!')
                start_thread = True
            else:
                start_thread = False

        elif which_models == 'selected':
            self.update_log.insert('exporting selected models -> please enter comma separated samplesIDs')
            msgBox = QtGui.QMessageBox()
            msgBoxLayout = msgBox.layout()
            msgBox.setText(
                "exporting selected models\nplease enter comma separated samplesIDs\ne.g. test-x001,test-x002...")
            vbox = QtGui.QVBoxLayout()
            sampleEntry = QtGui.QLineEdit()
            sampleEntry.setFixedWidth(200)
            vbox.addWidget(sampleEntry)
            msgBoxLayout.addLayout(vbox, 0, 0)
            msgBox.addButton(QtGui.QPushButton('Go'), QtGui.QMessageBox.YesRole)
            msgBox.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)
            reply = msgBox.exec_();
            if reply == 0:
                for item in str(sampleEntry.text()).split(','):
                    self.update_log.insert('trying to export: '+item.replace(' ',''))
                    modelList.append(item.replace(' ',''))
                start_thread = True
            else:
                start_thread = False



        else:
            self.update_log.insert('exporting new models only')
            start_thread = True

        if start_thread:
            self.work_thread = XChemPANDDA.run_pandda_export(self.panddas_directory,
                                                             os.path.join(self.database_directory,
                                                                          self.data_source_file),
                                                             self.initial_model_directory, self.xce_logfile,
                                                             update_datasource_only, which_models,modelList)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.work_thread.start()

    def show_pandda_html_summary(self):
        self.pandda_initial_html.load(QtCore.QUrl(self.pandda_initial_html_file))
        self.pandda_initial_html.show()
        self.pandda_analyse_html.load(QtCore.QUrl(self.pandda_analyse_html_file))
        self.pandda_analyse_html.show()
        self.add_map_html()
        self.pandda_inspect_html.load(QtCore.QUrl(self.pandda_inspect_html_file))
        self.pandda_inspect_html.show()

    def create_cif_pdb_png_files(self, todo):
        tmp = self.db.execute_statement(
            "select CrystalName,CompoundCode,CompoundSmiles from mainTable where CrystalName is not '' and CompoundSmiles is not '' and CompoundSmiles is not NULL;")
        compound_list = []
        for item in tmp:
            if str(item[1]) == '' or str(item[1]) == 'NULL':
                compoundID = 'compound'
            else:
                compoundID = str(item[1])

            if todo == 'ALL':
                compound_list.append([str(item[0]), compoundID, str(item[2])])
            elif todo == 'NEW':
                if not os.path.isfile(os.path.join(self.initial_model_directory, str(item[0]), compoundID + '.cif')):
                    compound_list.append([str(item[0]), compoundID, str(item[2])])
            elif todo == 'SELECTED':
                if str(item[0]) in self.initial_model_dimple_dict:
                    if self.initial_model_dimple_dict[str(item[0])][0].isChecked():
                        compound_list.append([str(item[0]), compoundID, str(item[2])])

        if compound_list != []:
            self.update_log.insert(
                'trying to create cif and pdb files for ' + str(len(compound_list)) + ' compounds using ACEDRG...')
            if self.external_software['qsub']:
                self.update_log.insert(
                    'will try sending ' + str(len(compound_list)) + ' jobs to your computer cluster!')
            elif self.external_software['qsub_array']:
                self.update_log.insert('will try sending ' + str(
                    len(compound_list)) + ' jobs as part of an ARRAY job to your computer cluster!')
            else:
                self.update_log.insert('apparently no cluster available, so will run ' + str(
                    len(compound_list)) + ' sequential jobs on one core of your local machine.')
                self.update_log.insert('this could take a while...')
            self.explorer_active = 1
            self.work_thread = XChemThread.create_png_and_cif_of_compound(self.external_software,
                                                                          self.initial_model_directory,
                                                                          compound_list,
                                                                          self.database_directory,
                                                                          self.data_source_file,
                                                                          todo,
                                                                          self.ccp4_scratch_directory,
                                                                          self.xce_logfile,
                                                                          self.max_queue_jobs,
                                                                          self.restraints_program)
            self.connect(self.work_thread, QtCore.SIGNAL("update_progress_bar"), self.update_progress_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("update_status_bar(QString)"), self.update_status_bar)
            self.connect(self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished)
            self.connect(self.work_thread, QtCore.SIGNAL("datasource_menu_reload_samples"),
                         self.datasource_menu_reload_samples)
            self.work_thread.start()

    def update_deposition_table(self):
        # check if PanDDA models are ready for deposition

        depositChecks = XChemDeposit.update_deposition_table(
            os.path.join(self.database_directory, self.data_source_file))

        toDeposit, mismatch = depositChecks.PanDDA_models_to_deposit()

        if mismatch != {}:
            self.update_log.insert('The following samples contain ligand that are not ready for deposition:')
            for entry in mismatch:
                self.update_log.insert(entry[0] + ' -> site: ' + entry[1] + ' @ ' + entry[2] + ' => ' + entry[4])
            self.update_log.insert('You need to change this before you can continue!')
            return None

        for xtal in toDeposit:
            self.db.update_insert_depositTable(xtal, {})

    def show_html_summary_and_diffraction_image(self):
        for key in self.albula_button_dict:
            if self.albula_button_dict[key][0] == self.sender():
                print('==> XCE: showing html summary in firefox')
                self.show_html_summary_in_firefox(key)

    def need_to_switch_main_tab(self, task_index):
        msgBox = QtGui.QMessageBox()
        msgBox.setText("Need to switch main tab before you can launch this job")
        msgBox.addButton(QtGui.QPushButton('Yes'), QtGui.QMessageBox.YesRole)
        msgBox.addButton(QtGui.QPushButton('No'), QtGui.QMessageBox.RejectRole)
        reply = msgBox.exec_();
        if reply == 0:
            self.main_tab_widget.setCurrentIndex(task_index)

    def check_write_permissions_of_data_source(self):
        write_enabled = True
        if not os.access(os.path.join(self.database_directory, self.data_source_file), os.W_OK):
            QtGui.QMessageBox.warning(self.window, "Data Source Problem",
                                      '\nData Source is Read-Only\n',
                                      QtGui.QMessageBox.Cancel, QtGui.QMessageBox.NoButton,
                                      QtGui.QMessageBox.NoButton)
            write_enabled = False
        return write_enabled

    def no_data_source_selected(self):
        QtGui.QMessageBox.warning(self.window, "Data Source Problem",
                                  ('Please set or create a data source file\n') +
                                  ('Options:\n') +
                                  ('1. Use an existing file:\n') +
                                  ('- Settings -> Select Data Source File\n') +
                                  ('2. Create a new file\n') +
                                  ('- Data Source -> Create New Data\nSource (SQLite)'),
                                  QtGui.QMessageBox.Cancel, QtGui.QMessageBox.NoButton,
                                  QtGui.QMessageBox.NoButton)

    def update_progress_bar(self, progress):
        self.progress_bar.setValue(progress)

    def update_status_bar(self, message):
        self.status_bar.showMessage(message)

    def thread_finished(self):
        self.explorer_active = 0
        self.update_progress_bar(0)
        self.update_status_bar('idle')

    def show_error_dict(self, errorDict):
        text = ''
        for key in errorDict:
            text += '{0!s}:\n'.format(key)
            for entry in errorDict[key]:
                text += '  - ' + entry + '\n'
        msgBox = QtGui.QMessageBox()
        msgBox.setText(text)
        msgBox.exec_()


    def find_suitable_reference_file(self, db_dict):
        reference_file = []
        dummy = ['...', '', '', '', 0, '0']
        reference_file.append([dummy, 999])
        suitable_reference = []
        for reference in self.reference_file_list:
            # first we need one in the same pointgroup
            if reference[5] == db_dict['DataProcessingPointGroup']:
                try:
                    difference = math.fabs(
                        1 - (float(db_dict['DataProcessingUnitCellVolume']) / float(reference[4]))) * 100
                    reference_file.append([reference, difference])
                except ValueError:
                    continue
        return reference_file

    def create_maps_table(self):
        column_name = self.db.translate_xce_column_list_to_sqlite(self.maps_table_columns)

        for xtal in sorted(self.xtal_db_dict):
            new_xtal = False
            db_dict = self.xtal_db_dict[xtal]
            if str(db_dict['DataCollectionOutcome']).lower().startswith('success'):
                reference_file = self.find_suitable_reference_file(db_dict)
                smallest_uc_difference = min(reference_file, key=lambda x: x[1])
                row = self.maps_table.rowCount()
                if xtal not in self.initial_model_dimple_dict:
                    self.maps_table.insertRow(row)
                    current_row = row
                    new_xtal = True
                else:
                    for table_row in range(row):
                        if self.maps_table.item(table_row, 0).text() == xtal:
                            current_row = table_row
                            break
                for column, header in enumerate(column_name):
                    if header[0] == 'Sample ID':
                        cell_text = QtGui.QTableWidgetItem()
                        cell_text.setText(str(xtal))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.maps_table.setItem(current_row, column, cell_text)
                    elif header[0] == 'Select':
                        if new_xtal:
                            run_dimple = QtGui.QCheckBox()
                            run_dimple.toggle()
                            self.maps_table.setCellWidget(current_row, column, run_dimple)
                            run_dimple.setChecked(False)
                    elif header[0] == 'Reference\nSpaceGroup':
                        cell_text = QtGui.QTableWidgetItem()
                        cell_text.setText(str(smallest_uc_difference[0][1]))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.maps_table.setItem(current_row, column, cell_text)
                    elif header[0] == 'Difference\nUC Volume (%)':
                        cell_text = QtGui.QTableWidgetItem()
                        smallest_uc_difference = min(reference_file, key=lambda x: x[1])
                        cell_text.setText(str(round(float(smallest_uc_difference[1]), 1)))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.maps_table.setItem(current_row, column, cell_text)
                    elif header[0] == 'Reference File':
                        if new_xtal:
                            reference_file_selection_combobox = QtGui.QComboBox()
                            self.populate_reference_combobox(reference_file_selection_combobox)
                            if float(smallest_uc_difference[1]) < self.allowed_unitcell_difference_percent:
                                index = reference_file_selection_combobox.findText(str(smallest_uc_difference[0][0]),
                                                                                   QtCore.Qt.MatchFixedString)
                                reference_file_selection_combobox.setCurrentIndex(index)
                            else:
                                reference_file_selection_combobox.setCurrentIndex(0)
                            self.maps_table.setCellWidget(current_row, column,
                                                          reference_file_selection_combobox)
                        else:
                            reference_file_selection_combobox = self.initial_model_dimple_dict[xtal][1]
                            self.populate_reference_combobox(reference_file_selection_combobox)
                            if float(smallest_uc_difference[1]) < self.allowed_unitcell_difference_percent:
                                index = reference_file_selection_combobox.findText(str(smallest_uc_difference[0][0]),
                                                                                   QtCore.Qt.MatchFixedString)
                                reference_file_selection_combobox.setCurrentIndex(index)
                            else:
                                reference_file_selection_combobox.setCurrentIndex(0)
                    else:
                        cell_text = QtGui.QTableWidgetItem()
                        cell_text.setText(str(db_dict[header[1]]))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        if header[0] == 'Dimple\nStatus':
                            if str(db_dict[header[1]]) == 'running':
                                cell_text.setBackground(QtGui.QColor(100, 230, 150))
                            elif str(db_dict[header[1]]) == 'pending':
                                cell_text.setBackground(QtGui.QColor(20, 100, 230))
                            elif str(db_dict[header[1]]) == 'started':
                                cell_text.setBackground(QtGui.QColor(230, 240, 110))
                            elif str(db_dict[header[1]]) == 'finished':
                                cell_text.setBackground(QtGui.QColor(255, 255, 255))
                        if header[0] == 'Compound\nStatus':
                            if str(db_dict[header[1]]) == 'running':
                                cell_text.setBackground(QtGui.QColor(100, 230, 150))
                            elif str(db_dict[header[1]]) == 'pending':
                                cell_text.setBackground(QtGui.QColor(20, 100, 230))
                            elif str(db_dict[header[1]]) == 'started':
                                cell_text.setBackground(QtGui.QColor(230, 240, 110))
                            elif str(db_dict[header[1]]) == 'restraints generated':
                                cell_text.setBackground(QtGui.QColor(255, 255, 255))
                            elif str(db_dict[header[1]]) == 'restraints failed':
                                cell_text.setBackground(QtGui.QColor(255, 0, 0))
                            elif str(db_dict[header[1]]) == 'missing smiles':
                                cell_text.setBackground(QtGui.QColor(240, 150, 20))
                        self.maps_table.setItem(current_row, column, cell_text)
            if new_xtal:
                self.initial_model_dimple_dict[xtal] = [run_dimple, reference_file_selection_combobox]

    def preferences_data_to_copy_combobox_changed(self, i):
        text = str(self.preferences_data_to_copy_combobox.currentText())
        for item in self.preferences_data_to_copy:
            if item[0] == text:
                self.preferences['processed_data_to_copy'] = item[1]
                break

    def preferences_selection_mechanism_combobox_changed(self, i):
        text = str(self.preferences_selection_mechanism_combobox.currentText())
        self.preferences['dataset_selection_mechanism'] = text

    def preferences_restraints_generation_combobox_changed(self):
        text = str(self.preferences_restraints_generation_combobox.currentText())
        self.restraints_program = text
        self.update_log.insert('will use {0!s} for generation of ligand coordinates and restraints'.format(text))

    def refinement_outcome_combobox_changed(self):
        for xtal in self.refinement_table_dict:
            if self.sender() == self.refinement_table_dict[xtal]:
                db_dict = {}
                db_dict['RefinementOutcome'] = str(self.sender().currentText())
                self.db.create_or_remove_missing_records_in_depositTable(self.xce_logfile, xtal, 'ligand_bound',
                                                                         db_dict)

    def get_reference_file_list(self, reference_root):
        # check available reference files
        reference_file_list = []
        dummy = ['...', '', '', '', 0, '0']
        reference_file_list.append(dummy)
        if os.path.isfile(os.path.join(self.reference_directory, reference_root + '.pdb')):
            pdb_reference = parse().PDBheader(os.path.join(self.reference_directory, reference_root + '.pdb'))
            spg_reference = pdb_reference['SpaceGroup']
            unitcell_reference = pdb_reference['UnitCell']
            lattice_reference = pdb_reference['Lattice']
            unitcell_volume_reference = pdb_reference['UnitCellVolume']
            pointgroup_reference = pdb_reference['PointGroup']
            reference_file_list.append([reference_root,
                                        spg_reference,
                                        unitcell_reference,
                                        lattice_reference,
                                        unitcell_volume_reference,
                                        pointgroup_reference])
        else:
            for files in glob.glob(self.reference_directory + '/*'):
                if files.endswith('.pdb'):
                    reference_root = files[files.rfind('/') + 1:files.rfind('.')]
                    if os.path.isfile(os.path.join(self.reference_directory, reference_root + '.pdb')):
                        pdb_reference = parse().PDBheader(
                            os.path.join(self.reference_directory, reference_root + '.pdb'))
                        spg_reference = pdb_reference['SpaceGroup']
                        unitcell_reference = pdb_reference['UnitCell']
                        lattice_reference = pdb_reference['Lattice']
                        unitcell_volume_reference = pdb_reference['UnitCellVolume']
                        pointgroup_reference = pdb_reference['PointGroup']
                        reference_file_list.append([reference_root,
                                                    spg_reference,
                                                    unitcell_reference,
                                                    lattice_reference,
                                                    unitcell_volume_reference,
                                                    pointgroup_reference])
        for n, file in enumerate(reference_file_list):
            self.update_log.insert('reference file {0!s}: {1!s}'.format(n, file))
        return reference_file_list

    def dataset_outcome_combobox_change_outcome(self, text):
        outcome = str(text)
        xtal = ''
        for key in self.dataset_outcome_combobox_dict:
            if self.dataset_outcome_combobox_dict[key] == self.sender():
                xtal = key
                self.update_log.insert('user changed data collection outcome of {0!s} to {1!s}'.format(xtal, outcome))
                self.update_log.insert('updating dataset outcome in datasource for {0!s}'.format(xtal))
                update_dict = {}
                update_dict['DataCollectionOutcome'] = outcome
                self.db.update_insert_data_source(xtal, update_dict)
                break

#
    def set_run_dimple_flag(self, state):
        if state == QtCore.Qt.Checked:
            for key in self.initial_model_dimple_dict:
                self.initial_model_dimple_dict[key][0].setChecked(True)
        else:
            for key in self.initial_model_dimple_dict:
                self.initial_model_dimple_dict[key][0].setChecked(False)

    def get_sample_list_from_table(self,table):
        sampleList = []
        allRows = table.rowCount()
        for row in xrange(0, allRows):
            sample_id = str(table.item(row, 0).text())
            sampleList.append(sample_id)
        return sorted(sampleList)

    def get_row_of_sample_in_table(self,table,xtal):
        allRows = table.rowCount()
        sampleRow = allRows
        for n,row in enumerate(xrange(0, allRows)):
            sample_id = str(table.item(row, 0).text())
            if sample_id == xtal:
                sampleRow = n
                break
        return sampleRow


    def update_row_in_table(self,sample,row,db_dict,table,columns_to_show):
        xtal = str(sample)
        column_name = self.db.translate_xce_column_list_to_sqlite(columns_to_show)

        for column, header in enumerate(column_name):

            if header[0] == 'Sample ID':
                cell_text = QtGui.QTableWidgetItem()
                cell_text.setText(str(xtal))
                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                table.setItem(row, column, cell_text)

            elif header[0] == 'DataCollection\nOutcome':
                if xtal not in self.dataset_outcome_combobox_dict:
                    dataset_outcome_combobox = QtGui.QComboBox()
                    for outcomeItem in self.dataset_outcome:
                        dataset_outcome_combobox.addItem(outcomeItem)
                    dataset_outcome_combobox.activated[str].connect(self.dataset_outcome_combobox_change_outcome)
                    self.dataset_outcome_combobox_dict[xtal] = dataset_outcome_combobox
                    table.setCellWidget(row, column, dataset_outcome_combobox)
                index = self.dataset_outcome_combobox_dict[xtal].findText(str(db_dict['DataCollectionOutcome']), QtCore.Qt.MatchFixedString)
                self.dataset_outcome_combobox_dict[xtal].setCurrentIndex(index)

            elif header[0].startswith('img'):
                if os.path.isfile(db_dict[header[1]]):
                    pixmap = QtGui.QPixmap(db_dict[header[1]])
                else:
                    pixmap = QtGui.QPixmap(
                        os.path.join(os.getenv('XChemExplorer_DIR'), 'image', 'IMAGE_NOT_AVAILABLE.png'))
                image = QtGui.QLabel()
                image.resize(128, 80)
                image.setPixmap(pixmap.scaled(image.size(), QtCore.Qt.KeepAspectRatio))
                table.setCellWidget(row, column, image)

            elif header[0] == 'Select':
                checkbox = QtGui.QCheckBox()
                checkbox.toggle()
                if table == self.deposition_table_apo:
                    if xtal not in self.deposition_table_apo_dict:
                        self.deposition_table_apo_dict[xtal] = checkbox
                if table == self.deposition_table_bound:
                    if xtal not in self.deposition_table_bound_dict:
                        self.deposition_table_bound_dict[xtal] = checkbox
                table.setCellWidget(row, column, checkbox)
                checkbox.setChecked(False)

            #elif header[0].startswith('SoakDB\nBarcode') or header[0].startswith('GDA\nBarcode'):
                #                        if new_xtal:
                #                            cell_text = QtGui.QTableWidgetItem()
                #                            if xtal in pinDict:
                #                                if header[0].startswith('SoakDB\nBarcode'):
                #                                    cell_text.setText(str(pinDict[xtal][0]))
                #                                elif header[0].startswith('GDA\nBarcode'):
                #                                    cell_text.setText(str(pinDict[xtal][1]))
                #                                if pinDict[xtal][0] == 'NULL' or pinDict[xtal][1] == 'NULL':
                #                                    cell_text.setBackground(QtGui.QColor(255, 215, 0))
                #                                elif pinDict[xtal][0] != pinDict[xtal][1]:
                #                                    cell_text.setBackground(QtGui.QColor(255, 0, 0))
                #                            else:
                #                                cell_text.setText('')
                #                            cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                #                            self.datasets_summary_table.setItem(current_row, column, cell_text)
            else:
                cell_text = QtGui.QTableWidgetItem()
                # in case data collection failed for whatever reason
                try:
                    cell_text.setText(str(db_dict[header[1]]))
                except KeyError:  # older pkl files may not have all the columns
                    cell_text.setText('n/a')
                    #                        else:
                    #                            if header[0].startswith('Resolution\n[Mn<I/sig(I)> = 1.5]'):
                    #                                cell_text.setText('999')
                    #                            elif header[0].startswith('DataProcessing\nRfree'):
                    #                                cell_text.setText('999')
                    #                            elif header[0].startswith('Rmerge\nLow'):
                    #                                cell_text.setText('999')
                    #                            else:
                    #                                cell_text.setText('')
                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                table.setItem(row, column, cell_text)

    def populate_datasets_summary_table(self):
        self.status_bar.showMessage(
            'Building summary table for data processing results; be patient this may take a while')

        # get information about all samples collected during the current visit
        visit, beamline = XChemMain.getVisitAndBeamline(self.beamline_directory)
        self.update_log.insert('reading information about collected crystals from database...')
        collectedXtalsDict = self.db.xtals_collected_during_visit_as_dict(visit)

        # instead of using dictionaries, query table of which crystals are in table
        samples_in_table = self.get_sample_list_from_table(self.datasets_summary_table)
        for xtal in sorted(collectedXtalsDict):
            if xtal not in samples_in_table:
                row = self.datasets_summary_table.rowCount()
                self.datasets_summary_table.insertRow(row)
            else:
                row = self.get_row_of_sample_in_table(self.datasets_summary_table,xtal)
            db_dict = collectedXtalsDict[xtal]
            self.update_row_in_table(xtal, row, db_dict, self.datasets_summary_table,
                                     self.datasets_summary_table_columns)

        self.datasets_summary_table.resizeRowsToContents()
        self.datasets_summary_table.resizeColumnsToContents()

        self.status_bar.showMessage('updating Overview table')

        self.status_bar.showMessage('idle')


    def get_selected_row(self,table):
        indexes = table.selectionModel().selectedRows()
        for index in sorted(indexes):
            selected_row = index.row()
        return selected_row

    def show_results_from_all_pipelines(self):
        selected_row=self.get_selected_row(self.datasets_summary_table)
        xtal = self.datasets_summary_table.item(selected_row, 0).text()
        # get details of currently selected autoprocessing result
        selectedResultDict = self.db.get_db_dict_for_sample(xtal)

        dbList=self.db.all_autoprocessing_results_for_xtal_as_dict(xtal)
        self.data_collection_table.setRowCount(0)       # reset table
        for db_dict in dbList:
            row = self.data_collection_table.rowCount()
            self.data_collection_table.insertRow(row)
            self.update_row_in_table(xtal, row, db_dict, self.data_collection_table, self.data_collection_table_columns)
            if selectedResultDict['DataCollectionVisit'] == db_dict['DataCollectionVisit'] \
                and selectedResultDict['DataCollectionRun'] == db_dict['DataCollectionRun'] \
                and selectedResultDict['DataProcessingProgram'] == db_dict['DataProcessingProgram']:
                self.current_row = row
                self.data_collection_table.selectRow(row)
        self.data_collection_table.cellClicked.connect(self.select_different_autoprocessing_result)
        self.data_collection_table_popup()

    def data_collection_table_popup(self):
        self.msgBox = QtGui.QMessageBox()
        msgBoxLayout = self.msgBox.layout()
        qWid = QtGui.QWidget()
        qWid.setFixedWidth(2000)
        qWid.setFixedHeight(500)
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.data_collection_table)
        qWid.setLayout(vbox)
#        msgBoxLayout.addLayout(vbox, 0, 0)
        msgBoxLayout.addWidget(qWid)
        self.msgBox.addButton(QtGui.QPushButton('Cancel'), QtGui.QMessageBox.RejectRole)
        self.msgBox.resize(1000,200)
        self.msgBox.exec_();

    def select_different_autoprocessing_result(self):
        selected_row=self.get_selected_row(self.data_collection_table)
        if selected_row != self.current_row:
            xtal =     self.data_collection_table.item(selected_row, 0).text()
            visit =    self.data_collection_table.item(selected_row, 1).text()
            run =      self.data_collection_table.item(selected_row, 2).text()
            autoproc = self.data_collection_table.item(selected_row, 3).text()
            # get db_dict from collectionTable for visit, run, autoproc
            dbDict = self.db.get_db_dict_for_visit_run_autoproc(xtal,visit,run,autoproc)
            dbDict['DataProcessingAutoAssigned'] = 'False'
            self.update_log.insert('%s: changing selected autoprocessing result to %s %s %s' %(xtal,visit,run,autoproc))
            self.update_log.insert('%s: updating row in Datasets table' %xtal)
            self.db.update_data_source(str(xtal),dbDict)
            self.update_log.insert('%s: getting updated information from DB mainTable' %xtal)
            dbDict = self.db.get_db_dict_for_sample(xtal)
            row = self.get_row_of_sample_in_table(self.datasets_summary_table,xtal)
            self.update_row_in_table(xtal, row, dbDict, self.datasets_summary_table,
                                     self.datasets_summary_table_columns)
        else:
            print 'nothing to change'
        self.msgBox.done(1)

    def update_outcome_datasets_summary_table(self, sample, outcome):
        rows_in_table = self.datasets_summary_table.rowCount()
        for row in range(rows_in_table):
            if self.datasets_summary_table.item(row, 0).text() == sample:
                cell_text = QtGui.QTableWidgetItem()
                cell_text.setText(outcome)
                self.datasets_summary_table.setItem(row, 3, cell_text)


    def populate_and_update_datasource_table(self):
        self.overview_datasource_table.setColumnCount(len(self.overview_datasource_table_columns))

        # first get a list of all the samples that are already in the table and which will be updated
        samples_in_table = []
        current_row = self.overview_datasource_table.rowCount()
        for row in range(current_row):
            sampleID = str(self.overview_datasource_table.item(row, 0).text())  # this must be the case
            samples_in_table.append(sampleID)

        columns_to_show = self.get_columns_to_show(self.overview_datasource_table_columns)
        n_rows = self.get_rows_with_sample_id_not_null_from_datasource()
        sample_id_column = self.get_columns_to_show(['Sample ID'])

        for row in self.data:
            if str(row[sample_id_column[0]]).lower() == 'none' or str(row[sample_id_column[0]]).replace(' ', '') == '':
                # do not show rows where sampleID is null
                continue
            else:
                if not str(row[sample_id_column[0]]) in samples_in_table:
                    # insert row, this is a new sample
                    x = self.overview_datasource_table.rowCount()
                    self.overview_datasource_table.insertRow(x)
                else:
                    # find row of this sample in data_source_table
                    for present_rows in range(self.overview_datasource_table.rowCount()):
                        if str(row[sample_id_column[0]]) == str(
                                self.overview_datasource_table.item(present_rows, 0).text()):
                            x = present_rows
                            break
            for y, item in enumerate(columns_to_show):
                cell_text = QtGui.QTableWidgetItem()
                if row[item] is None:
                    cell_text.setText('')
                else:
                    cell_text.setText(str(row[item]))
                if self.overview_datasource_table_columns[y] == 'Sample ID':  # assumption is that column 0 is always sampleID
                    cell_text.setFlags(QtCore.Qt.ItemIsEnabled)  # and this field cannot be changed
                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                self.overview_datasource_table.setItem(x, y, cell_text)
        self.overview_datasource_table.setHorizontalHeaderLabels(self.overview_datasource_table_columns)

    def kill_other_pandda_options(self):
        for i in range(0, self.pandda_analyse_data_table.rowCount()):
            checkbox1 = self.pandda_analyse_data_table.cellWidget(i,6)
            checkbox2 = self.pandda_analyse_data_table.cellWidget(i,7)
            checkbox3 = self.pandda_analyse_data_table.cellWidget(i,8)
            if checkbox1.isChecked():
                checkbox2.setChecked(False)
                checkbox3.setChecked(False)
            if checkbox1.isChecked() and checkbox2.isChecked() or checkbox3.isChecked():
                checkbox1.setChecked(False)
            if checkbox2.isChecked() or checkbox3.isChecked():
                checkbox1.setChecked(False)

    def populate_pandda_analyse_input_table(self):

        column_name = self.db.translate_xce_column_list_to_sqlite(self.pandda_table_columns)
        print(column_name)
        for xtal in sorted(self.xtal_db_dict):
            new_xtal = False
            db_dict = self.xtal_db_dict[xtal]
            if os.path.isfile(db_dict['DimplePathToPDB']):
                row = self.pandda_analyse_data_table.rowCount()
                if xtal not in self.pandda_analyse_input_table_dict:
                    self.pandda_analyse_data_table.insertRow(row)
                    current_row = row
                    new_xtal = True
                else:
                    for table_row in range(row):
                        if self.pandda_analyse_data_table.item(table_row, 0).text() == xtal:
                            current_row = table_row
                            break
                for column, header in enumerate(column_name):
                    if header[0]=='Exclude':
                        deselect_button = QtGui.QCheckBox()
                        deselect_button.stateChanged.connect(self.kill_other_pandda_options)
                        self.pandda_analyse_data_table.setCellWidget(current_row, column, deselect_button)

                    elif header[0]=='Ignore':
                        deselect_button = QtGui.QCheckBox()
                        deselect_button.stateChanged.connect(self.kill_other_pandda_options)
                        self.pandda_analyse_data_table.setCellWidget(current_row, column, deselect_button)

                    elif header[0] == 'Sample ID':
                        cell_text = QtGui.QTableWidgetItem()
                        cell_text.setText(str(xtal))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.pandda_analyse_data_table.setItem(current_row, column, cell_text)
                    else:
                        cell_text = QtGui.QTableWidgetItem()
                        cell_text.setText(str(db_dict[header[1]]))
                        if header[0] == 'PanDDA\nStatus':
                            if str(db_dict[header[1]]) == 'running':
                                cell_text.setBackground(QtGui.QColor(100, 230, 150))
                            elif str(db_dict[header[1]]) == 'pending':
                                cell_text.setBackground(QtGui.QColor(20, 100, 230))
                            elif str(db_dict[header[1]]) == 'started':
                                cell_text.setBackground(QtGui.QColor(230, 240, 110))
                            elif str(db_dict[header[1]]) == 'finished':
                                cell_text.setBackground(QtGui.QColor(255, 255, 255))
                            elif 'problem' in str(db_dict[header[1]]):
                                cell_text.setBackground(QtGui.QColor(255, 0, 0))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.pandda_analyse_data_table.setItem(current_row, column, cell_text)
            if new_xtal:
                self.pandda_analyse_input_table_dict[xtal] = []

    def populate_and_update_refinement_table(self):

        panddaList = self.db.execute_statement(
            "select CrystalName,PANDDA_site_index,PANDDA_site_name,RefinementOutcome "
            "from panddaTable where CrystalName is not '' and PANDDA_site_ligand_placed is 'True';")
        panddaDict = {}
        for item in panddaList:
            if str(item[0]) not in panddaDict:
                panddaDict[str(item[0])] = []
            panddaDict[str(item[0])].append([str(item[1]), str(item[2]), str(item[3])])

        column_name = self.db.translate_xce_column_list_to_sqlite(self.refinement_table_columns)
        for xtal in sorted(self.xtal_db_dict):
            new_xtal = False
            db_dict = self.xtal_db_dict[xtal]
            try:
                stage = int(str(db_dict['RefinementOutcome']).split()[0])
                refinementStage = db_dict['RefinementOutcome']
            except ValueError:
                stage = 0
            except IndexError:
                stage = 0

            if stage >= 3:
                row = self.refinement_table.rowCount()
                if xtal not in self.refinement_table_dict:
                    self.refinement_table.insertRow(row)
                    current_row = row
                    new_xtal = True
                else:
                    for table_row in range(row):
                        if self.refinement_table.item(table_row, 0).text() == xtal:
                            current_row = table_row
                            break
                for column, header in enumerate(column_name):
                    if header[0] == 'Sample ID':
                        cell_text = QtGui.QTableWidgetItem()
                        cell_text.setText(str(xtal))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.refinement_table.setItem(current_row, column, cell_text)

                    elif header[0] == 'Refinement\nOutcome':
                        if new_xtal:
                            refinement_outcome_combobox = QtGui.QComboBox()
                            self.populate_refinement_outcome_combobox(refinement_outcome_combobox)
                            self.refinement_table.setCellWidget(current_row, column, refinement_outcome_combobox)
                        else:
                            refinement_outcome_combobox = self.refinement_table_dict[xtal]
                        index = refinement_outcome_combobox.findText(refinementStage, QtCore.Qt.MatchFixedString)
                        refinement_outcome_combobox.setCurrentIndex(index)
                        refinement_outcome_combobox.currentIndexChanged.connect(
                            self.refinement_outcome_combobox_changed)

                    elif header[0] == 'PanDDA site details':
                        try:
                            panddaDict[xtal].insert(0, ['Index', 'Name', 'Status'])
                            outerFrame = QtGui.QFrame()
                            outerFrame.setFrameShape(QtGui.QFrame.Box)
                            grid = QtGui.QGridLayout()
                            for y, entry in enumerate(panddaDict[xtal]):
                                for x, info in enumerate(entry):
                                    frame = QtGui.QFrame()
                                    frame.setFrameShape(QtGui.QFrame.Box)
                                    vbox = QtGui.QVBoxLayout()
                                    vbox.addWidget(QtGui.QLabel(str(entry[x])))
                                    frame.setLayout(vbox)
                                    grid.addWidget(frame, y, x)
                            outerFrame.setLayout(grid)
                            self.refinement_table.setCellWidget(current_row, column, outerFrame)
                        except KeyError:
                            cell_text = QtGui.QTableWidgetItem()
                            cell_text.setText('*** N/A ***')
                            cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                            self.refinement_table.setItem(current_row, column, cell_text)
                    else:
                        cell_text = QtGui.QTableWidgetItem()
                        cell_text.setText(str(db_dict[header[1]]))
                        if header[0] == 'Refinement\nStatus':
                            if str(db_dict[header[1]]) == 'running':
                                cell_text.setBackground(QtGui.QColor(100, 230, 150))
                            elif str(db_dict[header[1]]) == 'pending':
                                cell_text.setBackground(QtGui.QColor(20, 100, 230))
                            elif str(db_dict[header[1]]) == 'started':
                                cell_text.setBackground(QtGui.QColor(230, 240, 110))
                            elif str(db_dict[header[1]]) == 'finished':
                                cell_text.setBackground(QtGui.QColor(255, 255, 255))
                            elif 'problem' in str(db_dict[header[1]]):
                                cell_text.setBackground(QtGui.QColor(255, 0, 0))
                        cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                        self.refinement_table.setItem(current_row, column, cell_text)
            if new_xtal:
                self.refinement_table_dict[xtal] = refinement_outcome_combobox

        self.refinement_table.resizeColumnsToContents()
        self.refinement_table.resizeRowsToContents()

    def get_columns_to_show(self, column_list):
        # maybe I coded some garbage before, but I need to find out which column name in the
        # data source corresponds to the actually displayed column name in the table
        # reason being that the unique column ID for DB may not be nice to look at
        columns_to_show = []
        for column in column_list:
            # first find out what the column name in the header is:
            column_name = ''
            for name in self.all_columns_in_data_source:
                if column == name[1]:
                    column_name = name[0]
            for n, all_column in enumerate(self.header):
                if column_name == all_column:
                    columns_to_show.append(n)
                    break
        return columns_to_show

    def get_rows_with_sample_id_not_null_from_datasource(self):
        sample_id_column = self.get_columns_to_show(['Sample ID'])
        n_rows = 0
        for row in self.data:
            if not str(row[sample_id_column[0]]).lower() != 'none' or not str(row[sample_id_column[0]]).replace \
                        (' ', '') == '':
                n_rows += 1
        return n_rows

    def update_data_source(self, sample, db_dict):
        data_source = XChemDB.data_source(os.path.join(self.database_directory, self.data_source_file))

    def quit_xce(self):
        # save pkl file
        if self.data_collection_dict != {}:
            if os.path.isfile(self.datasets_summary_file):
                self.update_log.insert('saving results to PKL file')
                pickle.dump(self.data_collection_dict, open(self.datasets_summary_file, 'wb'))
        self.update_log.insert('quitting XCE... bye,bye!')
        QtGui.qApp.quit()


if __name__ == "__main__":
    app = XChemExplorer(sys.argv[1:])


# "Debugging is twice as hard as writing the code in the first
#  place. Therefore, if you write the code as cleverly as
#  possible, you are, by definition, not smart enough to debug it."
#  -- Brian W. Kernighan
