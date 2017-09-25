import sys, os, multiprocessing
from PyQt4 import QtGui, QtCore, QtWebKit

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'gui_scripts'))

import layout


class PanddaTab():
    def __init__(self):
        self.layout_funcs = layout.LayoutFuncs()

    def setup(self, xce_object):
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
