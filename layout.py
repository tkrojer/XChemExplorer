import sys, os, subprocess
from PyQt4 import QtGui, QtCore

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'lib'))
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'web'))

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

    # function for opening help and tutorial files
    def openFile(self, file):
        if sys.platform == 'linux2':
            subprocess.call(["xdg-open", file])
        else:
            os.startfile(file)

    def add_to_box(self, frame, widgets_list):
        for widget in widgets_list:
            frame.addWidget(widget)
