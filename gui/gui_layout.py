from import_modules import *
import menu_bar as menbar
import task_containers as tcont
import main_tabs as tabs

#    ,----------------------------------------------,  <------
#   |            menu bar - menu_bar.py             |        |
#   -------------------------------------------------        |
#   |                  main tabs                    |<----   |
#   |     __________________________________        |    |   |
#   |    |                                  |       |    |   |
#   |    |                                  |       |    |   |
#   |    |              sub tabs            | main_tabs.py   |
#   |    |                                  |       |    |   |
#   |    |                                  |       |        self.window (layout: vbox_main)
#   |    |__________________________________|       |    |   |
#   |                                               |<----   |
#   -------------------------------------------------        |
#   |         bottom bar - task_containers.py       |        |
#   '-----------------------------------------------'  <------

def setup_layout(self):
    # Set up a widget to contain everything
    self.window=QtGui.QWidget()
    self.window.setWindowTitle("XChemExplorer")

    # define screen size for resizing
    self.screen = QtGui.QDesktopWidget().screenGeometry()

    # create a layout object to hold pieces
    self.vbox_main = QtGui.QVBoxLayout()

    # build menubar and task container (bottom bar)
    menbar.menu_bar(self)
    tcont.setup_for_layout(self)

    # check datasource permissions
    if self.data_source_file != '':
        write_enabled = self.check_write_permissions_of_data_source()
        if not write_enabled:
            self.data_source_set = False

    # set up main and sub tabs
    tabs.tabs_setup(self)

    # add all elements to main layout
    self.vbox_main.addWidget(self.menu_bar)
    self.vbox_main.addWidget(self.main_tab_widget)
    self.vbox_main.addLayout(self.hboxTaskFrames)
    self.vbox_main.addLayout(self.hbox_status)

    # put the main layout onto the main window
    self.window.setLayout(self.vbox_main)

    # display xce
    self.status_bar.showMessage('Ready')
    self.window.show()