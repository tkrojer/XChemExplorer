import sys
from PyQt4 import QtGui, QtCore


class Select(QtGui.QApplication):
    def __init__(self,args):
        QtGui.QApplication.__init__(self,args)
        self.select_columns_from_data_source_dialog()
        self.exec_()




    def select_columns_from_data_source_dialog(self):
        self.window=QtGui.QWidget()
        self.window.setGeometry(0,0, 800,600)
        self.window.setWindowTitle("Select Columns")
        self.center_dialog_window()

        column_checkbutton_hbox=QtGui.QHBoxLayout()
        select_column_button = QtGui.QCheckBox('XXX button name XXX')
        select_column_button.toggle()
#        select_column_button.connect(self.set_run_dimple_flag)
        column_checkbutton_hbox.addWidget(select_sample_for_dimple)

        self.window.addLayout(column_checkbutton_hbox)


        self.window.show()

    def center_dialog_window(self):
        screen = QtGui.QDesktopWidget().screenGeometry()
        size = self.window.geometry()
        self.window.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)


if __name__ == "__main__":
    app=Select(sys.argv)
