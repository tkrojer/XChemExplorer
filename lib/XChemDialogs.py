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


        self.window.setLayout(vbox_main)

        self.window.show()

    def center_dialog_window(self):
        screen = QtGui.QDesktopWidget().screenGeometry()
        size = self.window.geometry()
        self.window.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)


if __name__ == "__main__":
    app=Select(sys.argv)
