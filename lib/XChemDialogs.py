import sys
from PyQt4 import QtGui, QtCore

import XChemDB


class select_columns_to_show(QtGui.QDialog):
    def __init__(self,tmp, parent = None):
        super(select_columns_to_show, self).__init__(parent)

        layout = QtGui.QVBoxLayout(self)

        # nice widget for editing the date
        self.datetime = QtGui.QDateTimeEdit(self)
        self.datetime.setCalendarPopup(True)
        self.datetime.setDateTime(QtCore.QDateTime.currentDateTime())
        layout.addWidget(self.datetime)

        print tmp
        print len(tmp)


        #self.data_source=tmp[0]
        #self.active_columns=tmp[1]
        #print self.data_source
        #print self.active_columns
#        db_columns=XChemDB.data_source(self.data_source).return_column_list()
#        self.hallo='hallo'
#        layout = QtGui.QVBoxLayout(self)
#        for item in db_columns:
#            if item[1] in self.active_columns:
#                print 'present',item[1]
#            else:
#                print 'NOT present',item[1]
#            cb = QtGui.QCheckBox(item[1], self)


        # OK and Cancel buttons
        buttons = QtGui.QDialogButtonBox(
            QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
            QtCore.Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    # get current date and time from the dialog
    def dateTime(self):
        return self.datetime.dateTime()

    # static method to create the dialog and return (date, time, accepted)
    @staticmethod
    def show_columns(parent = None):
        dialog = select_columns_to_show(parent)
        result = dialog.exec_()
        date = dialog.dateTime()
        return (date.date(), date.time(), result == QtGui.QDialog.Accepted)





# stuff below works in principle
#class select_columns_to_show(QtGui.QDialog):
#    def __init__(self, parent = None):
#        super(select_columns_to_show, self).__init__(parent)
#
#        layout = QtGui.QVBoxLayout(self)
#
#        # nice widget for editing the date
#        self.datetime = QtGui.QDateTimeEdit(self)
#        self.datetime.setCalendarPopup(True)
#        self.datetime.setDateTime(QtCore.QDateTime.currentDateTime())
#        layout.addWidget(self.datetime)
#
#        # OK and Cancel buttons
#        buttons = QtGui.QDialogButtonBox(
#            QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
#            QtCore.Qt.Horizontal, self)
#        buttons.accepted.connect(self.accept)
#        buttons.rejected.connect(self.reject)
#        layout.addWidget(buttons)
#
#    # get current date and time from the dialog
#    def dateTime(self):
#        return self.datetime.dateTime()
#
#    # static method to create the dialog and return (date, time, accepted)
#    @staticmethod
#    def getDateTime(parent = None):
#        dialog = select_columns_to_show(parent)
#        result = dialog.exec_()
#        date = dialog.dateTime()
#        return (date.date(), date.time(), result == QtGui.QDialog.Accepted)


