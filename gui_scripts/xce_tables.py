import sys, os
from PyQt4 import QtGui, QtCore, QtWebKit

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'gui_scripts'))

import layout


class tables():
    def __init__(self,xce_object,table):
        self.layout_funcs = layout.LayoutFuncs()

        if table == xce_object.deposition_table_apo

    def update(self, xce_object,sample,row,db_dict,table,columns_to_show):
        ################################################################################################################
        #                                                                                                              #
        #                                                REFINEMENT TAB                                                #
        #                                                                                                              #
        ################################################################################################################
#        xce_object.summary_vbox_for_table = QtGui.QVBoxLayout()
#    def update_row_in_table(self,sample,row,db_dict,table,columns_to_show):
        xtal = str(sample)
        column_name = xce_object.db.translate_xce_column_list_to_sqlite(columns_to_show)

        for column, header in enumerate(column_name):

            if header[0] == 'Sample ID':
                cell_text = QtGui.QTableWidgetItem()
                cell_text.setText(str(xtal))
                cell_text.setTextAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter)
                table.setItem(row, column, cell_text)

            elif header[0] == 'DataCollection\nOutcome':
                if xtal not in xce_object.dataset_outcome_combobox_dict:
                    dataset_outcome_combobox = QtGui.QComboBox()
                    for outcomeItem in xce_object.dataset_outcome:
                        dataset_outcome_combobox.addItem(outcomeItem)
                    dataset_outcome_combobox.activated[str].connect(xce_object.dataset_outcome_combobox_change_outcome)
                    xce_object.dataset_outcome_combobox_dict[xtal] = dataset_outcome_combobox
                    table.setCellWidget(row, column, dataset_outcome_combobox)
                index = xce_object.dataset_outcome_combobox_dict[xtal].findText(str(db_dict['DataCollectionOutcome']), QtCore.Qt.MatchFixedString)
                xce_object.dataset_outcome_combobox_dict[xtal].setCurrentIndex(index)

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
                if table == xce_object.deposition_table_apo:
                    if xtal not in xce_object.deposition_table_apo_dict:
                        xce_object.deposition_table_apo_dict[xtal] = checkbox
                if table == xce_object.deposition_table_bound:
                    if xtal not in xce_object.deposition_table_bound_dict:
                        xce_object.deposition_table_bound_dict[xtal] = checkbox
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
                #                            xce_object.datasets_summary_table.setItem(current_row, column, cell_text)
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
