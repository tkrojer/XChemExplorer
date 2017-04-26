from import_modules import *

def set_run_dimple_flag(self, state):
    if state == QtCore.Qt.Checked:
        for key in self.initial_model_dimple_dict:
            self.initial_model_dimple_dict[key][0].setChecked(True)
    else:
        for key in self.initial_model_dimple_dict:
            self.initial_model_dimple_dict[key][0].setChecked(False)

def set_new_reference_if_applicable(self):
    print 'hallo'
    reference_root=str(self.reference_file_selection_combobox.currentText())
    pg_ref=''
    ucVol_ref=0.0
    for reference in self.reference_file_list:
        print reference[0],reference_root
        if reference[0]==reference_root:
            pg_ref=reference[5]
            ucVol_ref=reference[4]
            break
    if ucVol_ref==0.0:
        self.update_log.insert('cannot set reference file since unit cell volume of reference pdb is 0!')
        return

    for xtal in self.initial_model_dimple_dict:
        reference_file_selection_combobox=self.initial_model_dimple_dict[xtal][1]
        db_dict=self.xtal_db_dict[xtal]
        pg_xtal=db_dict['DataProcessingPointGroup']
        ucVol_xtal=db_dict['DataProcessingUnitCellVolume']

        try:
            difference=math.fabs(1-(float(ucVol_xtal)/float(ucVol_ref)))*100
        except ValueError:
            self.update_log.insert(xtal+' -> cannot calculate unit cell volume difference')
            continue

        if pg_xtal==pg_ref and difference < self.allowed_unitcell_difference_percent:
            print xtal,pg_xtal,ucVol_xtal
            index = reference_file_selection_combobox.findText(reference_root, QtCore.Qt.MatchFixedString)
            reference_file_selection_combobox.setCurrentIndex(index)
            self.update_log.insert(xtal+' -> setting '+reference_root+' as input PDB file for DIMPLE')

def on_context_menu_initial_model(self, point):
    # show context menu
    self.popMenu_for_initial_model_table.exec_(self.sender().mapToGlobal(point))
