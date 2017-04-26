from import_modules import *

def select_pandda_input_template(self):
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
    combo_text=str(self.pandda_reference_file_selection_combobox.currentText())
    for file in self.reference_file_list:
        if file[0] == combo_text:
            self.pandda_reference_file_spg_label.setText(file[1])
            break