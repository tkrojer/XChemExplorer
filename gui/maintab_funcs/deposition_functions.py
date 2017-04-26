from import_modules import *

def prepare_files_for_zenodo_upload(self):
    self.update_log.insert('preparing files for ZENODO upload...')
    os.system('ccp4-python ' + os.getenv(
        'XChemExplorer_DIR') + '/helpers/prepare_for_zenodo_upload.py ' + self.html_export_directory)