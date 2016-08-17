import getpass

def dataset_task_tip():
    tip =   (   'describes what you can do'   )
    return tip

def dataset_task_run_button_tip():
    tip =   (   'Run dataset'   )
    return tip

def dataset_task_status_button_tip():
    tip =   (   'Status dataset'   )
    return tip

def map_cif_file_task_tip():
    tip =   (   'describes what you can do'   )
    return tip

def map_cif_file_task_run_button_tip():
    tip =   (   'Run map_cif_file'   )
    return tip

def map_cif_file_task_status_button_tip():
    tip =   (   'Status map_cif_file'   )
    return tip

def panddas_file_task_tip():
    tip =   (   'describes what you can do'   )
    return tip

def panddas_file_task_run_button_tip():
    tip =   (   'Run panddas_file'   )
    return tip

def panddas_file_task_status_button_tip():
    tip =   (   'Status panddas_file'   )
    return tip

def refine_file_task_tip():
    tip =   (   'describes what you can do'   )
    return tip

def refine_file_task_run_button_tip():
    tip =   (   'Run refine_file'   )
    return tip

def refine_file_task_status_button_tip():
    tip =   (   'Status refine_file'   )
    return tip

def validation_file_task_tip():
    tip =   (   'describes what you can do'   )
    return tip

def validation_file_task_run_button_tip():
    tip =   (   'Run validation_file'   )
    return tip

def validation_file_task_status_button_tip():
    tip =   (   'Status validation_file'   )
    return tip

def update_from_datasource_button_tip():
    tip =   (   'Status validation_file'   )
    return tip




def load_samples_from_datasource():
    tip =   (   'loads all information\n'
                'which is currently in the\n'
                'data source'   )
    return tip

def save_samples_to_datasource():
    tip =   (   'saves all changes you made\n'
                'to any of the cells to the\n'
                'data source'   )
    return tip

def create_pdb_cif_png_files():
    tip =   (   'uses ACEDRG to create pdb,cif and png files for all compounds\n'
                'that have a smiles string assigned. The compounds ID\n'
                'serves as name if available.\n'
                'The files are written to initial_model/<sample ID>/compound\n'
                'data source'   )
    return tip

def check_status_create_pdb_cif_png_files():
    tip =   (   'the jobs are sent to a computer cluster if available\n'
                'gives you an overview about how it progresses')
    return tip

def run_pandda_inspect_at_home(pandda_directory):
    instruction =   (   ' 1. Make sure that the pandda package is installed at your home institution.\n'
                        '    For more details see: http://pandda.bitbucket.org\n'
                        ' 2. At home: go to the folder where you want the files from pandda.analyse to be;\n'
                        '    e.g. cd /home/tkrojer/fragment_screening\n'
                        ' 3. At home: run the following command:\n'
                        '    rsync -av %s@nx.diamond.ac.uk:%s .\n' %(getpass.getuser(),pandda_directory)+
                        ' 4. At home: go into pandda direcotory\n'
                        '    cd %s\n' %pandda_directory[pandda_directory.rfind('/')+1:]+
                        ' 5. At home: run pandda.inspect\n'
                        '    pandda.inspect\n'
                        ' 6. Once you finished inspecting your models, close pandda.inspect and copy the data back to DLS\n'
                        '    rsync -av * %s@nx.diamond.ac.uk:%s\n' %(getpass.getuser(),pandda_directory)+
                        ' 7. At DLS: continue using XChemExplorer; go to the PANDDA tab an run:\n'
                        '    "Export PANDDA models"\n'
                        '    This will trigger the follwing steps:\n'
                        '    - transform the pandda models back into the crystallographic unit cell\n'
                        '    - copy the transformed models into the respective folder of the project directory\n'
                        '    - launch an initial round of refinement with REFMAC\n'
                        ' 8. At DLS: go to the REFINEMENT tab and run "Open COOT" to review and further refine your models\n' )
    print instruction