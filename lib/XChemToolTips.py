# last edited: 01/03/2017, 17:00

import os
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
#    instruction =   (   '\n\n'
#                        'Please follow the instructions below to run pandda.inspect at your home institution:\n\n'
#                        ' 1. Make sure that the pandda package is installed at your home institution.\n'
#                        '    For more details see: http://pandda.bitbucket.org\n\n'
#                        ' 2. At home: go to the folder where you want the files from pandda.analyse to be;\n'
#                        '    e.g. cd /home/tkrojer/fragment_screening\n\n'
#                        ' 3. At home: run the following command:\n'
#                        '    rsync -av %s@nx.diamond.ac.uk:%s .\n\n' %(getpass.getuser(),pandda_directory)+
#                        ' 4. At home: go into the pandda directory\n'
#                        '    cd %s\n\n' %pandda_directory[pandda_directory.rfind('/')+1:]+
#                        ' 5. At home: run pandda.inspect\n'
#                        '    pandda.inspect\n\n'
#                        ' 6. Once you finished inspecting your models, close pandda.inspect and copy the data back to DLS\n'
#                        '    with the following command:\n'
#                        '    rsync -av * %s@nx.diamond.ac.uk:%s\n\n' %(getpass.getuser(),pandda_directory)+
#                        ' 7. At DLS: continue using XChemExplorer; go to the PANDDA tab and run:\n'
#                        '    "Export PANDDA models"\n'
#                        '    This will trigger the following steps:\n'
#                        '    - transform the pandda models back into the crystallographic unit cell\n'
#                        '    - copy the transformed models into the respective folder of the project directory\n'
#                        '    - launch an initial round of refinement with REFMAC\n'
#                        '      Note: All refinement jobs run on the cluster at DLS.\n'
#                        '            You can check the status of the jobs by typing:\n'
#                        '            module load global/cluster\n'
#                        '            qstat\n\n'
#                        ' 8. At DLS: go to the REFINEMENT tab and run "Open COOT" to review and further refine your models if necessary\n' )

    instruction =   (   '\n\n'
                        'Be sure to have pandda installed at home, and go to a clean subdirectory.\n'
                        'From that directory, do the steps below.\n'
                        'This moves the relevant files to your site so you can do the model building locally, and then moves the files back to Diamond.\n'
                        '1.	run:  rsync -av %s@nx.diamond.ac.uk:%s .\n' %(getpass.getuser(),pandda_directory)+
                        '2.	run: "pandda.inspect", and build all relevant models, etc.\n'
                        '3.	run:  rsync -av * %s@nx.diamond.ac.uk:%s\n' %(getpass.getuser(),pandda_directory)+
                        'Now proceed within XChemExplorer as before.\n' )

    print instruction

def deposition_interface_note():

    note = (    '\n\n'
                'Note: you can use this mask to save identical information for ALL structures to be deposited.\n'
                'However, this may not be suitable in cases where the information is different for certain samples.\n'
                'In such cases, please use for example SQLiteBrowser to edit the relevant fields in the depositTable.'   )

    return note

def html_summary_introduction():

    text = (    'This section describes how to generate an HTML summary page of the current experiment and \n'
                'how to share it with the scientific community via the ZENODO data repository.\n'
                'For example:')

    return text

def html_export_directory_background():

    text = (    'The HTML export directory will contain all files necessary to generate the summary file.\n'
                'It is also makes it easy to share the results with collaborators.')

def html_export_step():

    text = (    'This step will create an index.html file and it will copy pdb, mtz and cif files\n'
                'as well as PanDDA event maps into the html export directory.'  )

    return text


def icb_file_background():

    text = (    'This step will use ICM-pro to create interactive structure thumbnails in the html file (http://www.molsoft.com/activeicmjs.html)')

    return text



def prepare_ICB_files():

    instruction = (     '- Open file browser and navigate to HTML export directory.\n'
                        '- Open ICM (at DLS: select SERVER and enter diamvicmpro.diamond.ac.uk in case ICM asks).\n'
                        '- Drag-and-drop dsEvent_sqlite.icm file into ICM main window.\n'
                        '- In ICM workspace panel: right-click on dsEvent_sqlite.icm and select "RUN".')
    return instruction

def zenodo_upload_start(html_export_directory):

    instruction = (     'This will copy all files required for ZENODO upload to\n'
                        '%s'   %os.path.join(html_export_directory,'zenodo')    )
    return instruction

def zenodo_upload_part_one(html_export_directory):

    instruction = (     '- Register with ZENODO (www.zenodo.org)\n'
                        '- Select Upload -> New Upload\n'
                        '- upload all files from %s/zenodo, BUT NOT any of the .html files!!!' %html_export_directory  )
    return instruction

def zenodo_upload_part_two():

    instruction = (     'You might have noticed that once you started uploading files, the address bar in your browser will change to something like this:')

    return instruction

def zenodo_upload_part_three():

    instruction = (     'Please enter the upload ID, i.e. the number at the end of the line, into the field below')

    return instruction

def zenodo_upload_part_four(html_export_directory):

    instruction = (     'Upload ALL html files in %s to ZENODO and publish the page.' %os.path.join(html_export_directory,'zenodo'))

    return instruction

def deposition_pandda_site_not_ready(xtal):
    msg = (     'Please make sure that all placed ligands in '+xtal+' are ready for deposition.\n'
                'Go to the Refinement tab -> Open COOT and either delete the ligand or set the site '
                ' to "5 - Deposition ready"'    )
    return msg