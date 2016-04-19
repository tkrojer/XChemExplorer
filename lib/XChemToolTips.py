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
