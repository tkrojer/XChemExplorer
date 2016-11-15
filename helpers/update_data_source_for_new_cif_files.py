# last edited: 15/11/2016, 15:00

import os,sys
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

from XChemUtils import parse
import XChemDB

if __name__=='__main__':
    db_file=sys.argv[1]
    xtal=sys.argv[2]
    inital_model_directory=sys.argv[3]
    compoundID=sys.argv[4]

    db=XChemDB.data_source(db_file)
    db_dict={}
    if os.path.isfile(os.path.join(inital_model_directory,xtal,'compound',compoundID+'.cif')):
        db_dict['RefinementCIF']=os.path.join(inital_model_directory,xtal,'compound',compoundID+'.cif')
        db_dict['RefinementCIFStatus']='restraints\ngenerated'
    else:
        db_dict['RefinementCIF']=''
        db_dict['RefinementCIFStatus']='restraints\nfailed'
    db.update_data_source(xtal,db_dict)