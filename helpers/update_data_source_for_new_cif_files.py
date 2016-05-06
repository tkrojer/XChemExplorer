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
    if os.path.isfile(os.path.join(inital_model_directory,xtal,'compound',compoundID+'.cif')):
        db_dict={}
        db_dict['RefinementCIF']=os.path.join(inital_model_directory,xtal,'compound',compoundID+'.cif')
        db.update_data_source(xtal,db_dict)