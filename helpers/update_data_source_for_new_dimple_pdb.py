import os,sys
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

from XChemUtils import parse
import XChemDB

if __name__=='__main__':
    db_file=sys.argv[1]
    xtal=sys.argv[2]
    inital_model_directory=sys.argv[3]

    db=XChemDB.data_source(db_file)
    if os.path.isfile(os.path.join(inital_model_directory,xtal,'dimple.pdb')):
        db_dict={}
        db_dict['DimplePathToPDB']=os.path.join(inital_model_directory,xtal,'dimple.pdb')
        if os.path.isfile(os.path.join(inital_model_directory,xtal,'dimple.mtz')):
            db_dict['DimplePathToMTZ']=os.path.join(inital_model_directory,xtal,'dimple.mtz')
        print os.path.join(inital_model_directory,xtal,'dimple.pdb')
        pdb=parse().PDBheader(os.path.join(inital_model_directory,xtal,'dimple.pdb'))
        db_dict['DimpleRcryst']=pdb['Rcryst']
        db_dict['DimpleRfree']=pdb['Rfree']
        if os.path.isfile(os.path.join(inital_model_directory,xtal,'dimple','dimple_rerun_on_selected_file','dimple','prepared2.mtz')):
            os.chdir(os.path.join(inital_model_directory,xtal))
            if os.path.isfile(xtal+'.free.mtz')):
                os.system('/bin/rm '+xtal+'.free.mtz')
            os.symlink(os.path.join('dimple','dimple_rerun_on_selected_file','dimple','prepared2.mtz'),xtal+'.free.mtz')
            db_dict['RefinementMTZfree']=xtal+'.free.mtz'
        db.update_data_source(xtal,db_dict)