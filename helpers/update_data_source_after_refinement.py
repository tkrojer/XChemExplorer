import os,sys
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

from XChemUtils import parse
import XChemDB

if __name__=='__main__':
    db_file=sys.argv[1]
    xtal=sys.argv[2]
    inital_model_directory=sys.argv[3]
    refinement_directory=sys.argv[4]

    db=XChemDB.data_source(db_file)
    db_dict={}
    if os.path.isfile(os.path.join(inital_model_directory,xtal,'refine.pdb')):
        db_dict['RefinementPDB_latest']=os.path.realpath(os.path.join(inital_model_directory,xtal,'refine.pdb'))
        pdb=parse().PDBheader(os.path.join(inital_model_directory,xtal,'refine.pdb'))
        db_dict['RefinementRcryst']=pdb['Rcryst']
        db_dict['RefinementRfree']=pdb['Rfree']
        db_dict['RefinementRmsdBonds']=pdb['rmsdBonds']
        db_dict['RefinementRmsdAngles']=pdb['rmsdAngles']
        db_dict['RefinementSpaceGroup']=pdb['SpaceGroup']
    if os.path.isfile(os.path.join(inital_model_directory,xtal,'refine.mtz')):
        db_dict['RefinementMTZ_latest']=os.path.realpath(os.path.join(inital_model_directory,xtal,'refine.mtz'))

    if db_dict != {}:
        # finally, update data source
        print '==> xce: updating data source after DIMPLE run'
        db.update_data_source(xtal,db_dict)
        # update refinement outcome if necessary
        sqlite = (
            "update mainTable set RefinementOutcome = '3 - In Refinement' where CrystalName is '%s' " %xtal+
            "and RefinementOutcome is null or RefinementOutcome is '1 - Analysis Pending' or RefinementOutcome is '2 - PANDDA model'"
                )
        db.execute_statement(sqlite)

    site=db.execute_statement()