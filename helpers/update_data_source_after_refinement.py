import os,sys
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

from XChemUtils import parse
from XChemUtils import pdbtools
from XChemUtils import misc
import XChemDB
import csv

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

        if os.path.isfile(os.path.join(refinement_directory,'residue_scores.csv')):
            with open(os.path.join(refinement_directory,'residue_scores.csv'),'rb') as csv_import:
                csv_dict = csv.DictReader(csv_import)
                for i,line in enumerate(csv_dict):
                    residue=line['']
                    if len(residue.split('-'))==3:
                        residue_name=residue.split('-')[0]
                        residue_chain=residue.split('-')[1]
                        residue_number=residue.split('-')[2]
                        residue_xyz=pdbtools(os.path.join(inital_model_directory,xtal,'refine.pdb')).get_center_of_gravity_of_residue_ish(residue_chain,residue_number)

                        event=db.execute_statement("select PANDDA_site_x,PANDDA_site_y,PANDDA_site_z,PANDDA_site_index from panddaTable where CrystalName='%s'" %xtal)
                        for coord in event:
                            db_pandda_dict={}
                            event_x=float(str(coord[0]))
                            event_y=float(str(coord[1]))
                            event_z=float(str(coord[2]))
                            site_index=str(coord[3])
                            distance=misc().calculate_distance_between_coordinates(residue_xyz[0],residue_xyz[1],residue_xyz[2],event_x,event_y,event_z)
                            # if coordinate of ligand and event are closer than 5A, then we assume they belong together
                            if distance < 5:
                                db_pandda_dict['PANDDA_site_ligand_id']=residue
                                if os.path.isfile(os.path.join(refinement_directory,'residue_plots',residue+'.png')):
                                    db_pandda_dict['PANDDA_site_spider_plot']=os.path.join(refinement_directory,'residue_plots',residue+'.png')
                            if db_pandda_dict != {}:
                                db.update_panddaTable(xtal,site_index,db_pandda_dict)

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

#    site=db.execute_statement()