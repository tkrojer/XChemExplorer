import os,sys
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

from XChemUtils import parse
from XChemUtils import pdbtools
from XChemUtils import misc
import XChemDB
import csv


def parse_pdb(inital_model_directory,xtal,db_dict):
    if os.path.isfile(os.path.join(inital_model_directory,xtal,'refine.pdb')):
        db_dict['RefinementPDB_latest']=os.path.realpath(os.path.join(inital_model_directory,xtal,'refine.pdb'))
        pdb=parse().PDBheader(os.path.join(inital_model_directory,xtal,'refine.pdb'))
        db_dict['RefinementRcryst'] =               pdb['Rcryst']
        db_dict['RefinementRcrystTraficLight'] =    pdb['RcrystTL']
        db_dict['RefinementRfree']=                 pdb['Rfree']
        db_dict['RefinementRfreeTraficLight'] =     pdb['RfreeTL']
        db_dict['RefinementRmsdBonds']  =           pdb['rmsdBonds']
        db_dict['RefinementRmsdBondsTL'] =          pdb['rmsdBondsTL']
        db_dict['RefinementRmsdAngles'] =           pdb['rmsdAngles']
        db_dict['RefinementRmsdAnglesTL'] =         pdb['rmsdAnglesTL']
        db_dict['RefinementSpaceGroup'] =           pdb['SpaceGroup']
        db_dict['RefinementResolution'] =           pdb['ResolutionHigh']
        db_dict['RefinementResolutionTL'] =         pdb['ResolutionColor']

    return db_dict

def parse_mtz(inital_model_directory,xtal,db_dict):
    if os.path.isfile(os.path.join(inital_model_directory,xtal,'refine.mtz')):
        db_dict['RefinementMTZ_latest']=os.path.realpath(os.path.join(inital_model_directory,xtal,'refine.mtz'))
    return db_dict

def check_refmac_matrix_weight(refinement_directory,db_dict):
    if os.path.isfile(os.path.join(refinement_directory,'refmac.log')):
        for line in open(os.path.join(refinement_directory,'refmac.log')):
            if line.startswith(' Weight matrix') and len(line.split()) == 3:
                db_dict['RefinementMatrixWeight'] = line.split()[2]
    return db_dict

def parse_molprobity_output(inital_model_directory,xtal,db_dict):
    if os.path.isfile(os.path.join(inital_model_directory,xtal,'validation_summary.txt')):
        for line in open(os.path.join(inital_model_directory,xtal,'validation_summary.txt')):
            if 'molprobity score' in line.lower():
                if len(line.split()) >= 4:
                    db_dict['RefinementMolProbityScore'] = line.split()[3]
                    if float(line.split()[3]) < 2:
                        db_dict['RefinementMolProbityScoreTL'] = 'green'
                    if float(line.split()[3]) >= 2 and float(line.split()[3]) < 3:
                        db_dict['RefinementMolProbityScoreTL'] = 'orange'
                    if float(line.split()[3]) >= 3:
                        db_dict['RefinementMolProbityScoreTL'] = 'red'

            if 'ramachandran outliers' in line.lower():
                if len(line.split()) >= 4:
                    db_dict['RefinementRamachandranOutliers'] = line.split()[3]

                    if float(line.split()[3]) < 0.3:
                        db_dict['RefinementRamachandranOutliersTL'] = 'green'
                    if float(line.split()[3]) >= 0.3 and float(line.split()[3]) < 1:
                        db_dict['RefinementRamachandranOutliersTL'] = 'orange'
                    if float(line.split()[3]) >= 1:
                        db_dict['RefinementRamachandranOutliersTL'] = 'red'

            if 'favored' in line.lower():
                if len(line.split()) >= 3:
                    db_dict['RefinementRamachandranFavored'] = line.split()[2]
                    if float(line.split()[2]) < 90:
                        db_dict['RefinementRamachandranFavoredTL'] = 'red'
                    if float(line.split()[2]) >= 90 and float(line.split()[2]) < 98:
                        db_dict['RefinementRamachandranFavoredTL'] = 'orange'
                    if float(line.split()[2]) >= 98:
                        db_dict['RefinementRamachandranFavoredTL'] = 'green'

    return db_dict


def parse_ligand_validation(refinement_directory,xtal):
    if os.path.isfile(os.path.join(refinement_directory, 'residue_scores.csv')):
        with open(os.path.join(refinement_directory, 'residue_scores.csv'), 'rb') as csv_import:
            csv_dict = csv.DictReader(csv_import)
            for i, line in enumerate(csv_dict):
                db_pandda_dict = {}
                residue = line['']
                if len(residue.split('-')) == 3:
                    residue_name = residue.split('-')[0]
                    residue_chain = residue.split('-')[1]
                    residue_number = residue.split('-')[2]
                    residue_xyz = pdbtools(os.path.join(inital_model_directory, xtal, 'refine.pdb')).get_center_of_gravity_of_residue_ish(residue_chain, residue_number)
                    event = db.execute_statement("select PANDDA_site_x,PANDDA_site_y,PANDDA_site_z,PANDDA_site_index from panddaTable where CrystalName='%s'" % xtal)
                    print event
                    for coord in event:
                        db_pandda_dict = {}
                        event_x = float(str(coord[0]))
                        event_y = float(str(coord[1]))
                        event_z = float(str(coord[2]))
                        site_index = str(coord[3])
                        print event_x,event_y,event_z,site_index
                        distance = misc().calculate_distance_between_coordinates(residue_xyz[0], residue_xyz[1],residue_xyz[2],
                                                                                 event_x, event_y,event_z)
                        print 'distance',distance
                        # if coordinate of ligand and event are closer than 5A, then we assume they belong together
                        if distance < 5:
                            print 'hallo'
                            db_pandda_dict['PANDDA_site_ligand_id'] = residue
                            db_pandda_dict['PANDDA_site_occupancy'] = line['Occupancy']
                            db_pandda_dict['PANDDA_site_B_average'] = line['Average B-factor (Residue)']
                            db_pandda_dict['PANDDA_site_B_ratio_residue_surroundings'] = line['Surroundings B-factor Ratio']
                            db_pandda_dict['PANDDA_site_rmsd'] = line['Model RMSD']
                            db_pandda_dict['PANDDA_site_RSCC'] = line['RSCC']
                            db_pandda_dict['PANDDA_site_RSR'] = line['RSR']
                            db_pandda_dict['PANDDA_site_RSZD'] = line['RSZD']
                            if os.path.isfile(os.path.join(refinement_directory, 'residue_plots', residue + '.png')):
                                db_pandda_dict['PANDDA_site_spider_plot'] = os.path.join(refinement_directory,
                                                                                         'residue_plots',
                                                                                         residue + '.png')
                        print db_pandda_dict
                        if db_pandda_dict != {}:
                            print '==> XCE: updating pandda Table of data source'
                            db.update_panddaTable(xtal, site_index, db_pandda_dict)

def update_data_source(db_dict):
    if db_dict != {}:
        print '==> xce: updating mainTable of data source'
        db.update_data_source(xtal,db_dict)
        # update refinement outcome if necessary
        sqlite = (
            "update mainTable set RefinementOutcome = '3 - In Refinement' where CrystalName is '%s' " %xtal+
            "and (RefinementOutcome is null or RefinementOutcome is '1 - Analysis Pending' or RefinementOutcome is '2 - PANDDA model')"
                )
        db.execute_statement(sqlite)
        # now do the same for each site in the pandda table
        sqlite = (
            "update panddaTable set RefinementOutcome = '3 - In Refinement' where CrystalName is '%s' " %xtal+
            "and (RefinementOutcome is null or RefinementOutcome is '1 - Analysis Pending' or RefinementOutcome is '2 - PANDDA model')"
                )
        db.execute_statement(sqlite)

if __name__=='__main__':
    db_file=sys.argv[1]
    xtal=sys.argv[2]
    inital_model_directory=sys.argv[3]
    refinement_directory=sys.argv[4]

    db=XChemDB.data_source(db_file)
    db_dict={}

    db_dict=parse_pdb(inital_model_directory,xtal,db_dict)
    db_dict=parse_mtz(inital_model_directory,xtal,db_dict)
    db_dict=check_refmac_matrix_weight(refinement_directory,db_dict)
    db_dict=parse_molprobity_output(inital_model_directory,xtal,db_dict)

    parse_ligand_validation(refinement_directory,xtal)

    update_data_source(db_dict)
