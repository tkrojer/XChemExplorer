import sqlite3
import os,sys
import csv

# committed on: 05/12/2015

sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
from XChemUtils import parse

class data_source:

    def __init__(self,data_source_file):

        self.data_source_file=data_source_file

        #   [column_name in DB, column_name shown in XCE, SQLite column type (Integer,Text,PKEY...)]
        self.column_list=[
            # SQLite column name                    XCE column name                             SQLite type
            # from Lab36
            ['ID',                                   'ID',                                      'INTEGER PRIMARY KEY'],
            ['LabVisit',                             'LabVisit',                                'TEXT'],
            ['LibraryPlate',                         'LibraryPlate',                            'TEXT'],
            ['SourceWell',                           'SourceWell',                              'TEXT'],
            ['LibraryName',                          'LibraryName',                             'TEXT'],
            ['CompoundSMILES',                       'Smiles',                                  'TEXT'],
            ['CompoundCode',                         'Compound ID',                             'TEXT'],
            ['CrystalPlate',                         'CrystalPlate',                            'TEXT'],
            ['CrystalWell',                          'CrystalWell',                             'TEXT'],
            ['EchoX',                                'EchoX',                                   'TEXT'],
            ['EchoY',                                'EchoY',                                   'TEXT'],
            ['DropVolume',                           'DropVolume',                              'TEXT'],
            ['ProteinName',                          'ProteinName',                             'TEXT'],
            ['BatchNumber',                          'BatchNumber',                             'TEXT'],
            ['CompoundStockConcentration',           'CompoundStockConcentration',              'TEXT'],
            ['CompoundConcentration',                'CompoundConcentration',                   'TEXT'],
            ['SolventFraction',                      'SolventFraction',                         'TEXT'],
            ['SoakTransferVol',                      'SoakTransferVol',                         'TEXT'],
            ['SoakStatus',                           'SoakStatus',                              'TEXT'],
            ['SoakTimestamp',                        'SoakTimestamp',                           'TEXT'],
            ['CryoStockFraction',                    'CryoStockFraction',                       'TEXT'],
            ['CryoFraction',                         'CryoFraction',                            'TEXT'],
            ['CryoWell',                             'CryoWell',                                'TEXT'],
            ['CryoTransferVolume',                   'CryoTransferVolume',                      'TEXT'],
            ['CryoStatus',                           'CryoStatus',                              'TEXT'],
            ['CryoTimestamp',                        'CryoTimestamp',                           'TEXT'],
            ['SoakingTime',                          'SoakingTime',                             'TEXT'],
            ['HarvestStatus',                        'HarvestStatus',                           'TEXT'],
            ['CrystalName',                          'Sample ID',                               'TEXT'],
            ['Puck',                                 'Puck',                                    'TEXT'],
            ['PuckPosition',                         'PuckPosition',                            'TEXT'],
            ['PinBarcode',                           'PinBarcode',                              'TEXT'],
            ['MountingResult',                       'MountingResult',                          'TEXT'],
            ['MountingArrivalTime',                  'MountingArrivalTime',                     'TEXT'],
            ['MountedTimestamp',                     'MountedTimestamp',                        'TEXT'],
            ['MountingTime',                         'MountingTime',                            'TEXT'],
            ['ispybStatus',                          'ispybStatus',                             'TEXT'],
            ['DataCollectionVisit',                  'Visit',                     'TEXT'],

            # from XChemExplorer

            ['CompoundName',                         'Compound Name',                           'TEXT'],
            ['CrystalTag',                           'Tag',                                     'TEXT'],
#            ['CrystalForm',                          'Crystal\nForm',                           'TEXT'],
            ['CrystalFormName',                      'Crystal Form\nName',                         'TEXT'],
            ['CrystalFormSpaceGroup',                'Space\nGroup',                   'TEXT'],
            ['CrystalFormPointGroup',                'Point\nGroup',                   'TEXT'],
            ['CrystalFormA',                         'a',                            'TEXT'],
            ['CrystalFormB',                         'b',                            'TEXT'],
            ['CrystalFormC',                         'c',                            'TEXT'],
            ['CrystalFormAlpha',                     'alpha',                        'TEXT'],
            ['CrystalFormBeta',                      'beta',                         'TEXT'],
            ['CrystalFormGamma',                     'gamma',                        'TEXT'],
            ['CrystalFormVolume',                     'Crystal Form\nVolume',                        'TEXT'],


            ['DataCollectionBeamline',               'Beamline',                  'TEXT'],
            ['DataCollectionDate',                   'Data Collection\nDate',                      'TEXT'],
            ['DataCollectionOutcome',                'DataCollection\nOutcome',                   'TEXT'],
            ['DataCollectionRun',                    'Run',                       'TEXT'],
            ['DataCollectionComment',                'DataCollection\nComment',                   'TEXT'],
            ['DataCollectionWavelength',                'Wavelength',                   'TEXT'],

            ['DataProcessingProgram',                'Program',                   'TEXT'],
            ['DataProcessingSpaceGroup',             'DataProcessing\nSpaceGroup',                'TEXT'],
            ['DataProcessingUnitCell',               'DataProcessing\nUnitCell',                  'TEXT'],

            ['DataProcessingA',               'DataProcessing\nA',                  'TEXT'],
            ['DataProcessingB',               'DataProcessing\nB',                  'TEXT'],
            ['DataProcessingC',               'DataProcessing\nC',                  'TEXT'],
            ['DataProcessingAlpha',               'DataProcessing\nAlpha',                  'TEXT'],
            ['DataProcessingBeta',               'DataProcessing\nBeta',                  'TEXT'],
            ['DataProcessingGamma',               'DataProcessing\nGamma',                  'TEXT'],



            ['DataProcessingResolutionOverall',             'Resolution\nOverall',                          'TEXT'],
            ['DataProcessingResolutionLow',                 'Resolution\nLow',                              'TEXT'],
            ['DataProcessingResolutionLowInnerShell',       'Resolution\nLow (Inner Shell)',                'TEXT'],
            ['DataProcessingResolutionHigh',                'Resolution\nHigh',                             'TEXT'],
            ['DataProcessingResolutionHigh15sigma',         'Resolution\n[Mn<I/sig(I)> = 1.5]',             'TEXT'],
            ['DataProcessingResolutionHighOuterShell',      'Resolution\nHigh (Outer Shell)',               'TEXT'],
            ['DataProcessingRmergeOverall',                 'Rmerge\nOverall',                              'TEXT'],
            ['DataProcessingRmergeLow',                     'Rmerge\nLow',                                  'TEXT'],
            ['DataProcessingRmergeHigh',                    'Rmerge\nHigh',                                 'TEXT'],
            ['DataProcessingIsigOverall',                   'Mn<I/sig(I)>\nOverall',                        'TEXT'],
            ['DataProcessingIsigLow',                       'Mn<I/sig(I)>\nLow',                            'TEXT'],
            ['DataProcessingIsigHigh',                      'Mn<I/sig(I)>\nHigh',                           'TEXT'],
            ['DataProcessingCompletenessOverall',           'Completeness\nOverall',                        'TEXT'],
            ['DataProcessingCompletenessLow',               'Completeness\nLow',                            'TEXT'],
            ['DataProcessingCompletenessHigh',              'Completeness\nHigh',                           'TEXT'],
            ['DataProcessingMultiplicityOverall',           'Multiplicity\nOverall',                        'TEXT'],
            ['DataProcessingMultiplicityLow',               'Multiplicity\nLow',                            'TEXT'],
            ['DataProcessingMultiplicityHigh',              'Multiplicity\nHigh',                           'TEXT'],
            ['DataProcessingCChalfOverall',                 'CC(1/2)\nOverall',                             'TEXT'],
            ['DataProcessingCChalfLow',                     'CC(1/2)\nLow',                                 'TEXT'],
            ['DataProcessingCChalfHigh',                    'CC(1/2)\nHigh',                                'TEXT'],

            # the data source is a bit exploding with entries like the ones below,
            # but the many different filenames and folder structures of Diamond autoprocessing makes it necessary
            ['DataProcessingPathToLogfile',                 'DataProcessingPathToLogfile',                  'TEXT'],
            ['DataProcessingPathToMTZfile',                 'DataProcessingPathToMTZfile',                  'TEXT'],
            ['DataProcessingLOGfileName',                   'DataProcessingLOGfileName',                    'TEXT'],
            ['DataProcessingMTZfileName',                   'DataProcessingMTZfileName',                    'TEXT'],
            ['DataProcessingDirectoryOriginal',             'DataProcessingDirectoryOriginal',              'TEXT'],

            ['DataProcessingUniqueReflectionsOverall',      'Unique Reflections\nOverall',                  'TEXT'],
            ['DataProcessingLattice',                       'DataProcessing\nLattice',                      'TEXT'],
            ['DataProcessingPointGroup',                    'DataProcessing\nPointGroup',                   'TEXT'],
            ['DataProcessingUnitCellVolume',                'DataProcessing\nUnit Cell Volume',             'TEXT'],
            ['DataProcessingAlert',                         'DataProcessing\nAlert',                        'TEXT'],

            ['DataProcessingRcryst',                        'DataProcessing\nRcryst',                       'TEXT'],
            ['DataProcessingRfree',                         'DataProcessing\nRfree',                        'TEXT'],
            ['DataProcessingPathToDimplePDBfile',           'DataProcessingPathToDimplePDBfile',            'TEXT'],
            ['DataProcessingPathToDimpleMTZfile',           'DataProcessingPathToDimpleMTZfile',            'TEXT'],

            ['DimpleResolutionHigh',                        'Dimple\nResolution High',                      'TEXT'],
            ['DimpleRcryst',                                'Dimple\nRcryst',                               'TEXT'],
            ['DimpleRfree',                                 'Dimple\nRfree',                                'TEXT'],
            ['DimplePathToPDB',                             'Dimple\nPath to PDB file',                     'TEXT'],
            ['DimplePathToMTZ',                             'Dimple\nPath to MTZ file',                     'TEXT'],

            ['PANDDApath',                           'PANDDApath',                             'TEXT'],

            ['PANDDA_site_A_site_index',                 'PANDDA_site_A_site_index',                      'TEXT'],
            ['PANDDA_site_A_name',                  'PANDDA_site_A_name',                       'TEXT'],
            ['PANDDA_site_A_comment',               'PANDDA_site_A_comment',                    'TEXT'],
            ['PANDDA_site_A_event_index',                 'PANDDA_site_A_event_index',                      'TEXT'],
            ['PANDDA_site_A_event_comment',               'PANDDA_site_A_event_comment',                    'TEXT'],
            ['PANDDA_site_A_confidence',            'PANDDA_site_A_confidence',                 'TEXT'],
            ['PANDDA_site_A_ligand_placed',         'PANDDA_site_A_ligand_placed',              'TEXT'],
            ['PANDDA_site_A_viewed',                'PANDDA_site_A_viewed',                     'TEXT'],
            ['PANDDA_site_A_interesting',           'PANDDA_site_A_interesting',                'TEXT'],
            ['PANDDA_site_A_z_peak',                'PANDDA_site_A_z_peak',                     'TEXT'],
            ['PANDDA_site_A_x',           'PANDDA_site_A_x',                'TEXT'],
            ['PANDDA_site_A_y',           'PANDDA_site_A_y',                'TEXT'],
            ['PANDDA_site_A_z',           'PANDDA_site_A_z',                'TEXT'],
            ['PANDDA_site_A_ligand_id',           'PANDDA_site_A_ligand_id',                'TEXT'],
            ['PANDDA_site_A_event_map',           'PANDDA_site_A_event_map',                'TEXT'],

            ['PANDDA_site_B_site_index',                 'PANDDA_site_B_site_index',                      'TEXT'],
            ['PANDDA_site_B_name',                  'PANDDA_site_B_name',                       'TEXT'],
            ['PANDDA_site_B_comment',               'PANDDA_site_B_comment',                    'TEXT'],
            ['PANDDA_site_B_event_index',                 'PANDDA_site_B_event_index',                      'TEXT'],
            ['PANDDA_site_B_event_comment',               'PANDDA_site_B_event_comment',                    'TEXT'],
            ['PANDDA_site_B_confidence',            'PANDDA_site_B_confidence',                 'TEXT'],
            ['PANDDA_site_B_ligand_placed',         'PANDDA_site_B_ligand_placed',              'TEXT'],
            ['PANDDA_site_B_viewed',                'PANDDA_site_B_viewed',                     'TEXT'],
            ['PANDDA_site_B_interesting',           'PANDDA_site_B_interesting',                'TEXT'],
            ['PANDDA_site_B_z_peak',                'PANDDA_site_B_z_peak',                     'TEXT'],
            ['PANDDA_site_B_x',           'PANDDA_site_B_x',                'TEXT'],
            ['PANDDA_site_B_y',           'PANDDA_site_B_y',                'TEXT'],
            ['PANDDA_site_B_z',           'PANDDA_site_B_z',                'TEXT'],
            ['PANDDA_site_B_ligand_id',           'PANDDA_site_B_ligand_id',                'TEXT'],
            ['PANDDA_site_B_event_map',           'PANDDA_site_B_event_map',                'TEXT'],

            ['PANDDA_site_C_site_index',                 'PANDDA_site_C_site_index',                      'TEXT'],
            ['PANDDA_site_C_name',                  'PANDDA_site_C_name',                       'TEXT'],
            ['PANDDA_site_C_comment',               'PANDDA_site_C_comment',                    'TEXT'],
            ['PANDDA_site_C_event_index',                 'PANDDA_site_C_event_index',                      'TEXT'],
            ['PANDDA_site_C_event_comment',               'PANDDA_site_C_event_comment',                    'TEXT'],
            ['PANDDA_site_C_confidence',            'PANDDA_site_C_confidence',                 'TEXT'],
            ['PANDDA_site_C_ligand_placed',         'PANDDA_site_C_ligand_placed',              'TEXT'],
            ['PANDDA_site_C_viewed',                'PANDDA_site_C_viewed',                     'TEXT'],
            ['PANDDA_site_C_interesting',           'PANDDA_site_C_interesting',                'TEXT'],
            ['PANDDA_site_C_z_peak',                'PANDDA_site_C_z_peak',                     'TEXT'],
            ['PANDDA_site_C_x',           'PANDDA_site_C_x',                'TEXT'],
            ['PANDDA_site_C_y',           'PANDDA_site_C_y',                'TEXT'],
            ['PANDDA_site_C_z',           'PANDDA_site_C_z',                'TEXT'],
            ['PANDDA_site_C_ligand_id',           'PANDDA_site_C_ligand_id',                'TEXT'],
            ['PANDDA_site_C_event_map',           'PANDDA_site_C_event_map',                'TEXT'],

            ['PANDDA_site_D_site_index',                 'PANDDA_site_D_site_index',                      'TEXT'],
            ['PANDDA_site_D_name',                  'PANDDA_site_D_name',                       'TEXT'],
            ['PANDDA_site_D_comment',               'PANDDA_site_D_comment',                    'TEXT'],
            ['PANDDA_site_D_event_index',                 'PANDDA_site_D_event_index',                      'TEXT'],
            ['PANDDA_site_D_event_comment',               'PANDDA_site_D_event_comment',                    'TEXT'],
            ['PANDDA_site_D_confidence',            'PANDDA_site_D_confidence',                 'TEXT'],
            ['PANDDA_site_D_ligand_placed',         'PANDDA_site_D_ligand_placed',              'TEXT'],
            ['PANDDA_site_D_viewed',                'PANDDA_site_D_viewed',                     'TEXT'],
            ['PANDDA_site_D_interesting',           'PANDDA_site_D_interesting',                'TEXT'],
            ['PANDDA_site_D_z_peak',                'PANDDA_site_D_z_peak',                     'TEXT'],
            ['PANDDA_site_D_x',           'PANDDA_site_D_x',                'TEXT'],
            ['PANDDA_site_D_y',           'PANDDA_site_D_y',                'TEXT'],
            ['PANDDA_site_D_z',           'PANDDA_site_D_z',                'TEXT'],
            ['PANDDA_site_D_ligand_id',           'PANDDA_site_D_ligand_id',                'TEXT'],
            ['PANDDA_site_D_event_map',           'PANDDA_site_D_event_map',                'TEXT'],

            ['PANDDA_site_E_site_index',                 'PANDDA_site_E_site_index',                      'TEXT'],
            ['PANDDA_site_E_name',                  'PANDDA_site_E_name',                       'TEXT'],
            ['PANDDA_site_E_comment',               'PANDDA_site_E_comment',                    'TEXT'],
            ['PANDDA_site_E_event_index',                 'PANDDA_site_E_event_index',                      'TEXT'],
            ['PANDDA_site_E_event_comment',               'PANDDA_site_E_event_comment',                    'TEXT'],
            ['PANDDA_site_E_confidence',            'PANDDA_site_E_confidence',                 'TEXT'],
            ['PANDDA_site_E_ligand_placed',         'PANDDA_site_E_ligand_placed',              'TEXT'],
            ['PANDDA_site_E_viewed',                'PANDDA_site_E_viewed',                     'TEXT'],
            ['PANDDA_site_E_interesting',           'PANDDA_site_E_interesting',                'TEXT'],
            ['PANDDA_site_E_z_peak',                'PANDDA_site_E_z_peak',                     'TEXT'],
            ['PANDDA_site_E_x',           'PANDDA_site_E_x',                'TEXT'],
            ['PANDDA_site_E_y',           'PANDDA_site_E_y',                'TEXT'],
            ['PANDDA_site_E_z',           'PANDDA_site_E_z',                'TEXT'],
            ['PANDDA_site_E_ligand_id',           'PANDDA_site_E_ligand_id',                'TEXT'],
            ['PANDDA_site_E_event_map',           'PANDDA_site_E_event_map',                'TEXT'],

            ['PANDDA_site_F_site_index',                 'PANDDA_site_F_site_index',                      'TEXT'],
            ['PANDDA_site_F_name',                  'PANDDA_site_F_name',                       'TEXT'],
            ['PANDDA_site_F_comment',               'PANDDA_site_F_comment',                    'TEXT'],
            ['PANDDA_site_F_event_index',                 'PANDDA_site_F_event_index',                      'TEXT'],
            ['PANDDA_site_F_event_comment',               'PANDDA_site_F_event_comment',                    'TEXT'],
            ['PANDDA_site_F_confidence',            'PANDDA_site_F_confidence',                 'TEXT'],
            ['PANDDA_site_F_ligand_placed',         'PANDDA_site_F_ligand_placed',              'TEXT'],
            ['PANDDA_site_F_viewed',                'PANDDA_site_F_viewed',                     'TEXT'],
            ['PANDDA_site_F_interesting',           'PANDDA_site_F_interesting',                'TEXT'],
            ['PANDDA_site_F_z_peak',                'PANDDA_site_F_z_peak',                     'TEXT'],
            ['PANDDA_site_F_x',           'PANDDA_site_F_x',                'TEXT'],
            ['PANDDA_site_F_y',           'PANDDA_site_F_y',                'TEXT'],
            ['PANDDA_site_F_z',           'PANDDA_site_F_z',                'TEXT'],
            ['PANDDA_site_F_ligand_id',           'PANDDA_site_F_ligand_id',                'TEXT'],
            ['PANDDA_site_F_event_map',           'PANDDA_site_F_event_map',                'TEXT'],

            ['PANDDA_site_G_site_index',                 'PANDDA_site_G_site_index',                      'TEXT'],
            ['PANDDA_site_G_name',                  'PANDDA_site_G_name',                       'TEXT'],
            ['PANDDA_site_G_comment',               'PANDDA_site_G_comment',                    'TEXT'],
            ['PANDDA_site_G_event_index',                 'PANDDA_site_G_event_index',                      'TEXT'],
            ['PANDDA_site_G_event_comment',               'PANDDA_site_G_event_comment',                    'TEXT'],
            ['PANDDA_site_G_confidence',            'PANDDA_site_G_confidence',                 'TEXT'],
            ['PANDDA_site_G_ligand_placed',         'PANDDA_site_G_ligand_placed',              'TEXT'],
            ['PANDDA_site_G_viewed',                'PANDDA_site_G_viewed',                     'TEXT'],
            ['PANDDA_site_G_interesting',           'PANDDA_site_G_interesting',                'TEXT'],
            ['PANDDA_site_G_z_peak',                'PANDDA_site_G_z_peak',                     'TEXT'],
            ['PANDDA_site_G_x',           'PANDDA_site_G_x',                'TEXT'],
            ['PANDDA_site_G_y',           'PANDDA_site_G_y',                'TEXT'],
            ['PANDDA_site_G_z',           'PANDDA_site_G_z',                'TEXT'],
            ['PANDDA_site_G_ligand_id',           'PANDDA_site_G_ligand_id',                'TEXT'],
            ['PANDDA_site_G_event_map',           'PANDDA_site_G_event_map',                'TEXT'],

            ['PANDDA_site_H_site_index',                 'PANDDA_site_H_site_index',                      'TEXT'],
            ['PANDDA_site_H_name',                  'PANDDA_site_H_name',                       'TEXT'],
            ['PANDDA_site_H_comment',               'PANDDA_site_H_comment',                    'TEXT'],
            ['PANDDA_site_H_event_index',                 'PANDDA_site_H_event_index',                      'TEXT'],
            ['PANDDA_site_H_event_comment',               'PANDDA_site_H_event_comment',                    'TEXT'],
            ['PANDDA_site_H_confidence',            'PANDDA_site_H_confidence',                 'TEXT'],
            ['PANDDA_site_H_ligand_placed',         'PANDDA_site_H_ligand_placed',              'TEXT'],
            ['PANDDA_site_H_viewed',                'PANDDA_site_H_viewed',                     'TEXT'],
            ['PANDDA_site_H_interesting',           'PANDDA_site_H_interesting',                'TEXT'],
            ['PANDDA_site_H_z_peak',                'PANDDA_site_H_z_peak',                     'TEXT'],
            ['PANDDA_site_H_x',           'PANDDA_site_H_x',                'TEXT'],
            ['PANDDA_site_H_y',           'PANDDA_site_H_y',                'TEXT'],
            ['PANDDA_site_H_z',           'PANDDA_site_H_z',                'TEXT'],
            ['PANDDA_site_H_ligand_id',           'PANDDA_site_H_ligand_id',                'TEXT'],
            ['PANDDA_site_H_event_map',           'PANDDA_site_H_event_map',                'TEXT'],

            ['PANDDA_site_I_site_index',                 'PANDDA_site_I_site_index',                      'TEXT'],
            ['PANDDA_site_I_name',                  'PANDDA_site_I_name',                       'TEXT'],
            ['PANDDA_site_I_comment',               'PANDDA_site_I_comment',                    'TEXT'],
            ['PANDDA_site_I_event_index',                 'PANDDA_site_I_event_index',                      'TEXT'],
            ['PANDDA_site_I_event_comment',               'PANDDA_site_I_event_comment',                    'TEXT'],
            ['PANDDA_site_I_confidence',            'PANDDA_site_I_confidence',                 'TEXT'],
            ['PANDDA_site_I_ligand_placed',         'PANDDA_site_I_ligand_placed',              'TEXT'],
            ['PANDDA_site_I_viewed',                'PANDDA_site_I_viewed',                     'TEXT'],
            ['PANDDA_site_I_interesting',           'PANDDA_site_I_interesting',                'TEXT'],
            ['PANDDA_site_I_z_peak',                'PANDDA_site_I_z_peak',                     'TEXT'],
            ['PANDDA_site_I_x',           'PANDDA_site_I_x',                'TEXT'],
            ['PANDDA_site_I_y',           'PANDDA_site_I_y',                'TEXT'],
            ['PANDDA_site_I_z',           'PANDDA_site_I_z',                'TEXT'],
            ['PANDDA_site_I_ligand_id',           'PANDDA_site_I_ligand_id',                'TEXT'],
            ['PANDDA_site_I_event_map',           'PANDDA_site_I_event_map',                'TEXT'],

            ['PANDDA_site_J_site_index',                 'PANDDA_site_J_site_index',                      'TEXT'],
            ['PANDDA_site_J_name',                  'PANDDA_site_J_name',                       'TEXT'],
            ['PANDDA_site_J_comment',               'PANDDA_site_J_comment',                    'TEXT'],
            ['PANDDA_site_J_event_index',                 'PANDDA_site_J_event_index',                      'TEXT'],
            ['PANDDA_site_J_event_comment',               'PANDDA_site_J_event_comment',                    'TEXT'],
            ['PANDDA_site_J_confidence',            'PANDDA_site_J_confidence',                 'TEXT'],
            ['PANDDA_site_J_ligand_placed',         'PANDDA_site_J_ligand_placed',              'TEXT'],
            ['PANDDA_site_J_viewed',                'PANDDA_site_J_viewed',                     'TEXT'],
            ['PANDDA_site_J_interesting',           'PANDDA_site_J_interesting',                'TEXT'],
            ['PANDDA_site_J_z_peak',                'PANDDA_site_J_z_peak',                     'TEXT'],
            ['PANDDA_site_J_x',           'PANDDA_site_J_x',                'TEXT'],
            ['PANDDA_site_J_y',           'PANDDA_site_J_y',                'TEXT'],
            ['PANDDA_site_J_z',           'PANDDA_site_J_z',                'TEXT'],
            ['PANDDA_site_J_ligand_id',           'PANDDA_site_J_ligand_id',                'TEXT'],
            ['PANDDA_site_J_event_map',           'PANDDA_site_J_event_map',                'TEXT'],

            ['PANDDA_site_K_site_index',                 'PANDDA_site_K_site_index',                      'TEXT'],
            ['PANDDA_site_K_name',                  'PANDDA_site_K_name',                       'TEXT'],
            ['PANDDA_site_K_comment',               'PANDDA_site_K_comment',                    'TEXT'],
            ['PANDDA_site_K_event_index',                 'PANDDA_site_K_event_index',                      'TEXT'],
            ['PANDDA_site_K_event_comment',               'PANDDA_site_K_event_comment',                    'TEXT'],
            ['PANDDA_site_K_confidence',            'PANDDA_site_K_confidence',                 'TEXT'],
            ['PANDDA_site_K_ligand_placed',         'PANDDA_site_K_ligand_placed',              'TEXT'],
            ['PANDDA_site_K_viewed',                'PANDDA_site_K_viewed',                     'TEXT'],
            ['PANDDA_site_K_interesting',           'PANDDA_site_K_interesting',                'TEXT'],
            ['PANDDA_site_K_z_peak',                'PANDDA_site_K_z_peak',                     'TEXT'],
            ['PANDDA_site_K_x',           'PANDDA_site_K_x',                'TEXT'],
            ['PANDDA_site_K_y',           'PANDDA_site_K_y',                'TEXT'],
            ['PANDDA_site_K_z',           'PANDDA_site_K_z',                'TEXT'],
            ['PANDDA_site_K_ligand_id',           'PANDDA_site_K_ligand_id',                'TEXT'],
            ['PANDDA_site_K_event_map',           'PANDDA_site_K_event_map',                'TEXT'],

            ['PANDDA_site_L_site_index',                 'PANDDA_site_L_site_index',                      'TEXT'],
            ['PANDDA_site_L_name',                  'PANDDA_site_L_name',                       'TEXT'],
            ['PANDDA_site_L_comment',               'PANDDA_site_L_comment',                    'TEXT'],
            ['PANDDA_site_L_event_index',                 'PANDDA_site_L_event_index',                      'TEXT'],
            ['PANDDA_site_L_event_comment',               'PANDDA_site_L_event_comment',                    'TEXT'],
            ['PANDDA_site_L_confidence',            'PANDDA_site_L_confidence',                 'TEXT'],
            ['PANDDA_site_L_ligand_placed',         'PANDDA_site_L_ligand_placed',              'TEXT'],
            ['PANDDA_site_L_viewed',                'PANDDA_site_L_viewed',                     'TEXT'],
            ['PANDDA_site_L_interesting',           'PANDDA_site_L_interesting',                'TEXT'],
            ['PANDDA_site_L_z_peak',                'PANDDA_site_L_z_peak',                     'TEXT'],
            ['PANDDA_site_L_x',           'PANDDA_site_L_x',                'TEXT'],
            ['PANDDA_site_L_y',           'PANDDA_site_L_y',                'TEXT'],
            ['PANDDA_site_L_z',           'PANDDA_site_L_z',                'TEXT'],
            ['PANDDA_site_L_ligand_id',           'PANDDA_site_L_ligand_id',                'TEXT'],
            ['PANDDA_site_L_event_map',           'PANDDA_site_L_event_map',                'TEXT'],


            ['RefinementRcryst',                     'Refinement\nRcryst',                        'TEXT'],
            ['RefinementRfree',                      'Refinement\nRfree',                         'TEXT'],
            ['RefinementSpaceGroup',                      'Refinement\nSpace Group',                         'TEXT'],
            ['RefinementLigandCC',                   'RefinementLigandCC',                      'TEXT'],
            ['RefinementRmsdBonds',                  'RefinementRmsdBonds',                     'TEXT'],
            ['RefinementRmsdAngles',                 'RefinementRmsdAngles',                    'TEXT'],
            ['RefinementOutcome',                    'Refinement\nOutcome',                       'TEXT'],
            ['RefinementMTZfree',                    'RefinementMTZfree',                       'TEXT'],
            ['RefinementCIF',                        'RefinementCIF',                           'TEXT'],
            ['RefinementPDB_latest',                 'RefinementPDB_latest',                    'TEXT'],
            ['RefinementMTZ_latest',                 'RefinementMTZ_latest',                    'TEXT'],
            ['RefinementComment',                    'RefinementComment',                       'TEXT'],
            ['RefinementPathToRefinementFolder',     'RefinementPathToRefinementFolder',        'TEXT'],
            ['RefinementLigandConfidence',           'Ligand\nConfidence',                      'TEXT'],

            ['AssayIC50',                            'AssayIC50',                               'TEXT'],
            ['LastUpdated',                            'LastUpdated',                               'TEXT']
        ]

    def get_empty_db_dict(self):
        db_dict={}
        for column in self.column_list:
            if column[0] != 'ID':
                db_dict[column[0]]=''
        return db_dict

    def create_missing_columns(self):
        existing_columns=[]
        connect=sqlite3.connect(self.data_source_file)
        connect.row_factory = sqlite3.Row
        cursor = connect.cursor()
        cursor.execute("SELECT * FROM mainTable")
        for column in cursor.description:
            existing_columns.append(column[0])
        for column in self.column_list:
            if column[0] not in existing_columns:
                cursor.execute("alter table mainTable add column '"+column[0]+"' '"+column[2]+"'")
                connect.commit()

    def return_column_list(self):
        return self.column_list


    def create_empty_data_source_file(self):

        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        with connect:
            cursor = connect.cursor()
            cursor.execute("CREATE TABLE mainTable("+self.column_list[0][0]+' '+self.column_list[0][2]+")")
            for i in range(1,len(self.column_list)):
                cursor.execute("alter table mainTable add column '"+self.column_list[i][0]+"' '"+self.column_list[i][2]+"'")
            connect.commit()

    def get_all_samples_in_data_source_as_list(self):
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute("SELECT CrystalName FROM mainTable")
        existing_samples_in_db=[]
        samples = cursor.fetchall()
        for sample in samples:
            existing_samples_in_db.append(str(sample[0]))
        return existing_samples_in_db

    def execute_statement(self,cmd):
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute(cmd)
        output=cursor.fetchall()
        return output

    def get_db_dict_for_sample(self,sampleID):
        db_dict={}
        header=[]
        data=[]
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute("select * from mainTable where CrystalName='%s';" %sampleID)
        for column in cursor.description:
            header.append(column[0])
        data = cursor.fetchall()
        for n,item in enumerate(data[0]):
            db_dict[header[n]]=str(item)
        return db_dict


    def check_if_sample_exists_in_data_source(self,sampleID):
        sample_exists=False
        existing_samples_in_db=self.get_all_samples_in_data_source_as_list()
        if sampleID in existing_samples_in_db:
            sample_exists=True
        return sample_exists

    def import_csv_file(self,csv_file):
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        available_columns=[]
        cursor.execute("SELECT * FROM mainTable")
        for column in cursor.description:           # only update existing columns in data source
            available_columns.append(column[0])
        with open(csv_file,'rb') as csv_import: # `with` statement available in 2.5+
            # csv.DictReader uses first line in file for column headings by default
            csv_dict = csv.DictReader(csv_import) # comma is default delimiter
            for line in csv_dict:
                sampleID=line['CrystalName']
                if str(sampleID).replace(' ','')=='':
                    continue
                if self.check_if_sample_exists_in_data_source(sampleID):
                    update_string=''
                    for key,value in line.iteritems():
                        if key=='ID' or key=='CrystalName':
                            continue
                        if key not in available_columns:
                            continue
                        # this is how I had it originally, so it would ignore empty fields
                        # the problem is that if the user wants to set all values to Null,
                        # if will ignore it and leave the inital value in the datasource
#                        if not str(value).replace(' ','')=='':  # ignore if nothing in csv field
#                            update_string+=str(key)+'='+"'"+str(value)+"'"+','
#                            print "UPDATE mainTable SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"';"
#                            cursor.execute("UPDATE mainTable SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"';")
                        # now try this instead; not sure what will break now...
                        update_string+=str(key)+'='+"'"+str(value)+"'"+','
                        print "UPDATE mainTable SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"';"
                        cursor.execute("UPDATE mainTable SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"';")
                else:
                    column_string=''
                    value_string=''
                    for key,value in line.iteritems():
                        if key=='ID':
                            continue
                        if key not in available_columns:
                            continue
                        if not str(value).replace(' ','')=='':  # ignore if nothing in csv field
                            value_string+="'"+value+"'"+','
                            column_string+=key+','
#                    print sampleID
#                    print          "INSERT INTO mainTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");"
                    cursor.execute("INSERT INTO mainTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");")

        connect.commit()

    def update_data_source(self,sampleID,data_dict):
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        update_string=''
        for key in data_dict:
            value=data_dict[key]
            if key=='ID' or key=='CrystalName':
                continue
            if not str(value).replace(' ','')=='':  # ignore empty fields
                update_string+=str(key)+'='+"'"+str(value)+"'"+','
        if update_string != '':
            cursor.execute("UPDATE mainTable SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"'")
            connect.commit()

    def update_insert_data_source(self,sampleID,data_dict):
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute('Select CrystalName FROM mainTable')
        available_columns=[]
        cursor.execute("SELECT * FROM mainTable")
        for column in cursor.description:           # only update existing columns in data source
            available_columns.append(column[0])
        if self.check_if_sample_exists_in_data_source(sampleID):
            for key in data_dict:
                value=data_dict[key]
                if key=='ID' or key=='CrystalName':
                    continue
                if not str(value).replace(' ','')=='':  # ignore empty fields
                    update_string=str(key)+'='+"'"+str(value)+"'"
                    cursor.execute("UPDATE mainTable SET "+update_string+" WHERE CrystalName="+"'"+sampleID+"';")
        else:
            column_string='CrystalName'+','
            value_string="'"+sampleID+"'"+','
            for key in data_dict:
                value=data_dict[key]
                if key=='ID':
                    continue
                if key not in available_columns:
                    continue
                if not str(value).replace(' ','')=='':  # ignore if nothing in csv field
                    value_string+="'"+str(value)+"'"+','
                    column_string+=key+','
            cursor.execute("INSERT INTO mainTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");")
        connect.commit()



    def update_insert_not_null_fields_only(self,sampleID,data_dict):
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute('Select CrystalName FROM mainTable')
        available_columns=[]
        cursor.execute("SELECT * FROM mainTable")
        for column in cursor.description:           # only update existing columns in data source
            available_columns.append(column[0])
        if self.check_if_sample_exists_in_data_source(sampleID):
            for key in data_dict:
                value=data_dict[key]
#                print value
                if key=='ID' or key=='CrystalName':
                    continue
                if not str(value).replace(' ','')=='':  # ignore empty fields
                    update_string=str(key)+'='+"'"+str(value)+"'"
#                    cursor.execute("UPDATE mainTable SET "+update_string+" WHERE CrystalName="+"'"+sampleID+"' and "+str(key)+" is null;")
                    cursor.execute("UPDATE mainTable SET "+update_string+" WHERE CrystalName="+"'"+sampleID+"' and ("+str(key)+" is null or "+str(key)+"='');")
        else:
            column_string='CrystalName'+','
            value_string="'"+sampleID+"'"+','
            for key in data_dict:
                value=data_dict[key]
                if key=='ID':
                    continue
                if key not in available_columns:
                    continue
                if not str(value).replace(' ','')=='':  # ignore if nothing in csv field
                    value_string+="'"+str(value)+"'"+','
                    column_string+=key+','
            cursor.execute("INSERT INTO mainTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");")
        connect.commit()

    def get_value_from_field(self,sampleID,column):
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute("SELECT "+column+" FROM  mainTable WHERE CrystalName='"+sampleID+"';")
        return cursor.fetchone()

    def export_to_csv_file(self,csv_file):
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute("SELECT * FROM mainTable")
        header=()
        for column in cursor.description:
            header+=(column[0],)
        rows = cursor.fetchall()
        csvWriter = csv.writer(open(csv_file, "w"))
        csvWriter.writerows([header]+rows)


#    def load_samples_from_data_source(self):
#        header=[]
#        data=[]
#        connect=sqlite3.connect(self.data_source_file)
#        cursor = connect.cursor()
#        cursor.execute("SELECT * FROM mainTable")
#        for column in cursor.description:
#            header.append(column[0])
#        data = cursor.fetchall()
#        return ([header,data])

    def load_samples_from_data_source(self):
        header=[]
        data=[]
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute("SELECT * FROM mainTable")
        for column in cursor.description:
            header.append(column[0])
        data = cursor.fetchall()
        return header,data





    def get_samples_for_coot(self,RefinementOutcome):
        # 0: sampleID
        # 1: compoundID
        # 2: dataset outcome
        # 3: ligand confidence
        # 4: path to refinement
        sample_list_for_coot=[]
        colums_for_coot=[['RefinementOutcome',                  None],
                         ['CrystalName',                        None],
                         ['CompoundCode',                       None],
                         ['RefinementLigandConfidence',         None],
                         ['RefinementPathToRefinementFolder',   None]]
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()

        print 'data source',self.data_source_file
        print "SELECT CrystalName,CompoundCode,RefinementLigandConfidence,RefinementPathToRefinementFolder FROM mainTable WHERE %s;" % RefinementOutcome
        cursor.execute("SELECT CrystalName,CompoundCode,RefinementLigandConfidence,RefinementPathToRefinementFolder FROM mainTable WHERE %s;" % RefinementOutcome )

        sample_list_for_coot= cursor.fetchall()
        return sample_list_for_coot

    def translate_xce_column_list_to_sqlite(self,column_list):
        out_list=[]
        for item in column_list:
            if item.startswith('img'):
                out_list.append([item,item])
                continue
            if item.startswith('Show'):
                out_list.append([item,item])
                continue
            if item.startswith('Run\nDimple'):
                out_list.append([item,item])
                continue
            if item.startswith('Reference\nSpaceGroup'):
                out_list.append([item,item])
                continue
            if item.startswith('Difference\nUC Volume (%)'):
                out_list.append([item,item])
                continue
            if item.startswith('Reference File'):
                out_list.append([item,item])
                continue
            for entry in self.column_list:
                if entry[1]==item:
                    out_list.append([item,entry[0]])
                    break
        return out_list