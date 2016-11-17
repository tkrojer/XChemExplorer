# last edited: 15/11/2016, 15:00

import sqlite3
import os,sys
import csv
from datetime import datetime
import getpass

#sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
#from XChemUtils import parse

class data_source:

    def __init__(self,data_source_file):

        self.data_source_file=data_source_file

        #   [column_name in DB, column_name shown in XCE, SQLite column type (Integer,Text,PKEY...)]
        self.column_list=[
            # SQLite column name                    XCE column name                             SQLite type             display in overview tab
            # from Lab36
            ['ID',                                   'ID',                                      'INTEGER PRIMARY KEY',  0],
            ['LabVisit',                             'LabVisit',                                'TEXT',                 1],
            ['LibraryPlate',                         'LibraryPlate',                            'TEXT',                 0],
            ['SourceWell',                           'SourceWell',                              'TEXT',                 0],
            ['LibraryName',                          'LibraryName',                             'TEXT',                 1],
            ['CompoundSMILES',                       'Smiles',                                  'TEXT',                 1],
            ['CompoundCode',                         'Compound ID',                             'TEXT',                 1],
            ['CrystalPlate',                         'CrystalPlate',                            'TEXT',                 1],
            ['CrystalWell',                          'CrystalWell',                             'TEXT',                 1],
            ['EchoX',                                'EchoX',                                   'TEXT',                 0],
            ['EchoY',                                'EchoY',                                   'TEXT',                 0],
            ['DropVolume',                           'DropVolume',                              'TEXT',                 0],
            ['ProteinName',                          'ProteinName',                             'TEXT',                 1],
            ['BatchNumber',                          'BatchNumber',                             'TEXT',                 0],
            ['CompoundStockConcentration',           'CompoundStockConcentration',              'TEXT',                 0],
            ['CompoundConcentration',                'CompoundConcentration',                   'TEXT',                 1],
            ['SolventFraction',                      'SolventFraction',                         'TEXT',                 1],
            ['SoakTransferVol',                      'SoakTransferVol',                         'TEXT',                 0],
            ['SoakStatus',                           'SoakStatus',                              'TEXT',                 0],
            ['SoakTimestamp',                        'SoakTimestamp',                           'TEXT',                 0],
            ['CryoStockFraction',                    'CryoStockFraction',                       'TEXT',                 0],
            ['CryoFraction',                         'CryoFraction',                            'TEXT',                 0],
            ['CryoWell',                             'CryoWell',                                'TEXT',                 0],
            ['CryoTransferVolume',                   'CryoTransferVolume',                      'TEXT',                 0],
            ['CryoStatus',                           'CryoStatus',                              'TEXT',                 0],
            ['CryoTimestamp',                        'CryoTimestamp',                           'TEXT',                 0],
            ['SoakingTime',                          'SoakingTime',                             'TEXT',                 1],
            ['HarvestStatus',                        'HarvestStatus',                           'TEXT',                 0],
            ['CrystalName',                          'Sample ID',                               'TEXT',                 0],
            ['Puck',                                 'Puck',                                    'TEXT',                 1],
            ['PuckPosition',                         'PuckPosition',                            'TEXT',                 1],
            ['PinBarcode',                           'SoakDB\nBarcode',                         'TEXT',                 1],
            ['MountingResult',                       'MountingResult',                          'TEXT',                 0],
            ['MountingArrivalTime',                  'MountingArrivalTime',                     'TEXT',                 0],
            ['MountedTimestamp',                     'MountedTimestamp',                        'TEXT',                 0],
            ['MountingTime',                         'MountingTime',                            'TEXT',                 0],
            ['ispybStatus',                          'ispybStatus',                             'TEXT',                 0],
            ['DataCollectionVisit',                  'Visit',                                   'TEXT',                 1],

            # from XChemExplorer

            ['ProjectDirectory',                     'ProjectDirectory',                        'TEXT',                 0],

            ['CrystalTag',                           'Tag',                                     'TEXT',                 0],
            ['CrystalFormName',                      'Crystal Form\nName',                      'TEXT',                 0],
            ['CrystalFormSpaceGroup',                'Space\nGroup',                            'TEXT',                 0],
            ['CrystalFormPointGroup',                'Point\nGroup',                            'TEXT',                 0],
            ['CrystalFormA',                         'a',                                       'TEXT',                 0],
            ['CrystalFormB',                         'b',                                       'TEXT',                 0],
            ['CrystalFormC',                         'c',                                       'TEXT',                 0],
            ['CrystalFormAlpha',                     'alpha',                                   'TEXT',                 0],
            ['CrystalFormBeta',                      'beta',                                    'TEXT',                 0],
            ['CrystalFormGamma',                     'gamma',                                   'TEXT',                 0],
            ['CrystalFormVolume',                    'Crystal Form\nVolume',                    'TEXT',                 0],

            ['DataCollectionBeamline',               'Beamline',                                'TEXT',                 0],
            ['DataCollectionDate',                   'Data Collection\nDate',                   'TEXT',                 1],
            ['DataCollectionOutcome',                'DataCollection\nOutcome',                 'TEXT',                 1],
            ['DataCollectionRun',                    'Run',                                     'TEXT',                 0],
            ['DataCollectionComment',                'DataCollection\nComment',                 'TEXT',                 0],
            ['DataCollectionWavelength',             'Wavelength',                              'TEXT',                 0],
            ['DataCollectionPinBarcode',             'GDA\nBarcode',                            'TEXT',                 1],

            ['DataProcessingPathToImageFiles',       'Path to diffraction\nimage files',        'TEXT',                 1],
            ['DataProcessingProgram',                'Program',                                 'TEXT',                 1],
            ['DataProcessingSpaceGroup',             'DataProcessing\nSpaceGroup',              'TEXT',                 1],
            ['DataProcessingUnitCell',               'DataProcessing\nUnitCell',                'TEXT',                 0],
            ['DataProcessingAutoAssigned',           'auto-assigned',                           'TEXT',                 0],

            ['DataProcessingA',                      'DataProcessing\nA',                       'TEXT',                 0],
            ['DataProcessingB',                      'DataProcessing\nB',                       'TEXT',                 0],
            ['DataProcessingC',                      'DataProcessing\nC',                       'TEXT',                 0],
            ['DataProcessingAlpha',                  'DataProcessing\nAlpha',                   'TEXT',                 0],
            ['DataProcessingBeta',                   'DataProcessing\nBeta',                    'TEXT',                 0],
            ['DataProcessingGamma',                  'DataProcessing\nGamma',                   'TEXT',                 0],

            ['DataProcessingResolutionOverall',             'Resolution\nOverall',                          'TEXT',                 0],
            ['DataProcessingResolutionLow',                 'Resolution\nLow',                              'TEXT',                 0],
            ['DataProcessingResolutionLowInnerShell',       'Resolution\nLow (Inner Shell)',                'TEXT',                 0],
            ['DataProcessingResolutionHigh',                'Resolution\nHigh',                             'TEXT',                 1],
            ['DataProcessingResolutionHigh15sigma',         'Resolution\n[Mn<I/sig(I)> = 1.5]',             'TEXT',                 1],
            ['DataProcessingResolutionHighOuterShell',      'Resolution\nHigh (Outer Shell)',               'TEXT',                 0],
            ['DataProcessingRmergeOverall',                 'Rmerge\nOverall',                              'TEXT',                 1],
            ['DataProcessingRmergeLow',                     'Rmerge\nLow',                                  'TEXT',                 1],
            ['DataProcessingRmergeHigh',                    'Rmerge\nHigh',                                 'TEXT',                 1],
            ['DataProcessingIsigOverall',                   'Mn<I/sig(I)>\nOverall',                        'TEXT',                 1],
            ['DataProcessingIsigLow',                       'Mn<I/sig(I)>\nLow',                            'TEXT',                 1],
            ['DataProcessingIsigHigh',                      'Mn<I/sig(I)>\nHigh',                           'TEXT',                 1],
            ['DataProcessingCompletenessOverall',           'Completeness\nOverall',                        'TEXT',                 1],
            ['DataProcessingCompletenessLow',               'Completeness\nLow',                            'TEXT',                 1],
            ['DataProcessingCompletenessHigh',              'Completeness\nHigh',                           'TEXT',                 1],
            ['DataProcessingMultiplicityOverall',           'Multiplicity\nOverall',                        'TEXT',                 1],
            ['DataProcessingMultiplicityLow',               'Multiplicity\nLow',                            'TEXT',                 1],
            ['DataProcessingMultiplicityHigh',              'Multiplicity\nHigh',                           'TEXT',                 1],
            ['DataProcessingCChalfOverall',                 'CC(1/2)\nOverall',                             'TEXT',                 1],
            ['DataProcessingCChalfLow',                     'CC(1/2)\nLow',                                 'TEXT',                 1],
            ['DataProcessingCChalfHigh',                    'CC(1/2)\nHigh',                                'TEXT',                 1],

            # the data source is a bit exploding with entries like the ones below,
            # but the many different filenames and folder structures of Diamond autoprocessing makes it necessary
            ['DataProcessingPathToLogfile',                 'DataProcessingPathToLogfile',                  'TEXT',                 1],
            ['DataProcessingPathToMTZfile',                 'DataProcessingPathToMTZfile',                  'TEXT',                 1],
            ['DataProcessingLOGfileName',                   'DataProcessingLOGfileName',                    'TEXT',                 0],
            ['DataProcessingMTZfileName',                   'DataProcessingMTZfileName',                    'TEXT',                 0],
            ['DataProcessingDirectoryOriginal',             'DataProcessingDirectoryOriginal',              'TEXT',                 0],

            ['DataProcessingUniqueReflectionsOverall',      'Unique Reflections\nOverall',                  'TEXT',                 1],
            ['DataProcessingLattice',                       'DataProcessing\nLattice',                      'TEXT',                 0],
            ['DataProcessingPointGroup',                    'DataProcessing\nPointGroup',                   'TEXT',                 0],
            ['DataProcessingUnitCellVolume',                'DataProcessing\nUnit Cell Volume',             'TEXT',                 0],
            ['DataProcessingAlert',                         'DataProcessing\nAlert',                        'TEXT',                 0],
            ['DataProcessingScore',                         'DataProcessing\nScore',                        'TEXT',                 1],

            ['DataProcessingStatus',                        'DataProcessing\nStatus',                       'TEXT',                 1],

            ['DataProcessingRcryst',                        'DataProcessing\nRcryst',                       'TEXT',                 0],
            ['DataProcessingRfree',                         'DataProcessing\nRfree',                        'TEXT',                 0],
            ['DataProcessingPathToDimplePDBfile',           'DataProcessingPathToDimplePDBfile',            'TEXT',                 0],
            ['DataProcessingPathToDimpleMTZfile',           'DataProcessingPathToDimpleMTZfile',            'TEXT',                 0],
            ['DataProcessingDimpleSuccessful',              'DataProcessingDimpleSuccessful',               'TEXT',                 0],

            ['DimpleResolutionHigh',                        'Dimple\nResolution High',                      'TEXT',                 1],
            ['DimpleRcryst',                                'Dimple\nRcryst',                               'TEXT',                 1],
            ['DimpleRfree',                                 'Dimple\nRfree',                                'TEXT',                 1],
            ['DimplePathToPDB',                             'Dimple\nPath to PDB file',                     'TEXT',                 1],
            ['DimplePathToMTZ',                             'Dimple\nPath to MTZ file',                     'TEXT',                 1],
            ['DimpleReferencePDB',                          'Dimple\nReference PDB',                        'TEXT',                 0],
            ['DimplePANDDAwasRun',                          'PANDDA\nlaunched?',                            'TEXT',                 1],
            ['DimplePANDDAhit',                             'PANDDA\nhit?',                                 'TEXT',                 1],
            ['DimplePANDDAreject',                          'PANDDA\nreject?',                              'TEXT',                 1],
            ['DimplePANDDApath',                            'PANDDA\npath?',                                'TEXT',                 1],
            ['DimpleStatus',                                'Dimple\nStatus',                               'TEXT',                 1],

            ['RefinementResolution',                        'Refinement\nResolution',                       'TEXT',                 1],
            ['RefinementResolutionTL',                      'RefinementResolutionTL',                       'TEXT',                 0],
            ['RefinementRcryst',                            'Refinement\nRcryst',                           'TEXT',                 1],
            ['RefinementRcrystTraficLight',                 'RefinementRcrystTraficLight',                  'TEXT',                 0],
            ['RefinementRfree',                             'Refinement\nRfree',                            'TEXT',                 1],
            ['RefinementRfreeTraficLight',                  'RefinementRfreeTraficLight',                   'TEXT',                 0],
            ['RefinementSpaceGroup',                        'Refinement\nSpace Group',                      'TEXT',                 1],
            ['RefinementLigandCC',                          'RefinementLigandCC',                           'TEXT',                 0],
            ['RefinementRmsdBonds',                         'RefinementRmsdBonds',                          'TEXT',                 1],
            ['RefinementRmsdBondsTL',                       'RefinementRmsdBondsTL',                        'TEXT',                 0],
            ['RefinementRmsdAngles',                        'RefinementRmsdAngles',                         'TEXT',                 1],
            ['RefinementRmsdAnglesTL',                      'RefinementRmsdAnglesTL',                       'TEXT',                 0],
            ['RefinementOutcome',                           'Refinement\nOutcome',                          'TEXT',                 1],
            ['RefinementMTZfree',                           'RefinementMTZfree',                            'TEXT',                 1],
            ['RefinementCIF',                               'RefinementCIF',                                'TEXT',                 1],
            ['RefinementCIFStatus',                         'Compound\nStatus',                             'TEXT',                 1],
            ['RefinementCIFprogram',                        'RefinementCIFprogram',                         'TEXT',                 1],
            ['RefinementPDB_latest',                        'RefinementPDB_latest',                         'TEXT',                 1],
            ['RefinementMTZ_latest',                        'RefinementMTZ_latest',                         'TEXT',                 1],
            ['RefinementMatrixWeight',                      'RefinementMatrixWeight',                       'TEXT',                 0],
            ['RefinementComment',                           'RefinementComment',                            'TEXT',                 0],
            ['RefinementPathToRefinementFolder',            'RefinementPathToRefinementFolder',             'TEXT',                 0],
            ['RefinementLigandConfidence',                  'Ligand\nConfidence',                           'TEXT',                 0],
            ['RefinementLigandBoundConformation',           'RefinementLigandBoundConformation',            'TEXT',                 0],
            ['RefinementBoundConformation',                 'RefinementBoundConformation',                  'TEXT',                 0],
            ['RefinementMolProbityScore',                   'MolProbity Score',                             'TEXT',                 1],
            ['RefinementMolProbityScoreTL',                 'RefinementMolProbityScoreTL',                  'TEXT',                 0],
            ['RefinementRamachandranOutliers',              'Ramachandran\nOutliers',                       'TEXT',                 1],
            ['RefinementRamachandranOutliersTL',            'RefinementRamachandranOutliersTL',             'TEXT',                 0],
            ['RefinementRamachandranFavored',               'Ramachandran\nFavored',                        'TEXT',                 1],
            ['RefinementRamachandranFavoredTL',             'RefinementRamachandranFavoredTL',              'TEXT',                 0],
            ['RefinementStatus',                            'Refinement\nStatus',                           'TEXT',                 1],

            ['Deposition_PDB_ID',                           'Deposition_PDB_ID',                            'TEXT',                 1],
            ['Deposition_Date',                             'Deposition_Date',                              'TEXT',                 1],


            ['AssayIC50',                                   'AssayIC50',                                    'TEXT',                 0],
            ['LastUpdated',                                 'LastUpdated',                                  'TEXT',                 0],
            ['LastUpdated_by',                              'LastUpdated_by',                               'TEXT',                 0]
        ]

        self.pandda_table_columns = [
            ['ID',                                          'ID',                                       'INTEGER PRIMARY KEY'],
            ['CrystalName',                                 'Sample ID',                                'TEXT'],
            ['PANDDApath',                                  'PANDDApath',                               'TEXT'],
            ['PANDDA_site_index',                           'PANDDA_site_index',                        'TEXT'],
            ['PANDDA_site_name',                            'PANDDA_site_name',                         'TEXT'],
            ['PANDDA_site_comment',                         'PANDDA_site_comment',                      'TEXT'],
            ['PANDDA_site_event_index',                     'PANDDA_site_event_index',                  'TEXT'],
            ['PANDDA_site_event_comment',                   'PANDDA_site_event_comment',                'TEXT'],
            ['PANDDA_site_confidence',                      'PANDDA_site_confidence',                   'TEXT'],
            ['PANDDA_site_ligand_placed',                   'PANDDA_site_ligand_placed',                'TEXT'],
            ['PANDDA_site_viewed',                          'PANDDA_site_viewed',                       'TEXT'],
            ['PANDDA_site_interesting',                     'PANDDA_site_interesting',                  'TEXT'],
            ['PANDDA_site_z_peak',                          'PANDDA_site_z_peak',                       'TEXT'],
            ['PANDDA_site_x',                               'PANDDA_site_x',                            'TEXT'],
            ['PANDDA_site_y',                               'PANDDA_site_y',                            'TEXT'],
            ['PANDDA_site_z',                               'PANDDA_site_z',                            'TEXT'],
            ['PANDDA_site_ligand_id',                       'PANDDA_site_ligand_id',                    'TEXT'],

            ['PANDDA_site_ligand_resname',                  'PANDDA_site_ligand_resname',               'TEXT'],
            ['PANDDA_site_ligand_chain',                    'PANDDA_site_ligand_chain',                 'TEXT'],
            ['PANDDA_site_ligand_sequence_number',          'PANDDA_site_ligand_sequence_number',       'TEXT'],
            ['PANDDA_site_ligand_altLoc',                   'PANDDA_site_ligand_altLoc',                'TEXT'],


            ['PANDDA_site_event_map',                       'PANDDA_site_event_map',                    'TEXT'],
            ['PANDDA_site_event_map_mtz',                   'PANDDA_site_event_map_mtz',                'TEXT'],
            ['PANDDA_site_initial_model',                   'PANDDA_site_initial_model',                'TEXT'],
            ['PANDDA_site_initial_mtz',                     'PANDDA_site_initial_mtz',                  'TEXT'],
            ['PANDDA_site_spider_plot',                     'PANDDA_site_spider_plot',                  'TEXT'],
            ['PANDDA_site_occupancy',                       'PANDDA_site_occupancy',                    'TEXT'],
            ['PANDDA_site_B_average',                       'PANDDA_site_B_average',                    'TEXT'],
            ['PANDDA_site_B_ratio_residue_surroundings',    'PANDDA_site_B_ratio_residue_surroundings', 'TEXT'],
            ['PANDDA_site_RSCC',                            'PANDDA_site_RSCC',                         'TEXT'],
            ['PANDDA_site_RSR',                             'PANDDA_site_RSR',                          'TEXT'],
            ['PANDDA_site_RSZD',                            'PANDDA_site_RSZD',                         'TEXT'],
            ['PANDDA_site_rmsd',                            'PANDDA_site_rmsd',                         'TEXT'],
            ['RefinementOutcome',                           'RefinementOutcome',                        'TEXT'],
            ['LastUpdated',                                 'LastUpdated',                              'TEXT'],
            ['LastUpdated_by',                              'LastUpdated_by',                           'TEXT']
            ]

    def columns_not_to_display(self):
        do_not_display = []
        for column in self.column_list:
            if column[3]==0:
                do_not_display.append(column[1])
        return do_not_display
        
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
        # create PANDDA table if not exists
        cursor.execute("create table if not exists panddaTable (ID INTEGER);")
        existing_columns=[]
        cursor.execute("SELECT * FROM panddaTable")
        for column in cursor.description:
            existing_columns.append(column[0])
        for column in self.pandda_table_columns:
            if column[0] not in existing_columns:
                cursor.execute("alter table panddaTable add column '"+column[0]+"' '"+column[2]+"'")
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
        # Don't need to create panddaTable at this point, because table will be created by create_missing_columns
        # which is called the first time a data source in specified in XCE
#        with connect:
#            cursor = connect.cursor()
#            cursor.execute("CREATE TABLE panddaTable("+self.pandda_table_columns[0][0]+' '+self.pandda_table_columns[0][2]+")")
#            for i in range(1,len(self.pandda_table_columns)):
#                cursor.execute("alter table mainTable add column '"+self.pandda_table_columns[i][0]+"' '"+self.pandda_table_columns[i][2]+"'")
#            connect.commit()

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
        connect.commit()
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

    def get_db_pandda_dict_for_sample_and_site(self,sampleID,site_index):
        db_dict={}
        header=[]
        data=[]
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute("select * from panddaTable where CrystalName='%s' and PANDDA_site_index='%s';" %(sampleID,site_index))
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
#                        print "UPDATE mainTable SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"';"
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
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        update_string=''
        for key in data_dict:
            value=data_dict[key]
            if key=='ID' or key=='CrystalName':
                continue
            if not str(value).replace(' ','')=='':  # ignore empty fields
                update_string+=str(key)+'='+"'"+str(value)+"'"+','
            else:
                update_string+=str(key)+' = null,'
        if update_string != '':
#            print "UPDATE mainTable SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"'"
            cursor.execute("UPDATE mainTable SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"'")
            connect.commit()

    def update_panddaTable(self,sampleID,site_index,data_dict):
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        update_string=''
        for key in data_dict:
            value=data_dict[key]
            if key=='ID' or key=='CrystalName' or key=='PANDDA_site_index':
                continue
            if not str(value).replace(' ','')=='':  # ignore empty fields
                update_string+=str(key)+'='+"'"+str(value)+"'"+','
            else:
                update_string+=str(key)+' = null,'
        if update_string != '':
            cursor.execute("UPDATE panddaTable SET "+update_string[:-1]+" WHERE CrystalName='%s' and PANDDA_site_index='%s'" %(sampleID,site_index))
            connect.commit()

    def update_specified_table(self,sampleID,data_dict,table):
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()
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
            cursor.execute("UPDATE "+table+" SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"'")
            connect.commit()

    def update_insert_data_source(self,sampleID,data_dict):
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()
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



    def update_insert_panddaTable(self,sampleID,data_dict):
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute('Select CrystalName,PANDDA_site_index FROM panddaTable')
#        available_columns=[]
#        cursor.execute("SELECT * FROM panddaTable")
#        for column in cursor.description:           # only update existing columns in data source
#            available_columns.append(column[0])
        samples_sites_in_table=[]
        tmp=cursor.fetchall()
        for item in tmp:
            line=[x.encode('UTF8') for x in list(item)]
            samples_sites_in_table.append(line)

        found_sample_site=False
        for entry in samples_sites_in_table:
            if entry[0]==sampleID and entry[1]==data_dict['PANDDA_site_index']:
                found_sample_site=True

        if found_sample_site:
            for key in data_dict:
                value=data_dict[key]
                if key=='ID' or key=='CrystalName' or key=='PANDDA_site_index':
                    continue
                if not str(value).replace(' ','')=='':  # ignore empty fields
                    update_string=str(key)+'='+"'"+str(value)+"'"
#                    print "UPDATE panddaTable SET "+update_string+" WHERE CrystalName="+"'"+sampleID+"' and PANDDA_site_index is '"+data_dict['PANDDA_site_index']+"';"
                    cursor.execute("UPDATE panddaTable SET "+update_string+" WHERE CrystalName="+"'"+sampleID+"' and PANDDA_site_index is '"+data_dict['PANDDA_site_index']+"';")
        else:
            column_string=''
            value_string=''
            for key in data_dict:
                value=data_dict[key]
                if key=='ID':
                    continue
#                if key not in available_columns:
#                    continue
                if not str(value).replace(' ','')=='':  # ignore if nothing in csv field
                    value_string+="'"+str(value)+"'"+','
                    column_string+=key+','
            print "INSERT INTO panddaTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");"
            cursor.execute("INSERT INTO panddaTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");")
        connect.commit()





    def update_insert_not_null_fields_only(self,sampleID,data_dict):
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()
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





#    def get_samples_for_coot(self,RefinementOutcome,pandda_site):
#        sample_list_for_coot=[]
#        connect=sqlite3.connect(self.data_source_file)
#        cursor = connect.cursor()
#
#        if RefinementOutcome=='0 - All Datasets':
#            outcome = " not null "
#        else:
#            outcome = " '%s' " %RefinementOutcome
#
#        if int(pandda_site) > 0:
#            sqlite = (
#                "select"
#                " mainTable.CrystalName,"
#                " mainTable.CompoundCode,"
#                " mainTable.RefinementCIF,"
#                " mainTable.RefinementMTZfree,"
#                " mainTable.RefinementPathToRefinementFolder,"
#                " panddaTable.RefinementOutcome, "
#                " panddaTable.PANDDA_site_event_map,"
#                " panddaTable.PANDDA_site_confidence,"
#                " panddaTable.PANDDA_site_x,"
#                " panddaTable.PANDDA_site_y,"
#                " panddaTable.PANDDA_site_z,"
#                " panddaTable.PANDDA_site_initial_model,"
#                " panddaTable.PANDDA_site_initial_mtz,"
#                " panddaTable.PANDDA_site_spider_plot, "
#                " panddaTable.PANDDA_site_index "
#                "from mainTable inner join panddaTable on mainTable.CrystalName = panddaTable.CrystalName "
#                "where panddaTable.PANDDA_site_index is '%s'" %pandda_site+
#                " and panddaTable.PANDDA_site_ligand_placed is 'True'"
#                " and panddaTable.RefinementOutcome is %s;" %outcome
#                )
#
#        else:
#            sqlite = (
#                "select"
#                " CrystalName,"
#                " CompoundCode,"
#                " RefinementCIF,"
#                " RefinementMTZfree,"
#                " RefinementPathToRefinementFolder,"
#                " RefinementOutcome "
#                "from mainTable "
#                "where RefinementOutcome is %s;" %outcome
#                )
#
#        cursor.execute(sqlite)
#
#        tmp = cursor.fetchall()
#        for item in tmp:
#            tmpx=[]
#            for i in list(item):
#                if i==None:
#                    tmpx.append('None')
#                else:
#                    tmpx.append(i)
#            line=[x.encode('UTF8') for x in tmpx]
#            sample_list_for_coot.append(line)
#
#        return sample_list_for_coot

    def get_pandda_info_for_coot(self,xtalID,pandda_site):
        pandda_info_for_coot=[]
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()

        sqlite = (
            "select"
            " PANDDA_site_event_map,"
            " PANDDA_site_x,"
            " PANDDA_site_y,"
            " PANDDA_site_z,"
            " PANDDA_site_spider_plot "
            "from panddaTable  "
            "where "
            " CrystalName is '%s'" %xtalID+
            " and PANDDA_site_index is '%s';" %pandda_site
            )

        cursor.execute(sqlite)

        tmp = cursor.fetchall()
        for item in tmp:
            tmpx=[]
            for i in list(item):
                if i==None:
                    tmpx.append('None')
                else:
                    tmpx.append(i)
            line=[x.encode('UTF8') for x in tmpx]
            pandda_info_for_coot.append(line)

        return pandda_info_for_coot



    def get_todo_list_for_coot(self,RefinementOutcome,pandda_site):
        sample_list_for_coot=[]
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()

        if RefinementOutcome=='0 - All Datasets':
            outcome = " not null "
        else:
            outcome = " '%s' " %RefinementOutcome

        if int(pandda_site) > 0:
#            sqlite = (
#                "select"
#                " mainTable.CrystalName,"
#                " mainTable.CompoundCode,"
#                " mainTable.RefinementCIF,"
#                " mainTable.RefinementMTZfree,"
#                " mainTable.RefinementPathToRefinementFolder,"
#                " panddaTable.RefinementOutcome, "
#                " panddaTable.PANDDA_site_confidence "
#                "from mainTable inner join panddaTable on mainTable.CrystalName = panddaTable.CrystalName "
#                "where panddaTable.PANDDA_site_index is '%s'" %pandda_site+
#                " and panddaTable.PANDDA_site_ligand_placed is 'True'"
#                " and panddaTable.RefinementOutcome is %s;" %outcome
#                )
            sqlite = (
                "select"
                " mainTable.CrystalName,"
                " mainTable.CompoundCode,"
                " mainTable.RefinementCIF,"
                " mainTable.RefinementMTZfree,"
                " mainTable.RefinementPathToRefinementFolder,"
                " panddaTable.RefinementOutcome, "
                " panddaTable.PANDDA_site_confidence "
                "from mainTable inner join panddaTable on mainTable.CrystalName = panddaTable.CrystalName "
                "where panddaTable.PANDDA_site_index is '%s'" %pandda_site+
                " and panddaTable.PANDDA_site_ligand_placed is 'True'"
                " and panddaTable.RefinementOutcome like "+outcome.split()[0]+"%';"
                )
        else:
            sqlite = (
                "select"
                " CrystalName,"
                " CompoundCode,"
                " RefinementCIF,"
                " RefinementMTZfree,"
                " RefinementPathToRefinementFolder,"
                " RefinementOutcome,"
                " RefinementLigandConfidence "
                "from mainTable "
                "where RefinementOutcome is %s;" %outcome
                )
#            sqlite = (
#                "select"
#                " CrystalName,"
#                " CompoundCode,"
#                " RefinementCIF,"
#                " RefinementMTZfree,"
#                " RefinementPathToRefinementFolder,"
#                " RefinementOutcome,"
#                " RefinementLigandConfidence "
#                "from mainTable "
#                "where RefinementOutcome like "+outcome.split()[0]+"%';"
#                )

#        print sqlite
        cursor.execute(sqlite)

        tmp = cursor.fetchall()
        for item in tmp:
            tmpx=[]
            for i in list(item):
                if i==None:
                    tmpx.append('None')
                else:
                    tmpx.append(i)
            line=[x.encode('UTF8') for x in tmpx]
            sample_list_for_coot.append(line)

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
            if item.startswith('Select'):
                out_list.append([item,item])
                continue
            if item.startswith('Run\nxia2'):
                out_list.append([item,item])
                continue
            if item.startswith('Dataset ID'):
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
            if item.startswith('PanDDA site details'):
                out_list.append([item,item])
                continue
            for entry in self.column_list:
                if entry[1]==item:
                    out_list.append([item,entry[0]])
                    break
        return out_list

    def get_list_of_pandda_sites_for_coot(self):
        site_list=[ ['0','any site'] ]
        sqlite = (
            'select distinct'
            ' panddaTable.PANDDA_site_index,'
            ' panddaTable.PANDDA_site_name '
            'from panddaTable '
            'order by cast (panddaTable.PANDDA_site_index as integer) ASC;'
                  )

        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute(sqlite)
        tmp=cursor.fetchall()
        for item in tmp:
            line=[x.encode('UTF8') for x in list(item)]
            site_list.append(line)

        return site_list


