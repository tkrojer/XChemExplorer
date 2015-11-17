import sqlite3

class data_source:

    def __init__(self,data_source_file):

        self.data_source_file=data_source_file

        self.data_source_type=None
        if self.data_source_file.endswith('.sqlite'):
            self.data_source_type='sqlite'
        if self.data_source_file.endswith('.csv'):
            self.data_source_type='csv'


        self.column_list=[
            # from Lab36
            'ID',
            'LabVisit',
            'LibraryPlate',
            'SourceWell',
            'LibraryName',
            'CompoundSMILES',
            'CompoundCode',
            'CrystalPlate',
            'CrystalWell',
            'EchoX',
            'EchoY',
            'DropVolume',
            'ProteinName',
            'BatchNumber',
            'CompoundStockConcentration',
            'CompoundConcentration',
            'SolventFraction',
            'SoakTransferVol',
            'SoakStatus',
            'SoakTimestamp',
            'CryoStockFraction',
            'CryoFraction',
            'CryoWell',
            'CryoTransferVolume',
            'CryoStatus',
            'CryoTimestamp',
            'SoakingTime',
            'HarvestStatus',
            'CrystalName',
            'Puck',
            'PuckPosition',
            'PinBarcode',
            'MountingResult',
            'MountingArrivalTime',
            'MountedTimestamp',
            'MountingTime',
            'ispybStatus',
            'DataCollectionVisit',
            # from XChemExplorer


#            'SampleID',
#            'Visit',
            'DataCollectionBeamline'
            'DataCollectionDate',
            'DataCollectionOutcome',
            'DataCollectionRun',



            'DataProcessingProgram',
            'DataProcessingSpaceGroup',
            'DataProcessingUnitCell',
            'DataProcessingResolutionOverall',
            'DataProcessingResolutionLow',
            'DataProcessingResolutionHigh',
            'DataProcessingRmergeOverall',
            'DataProcessingRmergeLow',
            'DataProcessingRmergeHigh',
            'DataProcessingIsigOverall',
            'DataProcessingIsigLow',
            'DataProcessingIsigHigh',
            'DataProcessingCompletenessOverall',
            'DataProcessingCompletenessLow',
            'DataProcessingCompletenessHigh',
            'DataProcessingMultiplicityOverall',
            'DataProcessingMultiplicityLow',
            'DataProcessingMultiplicityHigh',
            'DataProcessingPathToLogfile'

            'RefinementRcryst',
            'RefinementRfree',
            'RefinementLigandCC',
            'RefinementRmsdBonds',
            'RefinementRmsdAngles',
            'RefinementOutcome'

        ]

        # 1. check data source if all columns exist
        # 2. create missing columns

        if self.data_source_type=='sqlite':
            sample_dict={}
con = sqlite3.connect('soakDBDataFile.sqlite')

con.row_factory = sqlite3.Row
cur = con.cursor()
cur.execute("SELECT * FROM mainTable")

#r=c.fetchone()
for column in cur.description:
    print column[0]


#rows = cur.fetchall()

#for row in rows:
#    print row["CompoundSMILES"]
#    print row
#    sample_dict[row["ID"]]=[row["CompoundCode"],row["CompoundSMILES"]]
#print sample_dict
#for key in sample_dict:
#    print key,sample_dict[key]




    def read_data_source_for_coot(self):
        print 'hallo'
    def update_data_source_from_coot(self):
        print 'hallo'
    def update_data_source(self):
        print 'hallo'
    def convert_sqlite_to_csv(self):
        print 'hallo'
    def convert_csv_to_sqlite(self):
        print 'hallo'

