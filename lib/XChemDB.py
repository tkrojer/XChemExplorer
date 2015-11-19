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

            'CompoundName',
            'CrystalTag',

            'DataCollectionBeamline',
            'DataCollectionDate',
            'DataCollectionOutcome',
            'DataCollectionRun',
            'DataCollectionComment',



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

            'PANDDAstuff'

            'RefinementRcryst',
            'RefinementRfree',
            'RefinementLigandCC',
            'RefinementRmsdBonds',
            'RefinementRmsdAngles',
            'RefinementOutcome',

            'AssayIC50'

        ]

        # 1. check data source if all columns exist
        # 2. create missing columns

        existing_columns=[]

        if self.data_source_type=='sqlite':
#            sample_dict={}
            self.connect=sqlite3.connect(self.data_source_file)
            self.connect.row_factory = sqlite3.Row
            self.cursor = self.connect.cursor()
            self.cursor.execute("SELECT * FROM mainTable")
            for column in self.cursor.description:
                existing_columns.append(column[0])
            for column in self.column_list:
                if column not in existing_columns:
                    self.cursor.execute("alter table mainTable add column '%s' 'TEXT'" % column)
                    self.connect.commit()

        if self.data_source_type=='csv':
            for n,line in enumerate(open(self.data_source_file)):
                if n==0:
                    for item in line.split(','):
                        existing_columns.append(item.replace('\r','').replace('\n',''))
                break
            new_columns=''
            for column in self.column_list:
                if column not in existing_columns:
                    new_columns+=column+','
            csv_content=''
            for n,line in enumerate(open(self.data_source_file)):
                if n==0:
                    csv_content+=line.replace('\r','').replace('\n','')+','+new_columns+'\n'
                else:
                    csv_content+=line
            f=open(self.data_source_file,'w')
            f.write(csv_content)
            f.close()



#rows = cur.fetchall()

#for row in rows:
#    print row["CompoundSMILES"]
#    print row
#    sample_dict[row["ID"]]=[row["CompoundCode"],row["CompoundSMILES"]]
#print sample_dict
#for key in sample_dict:
#    print key,sample_dict[key]


    def load_samples_from_data_source(self):

        data=[]

        columns_of_interest=[   'CrystalName',
                                'CompoundName',
                                'CompoundCode',
                                'CompoundSMILES',
                                'CrystalTag'        ]

        # first find the index of the columns of interest in the header
        if self.data_source_type=='csv':
            data=[]
            column_list=[]
            for n,line in enumerate(open(self.data_source_file)):
                if n==0:
                    for item in columns_of_interest:
                        for column,field in enumerate(line.split(',')):
                            if field==item:
                                column_list.append(column)
                else:
                    line_list=[]
                    for item in column_list:
                        found=False
                        for column,field in enumerate(line.split(',')):
                            if column==item:
                                line_list.append(field)
                                found=True
                        if not found:
                            line_list.append('')
                    data.append(tuple(line_list))
            return data
#
#                        if column in column_list:
#                    print len(line.split(','))

        if self.data_source_type=='sqlite':
            self.cursor.execute("SELECT "+str(columns_of_interest).translate(None,"'[]")+" FROM mainTable")
            data = self.cursor.fetchall()
            return data
#            for row in rows:
#                print row

#            for column in self.cursor.description:
#                print column


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

