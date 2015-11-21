import sqlite3
import os,sys

sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
from XChemUtils import parse

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
            'DataProcessingPathToLogfile',

            'PANDDAstuff',

            'RefinementRcryst',
            'RefinementRfree',
            'RefinementLigandCC',
            'RefinementRmsdBonds',
            'RefinementRmsdAngles',
            'RefinementOutcome',
            'RefinementMTZfree',
            'RefinementCIF',
            'RefinementPDB_latest',
            'RefinementMTZ_latest',
            'RefinementComment',

            'AssayIC50'

        ]


        #   [column_name in DB, column_name shown in XCE]

        self.column_list=[
            # from Lab36
            ['ID',                                   'ID'],
            ['LabVisit',                             'LabVisit'],
            ['LibraryPlate',                         'LibraryPlate'],
            ['SourceWell',                           'SourceWell'],
            ['LibraryName',                          'LibraryName'],
            ['CompoundSMILES',                       'Smiles'],
            ['CompoundCode',                         'Compound ID'],
            ['CrystalPlate',                         'CrystalPlate'],
            ['CrystalWell',                          'CrystalWell'],
            ['EchoX',                                'EchoX'],
            ['EchoY',                                'EchoY'],
            ['DropVolume',                           'DropVolume'],
            ['ProteinName',                          'ProteinName'],
            ['BatchNumber',                          'BatchNumber'],
            ['CompoundStockConcentration',           'CompoundStockConcentration'],
            ['CompoundConcentration',                'CompoundConcentration'],
            ['SolventFraction',                      'SolventFraction'],
            ['SoakTransferVol',                      'SoakTransferVol'],
            ['SoakStatus',                           'SoakStatus'],
            ['SoakTimestamp',                        'SoakTimestamp'],
            ['CryoStockFraction',                    'CryoStockFraction'],
            ['CryoFraction',                         'CryoFraction'],
            ['CryoWell',                             'CryoWell'],
            ['CryoTransferVolume',                   'CryoTransferVolume'],
            ['CryoStatus',                           'CryoStatus'],
            ['CryoTimestamp',                        'CryoTimestamp'],
            ['SoakingTime',                          'SoakingTime'],
            ['HarvestStatus',                        'HarvestStatus'],
            ['CrystalName',                          'Sample ID'],
            ['Puck',                                 'Puck'],
            ['PuckPosition',                         'PuckPosition'],
            ['PinBarcode',                           'PinBarcode'],
            ['MountingResult',                       'MountingResult'],
            ['MountingArrivalTime',                  'MountingArrivalTime'],
            ['MountedTimestamp',                     'MountedTimestamp'],
            ['MountingTime',                         'MountingTime'],
            ['ispybStatus',                          'ispybStatus'],
            ['DataCollectionVisit',                  'DataCollectionVisit'],
            # from XChemExplorer
            ['CompoundName',                         'Compound Name'],
            ['CrystalTag',                           'Tag'],
            ['DataCollectionBeamline',               'DataCollectionBeamline'],
            ['DataCollectionDate',                   'DataCollectionDate'],
            ['DataCollectionOutcome',                'DataCollectionOutcome'],
            ['DataCollectionRun',                    'DataCollectionRun'],
            ['DataCollectionComment',                'DataCollectionComment'],
            ['DataProcessingProgram',                'DataProcessingProgram'],
            ['DataProcessingSpaceGroup',             'DataProcessingSpaceGroup'],
            ['DataProcessingUnitCell',               'DataProcessingUnitCell'],
            ['DataProcessingResolutionOverall',      'DataProcessingResolutionOverall'],
            ['DataProcessingResolutionLow',          'DataProcessingResolutionLow'],
            ['DataProcessingResolutionHigh',         'DataProcessingResolutionHigh'],
            ['DataProcessingRmergeOverall',          'DataProcessingRmergeOverall'],
            ['DataProcessingRmergeLow',              'DataProcessingRmergeLow'],
            ['DataProcessingRmergeHigh',             'DataProcessingRmergeHigh'],
            ['DataProcessingIsigOverall',            'DataProcessingIsigOverall'],
            ['DataProcessingIsigLow',                'DataProcessingIsigLow'],
            ['DataProcessingIsigHigh',               'DataProcessingIsigHigh'],
            ['DataProcessingCompletenessOverall',    'DataProcessingCompletenessOverall'],
            ['DataProcessingCompletenessLow',        'DataProcessingCompletenessLow'],
            ['DataProcessingCompletenessHigh',       'DataProcessingCompletenessHigh'],
            ['DataProcessingMultiplicityOverall',    'DataProcessingMultiplicityOverall'],
            ['DataProcessingMultiplicityLow',        'DataProcessingMultiplicityLow'],
            ['DataProcessingMultiplicityHigh',       'DataProcessingMultiplicityHigh'],
            ['DataProcessingPathToLogfile',          'DataProcessingPathToLogfile'],
            ['PANDDAstuff',                          'PANDDAstuff'],
            ['RefinementRcryst',                     'RefinementRcryst'],
            ['RefinementRfree',                      'RefinementRfree'],
            ['RefinementLigandCC',                   'RefinementLigandCC'],
            ['RefinementRmsdBonds',                  'RefinementRmsdBonds'],
            ['RefinementRmsdAngles',                 'RefinementRmsdAngles'],
            ['RefinementOutcome',                    'RefinementOutcome'],
            ['RefinementMTZfree',                    'RefinementMTZfree'],
            ['RefinementCIF',                        'RefinementCIF'],
            ['RefinementPDB_latest',                 'RefinementPDB_latest'],
            ['RefinementMTZ_latest',                 'RefinementMTZ_latest'],
            ['RefinementComment',                    'RefinementComment'],
            ['AssayIC50',                            'AssayIC50']
        ]


    def return_column_list(self):
        return self.column_list


    def create_empty_data_source_file(self):

        if self.data_source_type=='csv':
            if not os.path.isfile(self.data_source_file):
                csv_header=''
                for column in self.column_list:
                    csv_header+=str(column[0])+','
                csv_header+='\n'
                f=open(self.data_source_file,'w')
                f.write(csv_header)
                f.close()



    def load_samples_from_data_source(self):

        # 1. check data source if all columns exist
        # 2. append missing columns

        existing_columns=[]

        if self.data_source_type=='sqlite':
            self.connect=sqlite3.connect(self.data_source_file)
            self.connect.row_factory = sqlite3.Row
            self.cursor = self.connect.cursor()
            self.cursor.execute("SELECT * FROM mainTable")
            for column in self.cursor.description:
                existing_columns.append(column[0])
            for column in self.column_list:
                if column[0] not in existing_columns:
                    self.cursor.execute("alter table mainTable add column '%s' 'TEXT'" % column[0])
                    self.connect.commit()

        if self.data_source_type=='csv':
            for n,line in enumerate(open(self.data_source_file)):
                if n==0:
                    for item in line.split(','):
                        existing_columns.append(item.replace('\r','').replace('\n',''))
                break
            new_columns=''
            for column in self.column_list:
                if column[0] not in existing_columns:
                    new_columns+=column[0]+','
            csv_content=''
            for n,line in enumerate(open(self.data_source_file)):
                if n==0:
                    csv_content+=line.replace('\r','').replace('\n','')+','+new_columns+'\n'
                else:
                    csv_content+=line
            f=open(self.data_source_file,'w')
            f.write(csv_content)
            f.close()

    # from now on we read all of them and decide in the main GUI which ones to display
#        columns_of_interest=[   'CrystalName',
#                                'CompoundName',
#                                'CompoundCode',
#                                'CompoundSMILES',
#                                'CrystalTag'        ]

        header=[]
        data=[]
        if self.data_source_type=='csv':
            for n,line in enumerate(open(self.data_source_file)):
                if n==0:
                    for field in line.split(','):
                        header.append(field)
                else:
                    tmp=[]
                    for field in line.split(','):
                        if '\n' in field or '\r' in field:
                            continue
                        else:
                            tmp.append(field)
                    if not all(element=='' for element in tmp):     # e.g. Mac OSX numbers adds lots of emtpy lines
                        data.append(tuple(tmp))



        if self.data_source_type=='sqlite':
            self.cursor.execute("SELECT * FROM mainTable")
            for column in self.cursor.description:
                header.append(column[0])
            data = self.cursor.fetchall()
#            print header
#            for i in data:
#                print data

        return ([header,data])


        # first find the index of the columns of interest in the header
#        if self.data_source_type=='csv':
#            data=[]
#            column_list=[]
#            for n,line in enumerate(open(self.data_source_file)):
#                if n==0:
#                    for item in columns_of_interest:
#                        for column,field in enumerate(line.split(',')):
#                            if field==item:
#                                column_list.append(column)
#                else:
#                    line_list=[]
#                    for item in column_list:
#                        found=False
#                        for column,field in enumerate(line.split(',')):
#                            if column==item:
#                                line_list.append(field)
#                                found=True
#                        if not found:
#                            line_list.append('')
#                    data.append(tuple(line_list))
#            return data
#
#                        if column in column_list:
#                    print len(line.split(','))

#        if self.data_source_type=='sqlite':
#            self.cursor.execute("SELECT "+str(columns_of_interest).translate(None,"'[]")+" FROM mainTable")
#            data = self.cursor.fetchall()
#            return data
#            for row in rows:
#                print row

#            for column in self.cursor.description:
#                print column

    def save_autoprocessing_results_to_data_source(self,sample,outcome,logfile):
        if logfile != None:
            aimless_results=parse().GetAimlessLog(logfile)
            columns_to_update=  [
                ['DataCollectionBeamline',              'n/a'],
                ['DataCollectionDate',                  'n/a'],
                ['DataCollectionOutcome',               outcome],
                ['DataCollectionRun',                   aimless_results['Run']],
                ['DataCollectionComment',               'n/a'],
                ['DataProcessingProgram',               aimless_results['AutoProc']],
                ['DataProcessingSpaceGroup',            aimless_results['SpaceGroup']],
                ['DataProcessingUnitCell',              aimless_results['UnitCell']],
                ['DataProcessingResolutionOverall',     aimless_results['ResolutionLow']+'-'+aimless_results['ResolutionHigh']],
                ['DataProcessingResolutionLow',         aimless_results['ResolutionLow']],
                ['DataProcessingResolutionHigh',        aimless_results['ResolutionHigh']],
                ['DataProcessingRmergeOverall',         aimless_results['RmergeOverall']],
                ['DataProcessingRmergeLow',             aimless_results['RmergeLow']],
                ['DataProcessingRmergeHigh',            aimless_results['RmergeHigh']],
                ['DataProcessingIsigOverall',           aimless_results['IsigOverall']],
                ['DataProcessingIsigLow',               aimless_results['IsigLow']],
                ['DataProcessingIsigHigh',              aimless_results['IsigHigh']],
                ['DataProcessingCompletenessOverall',   aimless_results['CompletenessOverall']],
                ['DataProcessingCompletenessLow',       aimless_results['CompletenessLow']],
                ['DataProcessingCompletenessHigh',      aimless_results['CompletenessHigh']],
                ['DataProcessingMultiplicityOverall',   aimless_results['MultiplicityOverall']],
                ['DataProcessingMultiplicityLow',       aimless_results['MultiplicityLow']],
                ['DataProcessingMultiplicityHigh',      aimless_results['MultiplicityHigh']],
                ['DataProcessingPathToLogfile',         logfile]    ]
        else:
            columns_to_update=  [
                ['DataCollectionBeamline',              'n/a'],
                ['DataCollectionDate',                  'n/a'],
                ['DataCollectionOutcome',               outcome],
                ['DataCollectionRun',                   'n/a'],
                ['DataCollectionComment',               'n/a'],
                ['DataProcessingProgram',               'n/a'],
                ['DataProcessingSpaceGroup',            'n/a'],
                ['DataProcessingUnitCell',              'n/a'],
                ['DataProcessingResolutionOverall',     'n/a'],
                ['DataProcessingResolutionLow',         'n/a'],
                ['DataProcessingResolutionHigh',        'n/a'],
                ['DataProcessingRmergeOverall',         'n/a'],
                ['DataProcessingRmergeLow',             'n/a'],
                ['DataProcessingRmergeHigh',            'n/a'],
                ['DataProcessingIsigOverall',           'n/a'],
                ['DataProcessingIsigLow',               'n/a'],
                ['DataProcessingIsigHigh',              'n/a'],
                ['DataProcessingCompletenessOverall',   'n/a'],
                ['DataProcessingCompletenessLow',       'n/a'],
                ['DataProcessingCompletenessHigh',      'n/a'],
                ['DataProcessingMultiplicityOverall',   'n/a'],
                ['DataProcessingMultiplicityLow',       'n/a'],
                ['DataProcessingMultiplicityHigh',      'n/a'],
                ['DataProcessingPathToLogfile',         'n/a']    ]



        if self.data_source_type=='csv':
            sample_column=None
            column_list=[]
            for n,line in enumerate(open(self.data_source_file)):
                if n==0:
                    for column,item in enumerate(line.split(',')):
                        if item=='CrystalName':
                            sample_column=column
                    for item in columns_to_update:
                        for column,field in enumerate(line.split(',')):
                            if field==item[0]:
                                column_list.append(column)

           # find sample line
            row_to_change=None
            if sample_column != None:
                for row,line in enumerate(open(self.data_source_file)):
                    if len(line.split(',')) >= sample_column:
                        if line.split(',')[sample_column].replace(' ','')==sample:
                            row_to_change=row

            csv_out=''
            if row_to_change==None:
                for line in open(self.data_source_file):
                    csv_out+=line
                csv_list=(','*max(column_list)).split(',')
                csv_list[sample_column]=sample
                for n,item in enumerate(column_list):
                    csv_list[item]=columns_to_update[n][1]
                csv_out+=str(csv_list).translate(None,"[']")+'\n'
            else:
                for row,line in enumerate(open(self.data_source_file)):
                    if row == row_to_change:
                        if len(line.split(',')) < max(column_list):
                            line+=(','*int(max(column_list)-len(x.split(','))))
                        csv_list=line.split(',')
                        for n,item in enumerate(column_list):
                            csv_list[item]=columns_to_update[n][1]
                        csv_out+=str(csv_list).translate(None,"[']")+'\n'
                    else:
                        csv_out+=line
            f=open(self.data_source_file,'w')
            f.write(csv_out)
            f.close()


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

