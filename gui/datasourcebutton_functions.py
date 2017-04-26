from import_modules import *

def update_header_and_data_from_datasource(self):
    self.update_log.insert('getting information for all samples from data source...')
    self.db=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
    self.update_log.insert('creating missing columns in data source')
    self.db.create_missing_columns()
    self.update_log.insert('load header and data from data source')
    self.header,self.data=self.db.load_samples_from_data_source()
    self.update_log.insert('get all samples in data source')
    all_samples_in_db=self.db.execute_statement("select CrystalName from mainTable where CrystalName is not '';")

    self.xtal_db_dict={}
    sampleID_column=0
    for n,entry in enumerate(self.header):
        if entry=='CrystalName':
            sampleID_column=n
            break
    for line in self.data:
        if str(line[sampleID_column]) != '':
            db_dict={}
            for n,entry in enumerate(line):
                if n != sampleID_column:
                    db_dict[str(self.header[n])]=str(entry)
            self.xtal_db_dict[str(line[sampleID_column])]=db_dict

    print '==> XCE: found '+str(len(self.xtal_db_dict))+' samples'

def datasource_menu_reload_samples(self):
    self.update_log.insert('reading samples from data source: '+os.path.join(self.database_directory,self.data_source_file))
    self.update_status_bar('reading samples from data source: '+os.path.join(self.database_directory,self.data_source_file))
    self.update_header_and_data_from_datasource()
    self.update_all_tables()