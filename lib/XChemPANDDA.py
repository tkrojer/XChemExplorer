# last edited: 12/01/2017, 15:00

import os, sys, glob
from datetime import datetime
from PyQt4 import QtGui, QtCore
import math
from XChemUtils import mtztools
import XChemDB
import XChemRefine
import XChemUtils
import XChemLog
import csv

def get_names_of_current_clusters(xce_logfile,panddas_directory):
    Logfile=XChemLog.updateLog(xce_logfile)
    Logfile.insert('parsing %s/cluster_analysis' %panddas_directory)
    os.chdir('%s/cluster_analysis' %panddas_directory)
    cluster_dict={}
    for out_dir in sorted(glob.glob('*')):
        if os.path.isdir(out_dir):
            cluster_dict[out_dir]=[]
            found_first_pdb=False
            for folder in glob.glob(os.path.join(out_dir,'pdbs','*')):
                xtal=folder[folder.rfind('/')+1:]
                if not found_first_pdb:
                    if os.path.isfile(os.path.join(panddas_directory,'cluster_analysis',out_dir,'pdbs',xtal,xtal+'.pdb') ):
                        cluster_dict[out_dir].append(os.path.join(panddas_directory,'cluster_analysis',out_dir,'pdbs',xtal,xtal+'.pdb'))
                        found_first_pdb=True
                cluster_dict[out_dir].append(xtal)
#    for key in cluster_dict:
#        Logfile.insert('cluster %s:   %s datasets' %(str(key),str(len(cluster_dict[key])-1)))
    return cluster_dict



class run_pandda_export(QtCore.QThread):

    def __init__(self,panddas_directory,datasource,initial_model_directory,xce_logfile,update_datasource_only,which_models):
        QtCore.QThread.__init__(self)
        self.panddas_directory=panddas_directory
        self.datasource=datasource
        self.initial_model_directory=initial_model_directory
        self.db=XChemDB.data_source(self.datasource)
        self.db.create_missing_columns()
        self.db_list=self.db.get_empty_db_dict()
        self.external_software=XChemUtils.external_software(xce_logfile).check()
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.update_datasource_only=update_datasource_only
        self.which_models=which_models
        self.already_exported_models=[]

        self.RefmacParams={ 'HKLIN':            '',                 'HKLOUT': '',
                            'XYZIN':            '',                 'XYZOUT': '',
                            'LIBIN':            '',                 'LIBOUT': '',
                            'TLSIN':            '',                 'TLSOUT': '',
                            'TLSADD':           '',
                            'NCYCLES':          '10',
                            'MATRIX_WEIGHT':    'AUTO',
                            'BREF':             '    bref ISOT\n',
                            'TLS':              '',
                            'NCS':              '',
                            'TWIN':             ''    }

    def run(self):

        self.import_samples_into_datasouce()

        if not self.update_datasource_only:
            self.export_models()

#        self.import_samples_into_datasouce()

        if not self.update_datasource_only:
            self.refine_exported_models()


    def refine_exported_models(self):

        if self.which_models=='new':
            sample_list=self.db.execute_statement("select CrystalName,CompoundCode from mainTable where RefinementOutcome='2 - PANDDA model';")
        elif self.which_models=='all':
            sample_list=self.db.execute_statement("select CrystalName,CompoundCode from mainTable where RefinementOutcome='2 - PANDDA model' or RefinementOutcome='3 - In Refinement';")
        for item in sample_list:
            xtal=str(item[0])
            compoundID=str(item[1])
            if os.path.isfile(os.path.join(self.initial_model_directory,xtal,xtal+'.free.mtz')):
                if os.path.isfile(os.path.join(self.initial_model_directory,xtal,xtal+'-ensemble-model.pdb')):
                    self.Logfile.insert('running inital refinement on PANDDA model of'+xtal)
                    Refine=XChemRefine.Refine(self.initial_model_directory,xtal,compoundID,self.datasource)
                    Serial=Refine.GetSerial()
                    os.mkdir(os.path.join(self.initial_model_directory,xtal,'Refine_'+str(Serial)))
                    os.chdir(os.path.join(self.initial_model_directory,xtal,'Refine_'+str(Serial)))
                    os.symlink(os.path.join(self.initial_model_directory,xtal,xtal+'-ensemble-model.pdb'),'in.pdb')
                    Refine.RunRefmac(Serial,self.RefmacParams,self.external_software,self.xce_logfile)


    def import_samples_into_datasouce(self):

        # first make a note of all the datasets which were used in pandda directory
        os.chdir(os.path.join(self.panddas_directory,'processed_datasets'))
        for xtal in glob.glob('*'):
#            self.Logfile.insert("update mainTable set DimplePANDDAwasRun = 'True',DimplePANDDAreject = 'False',DimplePANDDApath='%s' where CrystalName is '%s'" %(self.panddas_directory,xtal))
            self.db.execute_statement("update mainTable set DimplePANDDAwasRun = 'True',DimplePANDDAreject = 'False',DimplePANDDApath='%s' where CrystalName is '%s'" %(self.panddas_directory,xtal))
        # do the same as before, but look for rejected datasets
        os.chdir(os.path.join(self.panddas_directory,'rejected_datasets'))
        for xtal in glob.glob('*'):
#            self.Logfile.insert("update mainTable set DimplePANDDAwasRun = 'True',DimplePANDDAreject = 'True',DimplePANDDApath='%s',DimplePANDDAhit = 'False' where CrystalName is '%s'" %(self.panddas_directory,xtal))
            self.db.execute_statement("update mainTable set DimplePANDDAwasRun = 'True',DimplePANDDAreject = 'True',DimplePANDDApath='%s',DimplePANDDAhit = 'False' where CrystalName is '%s'" %(self.panddas_directory,xtal))

        site_list = []
        pandda_hit_list=[]

        with open(os.path.join(self.panddas_directory,'analyses','pandda_inspect_sites.csv'),'rb') as csv_import:
            csv_dict = csv.DictReader(csv_import)
            for i,line in enumerate(csv_dict):
                site_index=line['site_idx']
                name=line['Name'].replace("'","")
                comment=line['Comment']
                site_list.append([site_index,name,comment])


        progress_step=1
        for i,line in enumerate(open(os.path.join(self.panddas_directory,'analyses','pandda_inspect_events.csv'))):
            n_lines=i
        if n_lines != 0:
            progress_step=100/float(n_lines)
        else:
            progress_step=0
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        with open(os.path.join(self.panddas_directory,'analyses','pandda_inspect_events.csv'),'rb') as csv_import:
            csv_dict = csv.DictReader(csv_import)

            for i,line in enumerate(csv_dict):
                db_dict={}
                sampleID=line['dtag']
                if sampleID not in pandda_hit_list:
                    pandda_hit_list.append(sampleID)
                site_index=line['site_idx']
                event_index=line['event_idx']

                for entry in site_list:
                    if entry[0]==site_index:
                        site_name=entry[1]
                        site_comment=entry[2]
                        break

                # check if EVENT map exists in project directory
                event_map='event_map'
                for file in glob.glob(os.path.join(self.initial_model_directory,sampleID,'*ccp4')):
                    filename=file[file.rfind('/')+1:]
                    if filename.startswith(sampleID+'-event_'+event_index) and filename.endswith('map.native.ccp4'):
                        event_map=file
                        break

                # initial pandda model and mtz file
                pandda_model='pandda_model'
                for file in glob.glob(os.path.join(self.initial_model_directory,sampleID,'*pdb')):
                    filename=file[file.rfind('/')+1:]
                    if filename.endswith('-ensemble-model.pdb'):
                        pandda_model=file
                        if sampleID not in self.already_exported_models:
                            self.already_exported_models.append(sampleID)
                        break
                inital_mtz='initial_mtz'
                for file in glob.glob(os.path.join(self.initial_model_directory,sampleID,'*mtz')):
                    filename=file[file.rfind('/')+1:]
                    if filename.endswith('pandda-input.mtz'):
                        inital_mtz=file
                        break

                db_dict['CrystalName']                  =   sampleID
                db_dict['PANDDApath']                   =   self.panddas_directory
                db_dict['PANDDA_site_index']            =   site_index
                db_dict['PANDDA_site_name']             =   site_name
                db_dict['PANDDA_site_comment']          =   site_comment
                db_dict['PANDDA_site_event_index']      =   event_index
                db_dict['PANDDA_site_event_comment']    =   line['Comment'].replace("'","")
                db_dict['PANDDA_site_confidence']       =   line['Ligand Confidence']
                db_dict['PANDDA_site_ligand_placed']    =   line['Ligand Placed']
                db_dict['PANDDA_site_viewed']           =   line['Viewed']
                db_dict['PANDDA_site_interesting']      =   line['Interesting']
                db_dict['PANDDA_site_z_peak']           =   line['z_peak']
                db_dict['PANDDA_site_x']                =   line['x']
                db_dict['PANDDA_site_y']                =   line['y']
                db_dict['PANDDA_site_z']                =   line['z']
                db_dict['PANDDA_site_ligand_id']        =   ''
                db_dict['PANDDA_site_event_map']        =   ''
                db_dict['PANDDA_site_initial_model']    =   pandda_model
                db_dict['PANDDA_site_initial_mtz']      =   inital_mtz
                db_dict['PANDDA_site_spider_plot']      =   ''
#                db_dict['RefinementOutcome']            =   '2 - PANDDA model'

                # find apo structures which were used
                # XXX missing XXX

                self.db.update_insert_panddaTable(sampleID,db_dict)
                # this is necessary, otherwise RefinementOutcome will be reset for samples that are actually already in refinement
                self.db.execute_statement("update panddaTable set RefinementOutcome = '2 - PANDDA model' where CrystalName is '%s' and RefinementOutcome is null" %sampleID)
                self.db.execute_statement("update mainTable set RefinementOutcome = '2 - PANDDA model' where CrystalName is '%s' and (RefinementOutcome is null or RefinementOutcome is '1 - Analysis Pending')" %sampleID)
                self.db.execute_statement("update mainTable set DimplePANDDAhit = 'True' where CrystalName is '%s'" %sampleID)
                progress += progress_step
                self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        # finally find all samples which do not have a pandda hit
        os.chdir(os.path.join(self.panddas_directory,'processed_datasets'))
        for xtal in glob.glob('*'):
            if xtal not in pandda_hit_list:
                self.db.execute_statement("update mainTable set DimplePANDDAhit = 'False' where CrystalName is '%s'" %xtal)


    def export_models(self):

        self.Logfile.insert('finding out which PanDDA models need to be exported')

        # first find which samples are in interesting datasets and have a model
        # and determine the timestamp
        fileModelsDict={}
        queryModels=''
        for model in glob.glob(os.path.join(self.panddas_directory,'processed_datasets','*','modelled_structures','*-pandda-model.pdb')):
            sample=model[model.rfind('/')+1:].replace('-pandda-model.pdb','')
            timestamp=datetime.fromtimestamp(os.path.getmtime(model)).strftime('%Y-%m-%d %H:%M:%S')
            self.Logfile.insert(sample+'-pandda-model.pdb was created on '+str(timestamp))
            queryModels+="'"+sample+"',"
            fileModelsDict[sample]=timestamp

        # now get these models from the database and compare the datestamps
        # Note: only get the models that underwent some form of refinement,
        #       because only if the model was updated in pandda.inspect will it be exported and refined
        dbModelsDict={}
        if queryModels != '':
            print "select CrystalName,DatePanDDAModelCreated from mainTable where CrystalName in (%s) and (RefinementOutcome like '3%' or RefinementOutcome like '4%' or RefinementOutcome like '5%')" %queryModels[:-1]
            dbEntries=self.db.execute_statement("select CrystalName,DatePanDDAModelCreated from mainTable where CrystalName in (%s) and (RefinementOutcome like '3%' or RefinementOutcome like '4%' or RefinementOutcome like '5%')") %queryModels[:-1]
            for item in dbEntries:
                xtal=str(item[0])
                timestamp=str(item[1])
                dbModelsDict[xtal]=timestamp
                self.Logfile.insert('PanDDA model for '+xtal+' is in database and was created on '+str(timestamp))

        # compare timestamps and only export the ones where the timestamp of the file is newer than the one in the DB
        samples_to_export={}
        self.Logfile.insert('checking which PanDDA models were newly created or updated')
        for sample in fileModelsDict:
            if sample in dbModelsDict:
                try:
                    difference=(datetime.strptime(fileModelsDict[sample],'%Y-%m-%d %H:%M:%S') - datetime.strptime(dbModelsDict[sample],'%Y-%m-%d %H:%M:%S')  )
                    if difference.seconds != 0:
                        samples_to_export[sample]=fileModelsDict[sample]
                except ValueError:
                    # this will be raised if timestamp is not properly formatted;
                    # which will usually be the case when respective field in database is blank
                    # these are hopefully legacy cases which are from before this extensive check was introduced (13/01/2017)
                    advice = (  'The pandda model of '+xtal+' was changed, but it was already refined! '
                                'This is most likely because this was done with an older version of XCE. '
                                'If you really want to export and refine this model, you need to open the database '
                                'with DBbroweser (http://sqlitebrowser.org/); then change the RefinementOutcome field '
                                'of the respective sample to "2 - PANDDA model", save the database and repeat the export prodedure.'
                                )
                    self.Logfile.insert(advice)
            else:
                samples_to_export[sample]=fileModelsDict[sample]

        # update the DB:
        # set timestamp to current timestamp of file and set RefinementOutcome to '2-pandda...'

        if samples_to_export != {}:
            select_dir_string=''
            for sample in samples_to_export:
                db_dict={}
                db_dict['RefinementOutcome']='2 - PANDDA model'
                db_dict['DatePanDDAModelCreated']=samples_to_export[sample]
                select_dir_string+="select_dir=%s " %sample
                self.Logfile.insert('updating database for '+sample+' setting time model was created to '+db_dict['DatePanDDAModelCreated']+' and RefinementOutcome to '+db_dict['RefinementOutcome'])
                self.db.update_data_source(sample,db_dict)

            Cmds = (
                'source '+os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-sh')+'\n'
                '\n'
                '/dls/science/groups/i04-1/software/pandda-install/ccp4-pandda/bin/pandda.export'
                ' pandda_dir=%s' %self.panddas_directory+
                ' export_dir=%s' %self.initial_model_directory+
                ' %s' %select_dir_string+
                ' export_ligands=False'
                ' generate_occupancy_groupings=True\n'
                )

            self.Logfile.insert('running pandda.export with the following settings:\n'+Cmds)
            os.system(Cmds)

#        Cmds = (
#                '#!'+os.getenv('SHELL')+'\n'
#                'unset PYTHONPATH\n'
#                'module load ccp4\n'
#                '$CCP4/bin/pandda.export'
#                ' pandda_dir=%s' %self.panddas_directory+
#                ' export_dir=%s' %self.initial_model_directory+
#                ' export_ligands=False'
#                ' generate_occupancy_groupings=True\n'
#                )
#        os.system(Cmds)
#        os.system('pandda.export pandda_dir=%s export_dir=%s export_ligands=False generate_occupancy_groupings=True' %(self.panddas_directory,self.initial_model_directory))
#        self.Logfile.insert('ran pandda.export with the following command:\n'+Cmds)
#        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'running pandda.export: check terminal for details')





class run_pandda_analyse(QtCore.QThread):

    def __init__(self,pandda_params,xce_logfile,dataset_list,datasource):
        QtCore.QThread.__init__(self)
        self.data_directory=pandda_params['data_dir']
        self.panddas_directory=pandda_params['out_dir']
        self.submit_mode=pandda_params['submit_mode']
#        if self.submit_mode == 'local machine':
#            self.nproc=pandda_params['nproc']
#        else:
#            self.nproc='7'
        self.nproc=pandda_params['nproc']
        self.min_build_datasets=pandda_params['min_build_datasets']
        self.pdb_style=pandda_params['pdb_style']
        self.mtz_style=pandda_params['mtz_style']
        self.sort_event=pandda_params['sort_event']
        self.number_of_datasets=pandda_params['N_datasets']
        self.max_new_datasets=pandda_params['max_new_datasets']
        self.grid_spacing=pandda_params['grid_spacing']
        self.filter_pdb=pandda_params['filter_pdb']
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.dataset_list=dataset_list
        self.datasource=datasource
        self.db=XChemDB.data_source(datasource)

    def run(self):

        # how to run pandda.analyse on large datasets
        #
        # 1) Run the normal pandda command, with the new setting, e.g.
        # pandda.analyse data_dirs=... max_new_datasets=500
        # This will do the analysis on the first 500 datasets and build the statistical maps - just as normal.
        #
        # 2) Run pandda with the same command:
        # pandda.analyse data_dirs=... max_new_datasets=500
        # This will add 500 new datasets, and process them using the existing statistical maps
        # (this will be quicker than the original analysis). It will then merge the results of the two analyses.
        #
        # 3) Repeat 2) until you don't add any "new" datasets. Then you can build the models as normal.

        crystalString=''
        for n,dataset in enumerate(self.dataset_list):
            if n > 0:       # first entry is reference file!
                crystalString+="'"+dataset+"',"
        print ("update mainTable set PANDDAStatus = 'started' where CrystalName in (%s)" %crystalString[:-1])
        self.db.execute_statement("update mainTable set PANDDAStatus = 'started' where CrystalName in (%s)" %crystalString[:-1])

        number_of_cyles=int(self.number_of_datasets)/int(self.max_new_datasets)
        if int(self.number_of_datasets) % int(self.max_new_datasets) != 0:  # modulo gives remainder after integer division
            number_of_cyles+=1

        if os.path.isfile(os.path.join(self.panddas_directory,'pandda.running')):
            self.Logfile.insert('it looks as if a pandda.analyse job is currently running in: '+self.panddas_directory)
            msg = ( 'there are three possibilities:\n'
                    '1.) choose another PANDDA directory\n'
                    '2.) - check if the job is really running either on the cluster (qstat) or on your local machine\n'
                    '    - if so, be patient and wait until the job has finished\n'
                    '3.) same as 2., but instead of waiting, kill the job and remove at least the pandda.running file\n'
                    '   (or all the contents in the directory if you want to start from scratch)\n'   )
            self.Logfile.insert(msg)
            return None
        else:
            if os.getenv('SHELL') == '/bin/tcsh' or os.getenv('SHELL') == '/bin/csh':
                source_file=os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-csh')
            elif os.getenv('SHELL') == '/bin/bash':
                source_file=os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-sh')
            else:
                source_file=''

            if os.path.isfile(self.filter_pdb):
                filter_pdb=' filter.pdb='+self.filter_pdb
            else:
                filter_pdb=''

            os.chdir(self.panddas_directory)
            Cmds = (
                '#!'+os.getenv('SHELL')+'\n'
                '\n'
                'source '+source_file+'\n'
                '\n'
                'cd '+self.panddas_directory+'\n'
                '\n'
                '$CCP4/bin/ccp4-python $XChemExplorer_DIR/helpers/update_pandda_status_flag.py %s %s %s\n' %(self.datasource,crystalString[:-1],'running') +
                '\n'
                )

            for i in range(number_of_cyles):
                Cmds += (
                    '\n'
                    'pandda.analyse '
                    ' data_dirs="'+self.data_directory+'"'
                    ' out_dir='+self.panddas_directory+
                    ' min_build_datasets='+self.min_build_datasets+
                    ' maps.ampl_label=FWT maps.phas_label=PHWT'
                    ' max_new_datasets='+self.max_new_datasets+
                    ' grid_spacing='+self.grid_spacing+
                    ' cpus='+self.nproc+
                    ' events.order_by='+self.sort_event+
                    filter_pdb+
                    ' pdb_style='+self.pdb_style+
                    ' mtz_style='+self.mtz_style+'\n'
                    '\n'
                    )

            Cmds += '$CCP4/bin/ccp4-python $XChemExplorer_DIR/helpers/update_pandda_status_flag.py %s %s %s\n' %(self.datasource,crystalString[:-1],'finished')

            self.Logfile.insert('running pandda.analyse with the following command:\n'+Cmds)

            f = open('pandda.sh','w')
            f.write(Cmds)
            f.close()
            if self.submit_mode=='local machine':
                self.Logfile.insert('running PANDDA on local machine')
                os.system('chmod +x pandda.sh')
                os.system('./pandda.sh &')
            else:
                self.Logfile.insert('running PANDDA on cluster, using qsub...')
                os.system('qsub -P labxchem pandda.sh')

        self.emit(QtCore.SIGNAL('datasource_menu_reload_samples'))



class giant_cluster_datasets(QtCore.QThread):

    def __init__(self,initial_model_directory,pandda_params,xce_logfile,datasource,run_pandda_analyse):
        QtCore.QThread.__init__(self)
        self.panddas_directory=pandda_params['out_dir']
        self.pdb_style=pandda_params['pdb_style']
        self.mtz_style=pandda_params['mtz_style']
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.initial_model_directory=initial_model_directory
        self.db=XChemDB.data_source(datasource)
        self.run_pandda_analyse=run_pandda_analyse

    def run(self):

        self.emit(QtCore.SIGNAL('update_progress_bar'), 0)

        # 1.) prepare output directory
        os.chdir(self.panddas_directory)
        if os.path.isdir('cluster_analysis'):
            self.Logfile.insert('removing old cluster_analysis directory in %s' %self.panddas_directory)
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'removing old cluster_analysis directory in %s' %self.panddas_directory)
            os.system('/bin/rm -fr cluster_analysis 2> /dev/null')
        self.Logfile.insert('creating cluster_analysis directory in %s' %self.panddas_directory)
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'creating cluster_analysis directory in %s' %self.panddas_directory)
        os.mkdir('cluster_analysis')
        self.emit(QtCore.SIGNAL('update_progress_bar'), 10)

        # 2.) go through project directory and make sure that all pdb files really exist
        # broken links derail the giant.cluster_mtzs_and_pdbs script
        self.Logfile.insert('cleaning up broken links of %s and %s in %s' %(self.pdb_style,self.mtz_style,self.initial_model_directory))
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'cleaning up broken links of %s and %s in %s' %(self.pdb_style,self.mtz_style,self.initial_model_directory))
        os.chdir(self.initial_model_directory)
        for xtal in glob.glob('*'):
            if not os.path.isfile(os.path.join(xtal,self.pdb_style)):
                self.Logfile.insert('missing %s and %s for %s' %(self.pdb_style,self.mtz_style,xtal))
                os.system('/bin/rm %s/%s 2> /dev/null' %(xtal,self.pdb_style))
                os.system('/bin/rm %s/%s 2> /dev/null' %(xtal,self.mtz_style))
        self.emit(QtCore.SIGNAL('update_progress_bar'), 20)

        # 3.) giant.cluster_mtzs_and_pdbs
        self.Logfile.insert("running giant.cluster_mtzs_and_pdbs %s/*/%s pdb_regex='%s/(.*)/%s' out_dir='%s/cluster_analysis'" %(self.initial_model_directory,self.pdb_style,self.initial_model_directory,self.pdb_style,self.panddas_directory))
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'running giant.cluster_mtzs_and_pdbs')

        if os.getenv('SHELL') == '/bin/tcsh' or os.getenv('SHELL') == '/bin/csh':
            source_file=os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-csh')
        elif os.getenv('SHELL') == '/bin/bash':
            source_file=os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-sh')
        else:
            source_file=''

        Cmds = (
                '#!'+os.getenv('SHELL')+'\n'
                'unset PYTHONPATH\n'
                'source '+source_file+'\n'
                "giant.cluster_mtzs_and_pdbs %s/*/%s pdb_regex='%s/(.*)/%s' out_dir='%s/cluster_analysis'" %(self.initial_model_directory,self.pdb_style,self.initial_model_directory,self.pdb_style,self.panddas_directory)
            )

#        os.system("giant.cluster_mtzs_and_pdbs %s/*/%s pdb_regex='%s/(.*)/%s' out_dir='%s/cluster_analysis'" %(self.initial_model_directory,self.pdb_style,self.initial_model_directory,self.pdb_style,self.panddas_directory))
        os.system(Cmds)
        self.emit(QtCore.SIGNAL('update_progress_bar'), 80)

        # 4.) analyse output
        self.Logfile.insert('parsing %s/cluster_analysis' %self.panddas_directory)
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'parsing %s/cluster_analysis' %self.panddas_directory)
        os.chdir('%s/cluster_analysis' %self.panddas_directory)
        cluster_dict={}
        for out_dir in sorted(glob.glob('*')):
            if os.path.isdir(out_dir):
                cluster_dict[out_dir]=[]
                for folder in glob.glob(os.path.join(out_dir,'pdbs','*')):
                    xtal=folder[folder.rfind('/')+1:]
                    cluster_dict[out_dir].append(xtal)
        self.emit(QtCore.SIGNAL('update_progress_bar'), 90)

        # 5.) update datasource
        self.Logfile.insert('updating datasource with results from giant.cluster_mtzs_and_pdbs')
        if cluster_dict != {}:
            for key in cluster_dict:
                for xtal in cluster_dict[key]:
                    db_dict={}
                    db_dict['CrystalFormName']=key
                    self.db.update_data_source(xtal,db_dict)
        self.emit(QtCore.SIGNAL('update_progress_bar'), 100)
        self.Logfile.insert('finished giant.cluster_mtzs_and_pdbs')
        self.emit(QtCore.SIGNAL('datasource_menu_reload_samples'))
        if self.run_pandda_analyse:
            self.emit(QtCore.SIGNAL('run_pandda_analyse'))

class check_if_pandda_can_run:

    # reasons why pandda cannot be run
    # - there is currently a job running in the pandda directory
    # - min datasets available is too low
    # - required input paramters are not complete
    # - map amplitude and phase labels don't exist

    def __init__(self,pandda_params,xce_logfile,datasource):
        self.data_directory=pandda_params['data_dir']
        self.panddas_directory=pandda_params['out_dir']
        self.min_build_datasets=pandda_params['min_build_datasets']
        self.pdb_style=pandda_params['pdb_style']
        self.mtz_style=pandda_params['mtz_style']
        self.input_dir_structure=pandda_params['pandda_dir_structure']
        self.problem_found=False
        self.error_code=-1
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.db=XChemDB.data_source(datasource)

    def number_of_available_datasets(self):
        counter=0
        for file in glob.glob(os.path.join(self.input_dir_structure,self.pdb_style)):
            if os.path.isfile(file):
                counter+=1
        self.Logfile.insert('pandda.analyse: found %s useable datasets' %counter)
        return counter

    def get_first_dataset_in_project_directory(self):
        first_dataset=''
        for file in glob.glob(os.path.join(self.input_dir_structure,self.pdb_style)):
            if os.path.isfile(file):
                first_dataset=file
                break
        return first_dataset

    def compare_number_of_atoms_in_reference_vs_all_datasets(self,refData,dataset_list):
        mismatched_datasets=[]
        pdbtools=XChemUtils.pdbtools(refData)
        refPDB=refData[refData.rfind('/')+1:]
        refPDBlist=pdbtools.get_init_pdb_as_list()
        n_atom_ref=len(refPDBlist)
        for n_datasets,dataset in enumerate(dataset_list):
            if os.path.isfile(os.path.join(self.data_directory.replace('*',''),dataset,self.pdb_style)):
                n_atom=len(pdbtools.get_pdb_as_list(os.path.join(self.data_directory.replace('*',''),dataset,self.pdb_style)))
                if n_atom_ref == n_atom:
                    self.Logfile.insert('%s: atoms in PDB file (%s): %s; atoms in Reference file: %s ===> OK' %(dataset,self.pdb_style,str(n_atom),str(n_atom_ref)))
                if n_atom_ref != n_atom:
                    self.Logfile.insert('%s: atoms in PDB file (%s): %s; atoms in Reference file: %s ===> ERROR' %(dataset,self.pdb_style,str(n_atom),str(n_atom_ref)))
                    mismatched_datasets.append(dataset)
        return n_datasets,mismatched_datasets

    def get_datasets_which_fit_to_reference_file(self,ref,reference_directory,cluster_dict,allowed_unitcell_difference_percent):
        refStructure=XChemUtils.pdbtools(os.path.join(reference_directory,ref+'.pdb'))
        symmRef=refStructure.get_spg_number_from_pdb()
        ucVolRef=refStructure.calc_unitcell_volume_from_pdb()
        cluster_dict[ref]=[]
        cluster_dict[ref].append(os.path.join(reference_directory,ref+'.pdb'))
        for dataset in glob.glob(os.path.join(self.data_directory,self.pdb_style)):
            datasetStructure=XChemUtils.pdbtools(dataset)
            symmDataset=datasetStructure.get_spg_number_from_pdb()
            ucVolDataset=datasetStructure.calc_unitcell_volume_from_pdb()
            if symmDataset == symmRef:
                try:
                    difference=math.fabs(1-(float(ucVolRef)/float(ucVolDataset)))*100
                    if difference < allowed_unitcell_difference_percent:
                        sampleID=dataset.replace('/'+self.pdb_style,'')[dataset.replace('/'+self.pdb_style,'').rfind('/')+1:]
                        cluster_dict[ref].append(sampleID)
                except ZeroDivisionError:
                    continue
        return cluster_dict

    def remove_dimple_files(self,dataset_list):
        for n_datasets,dataset in enumerate(dataset_list):
            db_dict={}
            if os.path.isfile(os.path.join(self.data_directory.replace('*',''),dataset,self.pdb_style)):
                os.system('/bin/rm '+os.path.join(self.data_directory.replace('*',''),dataset,self.pdb_style))
                self.Logfile.insert('%s: removing %s' %(dataset,self.pdb_style))
                db_dict['DimplePathToPDB']=''
                db_dict['DimpleRcryst']=''
                db_dict['DimpleRfree']=''
                db_dict['DimpleResolutionHigh']=''
                db_dict['DimpleStatus']='pending'
            if os.path.isfile(os.path.join(self.data_directory.replace('*',''),dataset,self.mtz_style)):
                os.system('/bin/rm '+os.path.join(self.data_directory.replace('*',''),dataset,self.mtz_style))
                self.Logfile.insert('%s: removing %s' %(dataset,self.mtz_style))
                db_dict['DimplePathToMTZ']=''
            if db_dict != {}:
                self.db.update_data_source(dataset,db_dict)


    def analyse_pdb_style(self):
        pdb_found=False
        for file in glob.glob(os.path.join(self.data_directory,self.pdb_style)):
            if os.path.isfile(file):
                pdb_found=True
                break
        if not pdb_found:
            self.error_code=1
            message=self.warning_messages()
        return message

    def analyse_mtz_style(self):
        mtz_found=False
        for file in glob.glob(os.path.join(self.data_directory,self.mtz_style)):
            if os.path.isfile(file):
                mtz_found=True
                break
        if not mtz_found:
            self.error_code=2
            message=self.warning_messages()
        return message

    def analyse_min_build_dataset(self):
        counter=0
        for file in glob.glob(os.path.join(self.data_directory,self.mtz_style)):
            if os.path.isfile(file):
                counter+=1
        if counter <= self.min_build_datasets:
            self.error_code=3
            message=self.warning_messages()
        return message

#    def analyse_amplitude_and_phase_labels(self):


#    def analyse_all_input_parameter(self):
#        print 'hallo'

    def warning_messages(self):
        message=''
        if self.error_code==1:
            message='PDB file does not exist'
        if self.error_code==2:
            message='MTZ file does not exist'
        if self.error_code==3:
            message='Not enough datasets available'

        return message


class convert_all_event_maps_in_database(QtCore.QThread):

    def __init__(self,initial_model_directory,xce_logfile,datasource):
        QtCore.QThread.__init__(self)
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.initial_model_directory=initial_model_directory
        self.datasource=datasource
        self.db=XChemDB.data_source(datasource)

    def run(self):
        sqlite = (
            'select'
            ' CrystalName,'
            ' PANDDA_site_event_map,'
            ' PANDDA_site_ligand_resname,'
            ' PANDDA_site_ligand_chain,'
            ' PANDDA_site_ligand_sequence_number,'
            ' PANDDA_site_ligand_altLoc '
            'from panddaTable '
            'where PANDDA_site_event_map not like "event%"'
        )

        query=self.db.execute_statement(sqlite)

        progress_step=1
        if len(query) != 0:
            progress_step=100/float(len(query))
        else:
            progress_step=1
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        for item in query:
            xtalID=str(item[0])
            event_map=str(item[1])
            resname=str(item[2])
            chainID=str(item[3])
            resseq=str(item[4])
            altLoc=str(item[5])

            if os.path.isfile(os.path.join(self.initial_model_directory,xtalID,'refine.pdb')):
                os.chdir(os.path.join(self.initial_model_directory,xtalID))
                self.Logfile.insert('extracting ligand (%s,%s,%s,%s) from refine.pdb' %(str(resname),str(chainID),str(resseq),str(altLoc)))
                XChemUtils.pdbtools('refine.pdb').save_specific_ligands_to_pdb(resname,chainID,resseq,altLoc)
                if os.path.isfile('ligand_%s_%s_%s_%s.pdb' %(str(resname),str(chainID),str(resseq),str(altLoc))):
                    ligand_pdb='ligand_%s_%s_%s_%s.pdb' %(str(resname),str(chainID),str(resseq),str(altLoc))
                    print os.path.join(self.initial_model_directory,xtalID,ligand_pdb)
                else:
                    self.Logfile.insert('could not extract ligand; trying next...')
                    continue
            else:
                self.Logfile.insert('directory: '+os.path.join(self.initial_model_directory,xtalID)+' -> cannot find refine.pdb; trying next')
                continue

            if os.path.isfile(os.path.join(self.initial_model_directory,xtalID,'refine.mtz')):
                resolution=XChemUtils.mtztools(os.path.join(self.initial_model_directory,xtalID,'refine.mtz')).get_high_resolution_from_mtz()
            else:
                self.Logfile.insert('directory: '+os.path.join(self.initial_model_directory,xtalID)+' -> cannot find refine.mtz; trying next')
                continue

            convert_event_map_to_SF(self.initial_model_directory,xtalID,event_map,ligand_pdb,self.xce_logfile,self.datasource,resolution).run()

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)



class convert_event_map_to_SF:
    def __init__(self,project_directory,xtalID,event_map,ligand_pdb,xce_logfile,db_file,resolution):
        self.Logfile=XChemLog.updateLog(xce_logfile)

        self.event_map=event_map
        if not os.path.isfile(self.event_map):
            self.Logfile.insert('cannot find Event map: '+self.event_map)
            self.Logfile.insert('cannot convert event_map to structure factors!')
            return None

        self.project_directory=project_directory
        self.xtalID=xtalID
        self.event_map=event_map
        self.ligand_pdb=ligand_pdb
        self.event=event_map[event_map.rfind('/')+1:].replace('.map','').replace('.ccp4','')
#        self.resolution=resolution
        self.db=XChemDB.data_source(db_file)
#        self.db_file=db_file
        self.resolution=resolution

    def run(self):
        os.chdir(os.path.join(self.project_directory,self.xtalID))

#        if not os.path.isfile(os.path.join(self.project_directory,self.xtalID,'2fofc.map')):
#            self.Logfile.insert('cannot find 2fofc.map in '+os.path.join(self.project_directory,self.xtalID))
#            self.Logfile.insert('--> need 2fofc.map to determine grid')
#            mtzin=''
#            if os.path.isfile(os.path.join(self.project_directory,self.xtalID,'refine.mtz')):
#                mtzin='refine.mtz'
#            elif os.path.isfile(os.path.join(self.project_directory,self.xtalID,'dimple.mtz')):
#                mtzin='dimple.mtz'
#            if mtzin != '':
#                self.calculate_electron_density_map(mtzin)
#            else:
#                self.Logfile.insert('cannot find refine.mtz or dimple.mtz in '+os.path.join(self.project_directory,self.xtalID))
#                self.Logfile.insert('cannot calculate structure factors for '+self.event_map)
#                self.Logfile.insert('stopping')
#                return None
#
#
#        ElectronDensityMap=XChemUtils.maptools(os.path.join(self.project_directory,self.xtalID,'2fofc.map'))
#
#        self.gridElectronDensityMap=ElectronDensityMap.grid_sampling
#        self.Logfile.insert('using '+str(self.gridElectronDensityMap)+' as grid')
#
#        self.space_group_numberElectronDensityMap=ElectronDensityMap.space_group_number
#        self.Logfile.insert('using '+str(self.space_group_numberElectronDensityMap)+' as space group')
#
#        self.space_group=ElectronDensityMap.space_group()
#        self.Logfile.insert('using '+str(self.space_group)+' as space group')
#
#        self.unit_cell=ElectronDensityMap.cell_dimensions
#        self.Logfile.insert('using '+str(self.unit_cell)+' as cell dimensions')
#
#        if not os.path.isfile(self.ligand_pdb):
#            self.Logfile.insert('cannot find '+self.ligand_pdb)
#            self.Logfile.insert('stopping')
#            return None
#
#        # prepare input script
#        self.prepare_conversion_script()
#
#        # run script
#        self.run_conversion_script()

        # run phenix.map_to_structure_factors
        self.run_phenix_map_to_structure_factors()


        # check if output files exist
        if not os.path.isfile('%s.mtz' %self.event):
            self.Logfile.insert('cannot find %s.mtz' %self.event)
        else:
            self.Logfile.insert('conversion successful, %s.mtz exists' %self.event)
            # update datasource with event_map_mtz information
            self.update_database()

#        # remove all temporary files
#        self.clean_output_directory()


    def calculate_electron_density_map(self,mtzin):
        missing_columns=False
        column_dict=XChemUtils.mtztools(mtzin).get_all_columns_as_dict()
        if 'FWT' in column_dict['F'] and 'PHWT' in column_dict['PHS']:
            labin=' labin F1=FWT PHI=PHWT\n'
        elif '2FOFCWT' in column_dict['F'] and 'PH2FOFCWT' in column_dict['PHS']:
            labin=' labin F1=2FOFCWT PHI=PH2FOFCWT\n'
        else:
            missing_columns=True

        if not missing_columns:
            os.chdir(os.path.join(self.project_directory,self.xtalID))
            cmd = (
                'fft hklin '+mtzin+' mapout 2fofc.map << EOF\n'
                +labin+
                'EOF\n'
                    )
            self.Logfile.insert('calculating 2fofc map from '+mtzin)
            os.system(cmd)
        else:
            self.Logfile.insert('cannot calculate 2fofc.map; missing map coefficients')

    def prepare_conversion_script(self):

        os.chdir(os.path.join(self.project_directory,self.xtalID))

        # see also:
        # http://www.phaser.cimr.cam.ac.uk/index.php/Using_Electron_Density_as_a_Model

        if os.getcwd().startswith('/dls'):
            phenix_module='module_load_phenix\n'
        else:
            phenix_module=''

        cmd = (
            '#!'+os.getenv('SHELL')+'\n'
            '\n'
            +phenix_module+
            '\n'
            'pdbset XYZIN %s XYZOUT mask_ligand.pdb << eof\n' %self.ligand_pdb+
            ' SPACEGROUP %s\n' %self.space_group+
            ' CELL %s\n' %(' '.join(self.unit_cell))+
            ' END\n'
            'eof\n'
            '\n'
            'ncsmask XYZIN mask_ligand.pdb MSKOUT mask_ligand.msk << eof\n'
            ' GRID %s\n' %(' '.join(self.gridElectronDensityMap))+
#            ' AXIS Z    X    Y\n'
            ' RADIUS 10\n'
            ' PEAK 1\n'
            'eof\n'
            '\n'
            'mapmask MAPIN %s MAPOUT onecell_event_map.map << eof\n' %self.event_map+
            ' XYZLIM CELL\n'
#            ' AXIS Z    X    Y\n'
            'eof\n'
            '\n'
            'maprot MAPIN onecell_event_map.map MSKIN mask_ligand.msk WRKOUT masked_event_map.map << eof\n'
            ' MODE FROM\n'
            ' SYMMETRY WORK %s\n' %self.space_group_numberElectronDensityMap+
            ' AVERAGE\n'
            ' ROTATE EULER 0 0 0\n'
            ' TRANSLATE 0 0 0\n'
            'eof\n'
            '\n'
            'mapmask MAPIN masked_event_map.map MAPOUT masked_event_map_fullcell.map << eof\n'
            ' XYZLIM CELL\n'
#            ' AXIS Z    X    Y\n'
            ' PAD 0.0\n'
            'eof\n'
            '\n'
            'sfall HKLOUT %s.mtz MAPIN masked_event_map_fullcell.map << eof\n' %self.event+
            ' LABOUT FC=FC_event PHIC=PHIC_event\n'
            ' MODE SFCALC MAPIN\n'
            ' RESOLUTION %s\n' %self.resolution+
            ' END\n'
        )

        self.Logfile.insert('preparing script for conversion of Event map to SF')
        f = open('eventMap2sf.sh','w')
        f.write(cmd)
        f.close()
        os.system('chmod +x eventMap2sf.sh')


    def run_conversion_script(self):
        self.Logfile.insert('running conversion script...')
        os.system('./eventMap2sf.sh')

    def run_phenix_map_to_structure_factors(self):
        if float(self.resolution) < 1.21:   # program complains if resolution is 1.2 or higher
            self.resolution='1.21'
        self.Logfile.insert('running phenix.map_to_structure_factors %s d_min=%s output_file_name=%s.mtz' %(self.event_map,self.resolution,self.event))
        os.system('phenix.map_to_structure_factors %s d_min=%s output_file_name=%s.mtz' %(self.event_map,self.resolution,self.event))

    def update_database(self):
        sqlite = ( "update panddaTable set "
                   " PANDDA_site_event_map_mtz = '%s' " %os.path.join(self.project_directory,self.event+'.mtz')+
                   " where PANDDA_site_event_map is '%s' " %self.event_map
                    )
        self.db.execute_statement(sqlite)
        self.Logfile.insert('updating data source: '+sqlite)

    def clean_output_directory(self):
        os.system('/bin/rm mask_targetcell.pdb')
        os.system('/bin/rm mask_targetcell.msk')
        os.system('/bin/rm onecell.map')
        os.system('/bin/rm masked_targetcell.map')
        os.system('/bin/rm masked_fullcell.map')
        os.system('/bin/rm eventMap2sf.sh')
        os.system('/bin/rm '+self.ligand_pdb)

