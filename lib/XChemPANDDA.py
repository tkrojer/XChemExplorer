# last edited: 28/07/2016, 16:50

import os, sys, glob
from datetime import datetime
from PyQt4 import QtGui, QtCore
from XChemUtils import mtztools
import XChemDB
import XChemRefine
import XChemUtils
import XChemLog
import csv

class run_pandda_export(QtCore.QThread):

    def __init__(self,panddas_directory,datasource,initial_model_directory,xce_logfile,update_datasource_only):
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

        if not self.update_datasource_only:
            self.export_models()

        self.import_samples_into_datasouce()

        if not self.update_datasource_only:
            self.refine_exported_models()


    def refine_exported_models(self):

        sample_list=self.db.execute_statement("select CrystalName,CompoundCode from mainTable where RefinementOutcome='2 - PANDDA model';")
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





 #       progress_step=1
 #       if len(db_dict) != 0:
 #           progress_step=100/float(len(db_dict))
 #       else:
 #           progress_step=0
 #       progress=0
#
#        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
#
#        for xtal in db_dict:
##            print '==> XCE: updating panddaTable of data source with PANDDA site information for',xtal
#            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'updating data source with PANDDA site information for '+xtal)
#            self.db.update_insert_panddaTable(xtal,db_dict[xtal])
#            self.db.execute_statement("update mainTable set RefinementOutcome = '2 - PANDDA model' where CrystalName is '%s' and RefinementOutcome is null or RefinementOutcome is '1 - Analysis Pending'" %xtal)
#            os.chdir(os.path.join(self.initial_model_directory,xtal))
#            if os.path.isfile(xtal+'-ensemble-model.pdb'):
#                if os.path.isfile('refine.pdb'):
#                    os.system('/bin/rm refine.pdb')
#                os.symlink(xtal+'-ensemble-model.pdb','refine.pdb')
#            if os.path.isfile(xtal+'-pandda-input.mtz'):
#                if os.path.isfile('refine.mtz'):
#                    os.system('/bin/rm refine.mtz')
#                os.symlink(xtal+'-pandda-input.mtz','refine.mtz')
#            progress += progress_step
#            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)


    def import_samples_into_datasouce(self):

        # first make a note of all the datasets which were used in pandda directory
        os.chdir(os.path.join(self.panddas_directory,'processed_datasets'))
        for xtal in glob.glob('*'):
            self.Logfile.insert("update mainTable set DimplePANDDAwasRun = 'True',DimplePANDDAreject = 'False',DimplePANDDApath='%s' where CrystalName is '%s'" %(self.panddas_directory,xtal))
            self.db.execute_statement("update mainTable set DimplePANDDAwasRun = 'True',DimplePANDDAreject = 'False',DimplePANDDApath='%s' where CrystalName is '%s'" %(self.panddas_directory,xtal))
        # do the same as before, but look for rejected datasets
        os.chdir(os.path.join(self.panddas_directory,'rejected_datasets'))
        for xtal in glob.glob('*'):
            self.Logfile.insert("update mainTable set DimplePANDDAwasRun = 'True',DimplePANDDAreject = 'True',DimplePANDDApath='%s',DimplePANDDAhit = 'False' where CrystalName is '%s'" %(self.panddas_directory,xtal))
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
                db_dict['PANDDA_site_ligand_id']        =   'LIG'
                db_dict['PANDDA_site_event_map']        =   event_map
                db_dict['PANDDA_site_initial_model']    =   pandda_model
                db_dict['PANDDA_site_initial_mtz']      =   inital_mtz
                db_dict['PANDDA_site_spider_plot']      =   ''
#                db_dict['RefinementOutcome']            =   '2 - PANDDA model'

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
        Cmds = (
                'source '+os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-sh')+'\n'
                '\n'
                'pandda.export'
                ' pandda_dir=%s' %self.panddas_directory+
                ' export_dir=%s' %self.initial_model_directory+
                ' export_ligands=False'
                ' generate_occupancy_groupings=True\n'
                )
        self.Logfile.insert('running pandda.export with the following command:\n'+Cmds)
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'running pandda.export: check terminal for details')
        os.system(Cmds)




class run_pandda_analyse(QtCore.QThread):

    def __init__(self,pandda_params,xce_logfile):
        QtCore.QThread.__init__(self)
        self.data_directory=pandda_params['data_dir']
        self.panddas_directory=pandda_params['out_dir']
        self.submit_mode=pandda_params['submit_mode']
        if self.submit_mode == 'local machine':
            self.nproc=pandda_params['nproc']
        else:
            self.nproc='7'
        self.min_build_datasets=pandda_params['min_build_datasets']
        self.pdb_style=pandda_params['pdb_style']
        self.mtz_style=pandda_params['mtz_style']
        self.sort_event=pandda_params['sort_event']
        self.number_of_datasets=pandda_params['N_datasets']
        self.max_new_datasets=pandda_params['max_new_datasets']
        self.grid_spacing=pandda_params['grid_spacing']
        self.filter_pdb=pandda_params['filter_pdb']
        self.Logfile=XChemLog.updateLog(xce_logfile)

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
                os.system('qsub pandda.sh')

class giant_cluster_datasets(QtCore.QThread):

    def __init__(self,initial_model_directory,pandda_params,xce_logfile):
        QtCore.QThread.__init__(self)
        self.panddas_directory=pandda_params['out_dir']
        self.pdb_style=pandda_params['pdb_style']
        self.mtz_style=pandda_params['mtz_style']
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.initial_model_directory=initial_model_directory

    def run(self):

        # first go through project directory and make sure that all pdb files really exist
        # broken links derail the giant.cluster_mtzs_and_pdbs script

        self.Logfile.insert('cleaning up broken links of %s and %s in %s' %(self.pdb_style,self.mtz_style,self.initial_model_directory))
        os.chdir(self.initial_model_directory)
        for xtal in glob.glob('*'):
            if not os.path.isfile(os.path.join(xtal,self.pdb_style)):
                self.Logfile.insert('missing %s and %s for %s' %(self.pdb_style,self.mtz_style,xtal))
                os.system('/bin/rm %s/%s 2> /dev/null' %(xtal,self.pdb_style))
                os.system('/bin/rm %s/%s 2> /dev/null' %(xtal,self.mtz_style))

        self.Logfile.insert("running giant.cluster_mtzs_and_pdbs %s/*/%s pdb_regex='%s/(.*)/%s' out_dir='%s'" %(self.initial_model_directory,self.pdb_style,self.initial_model_directory,self.pdb_style,self.panddas_directory))
        os.system("giant.cluster_mtzs_and_pdbs %s/*/%s pdb_regex='%s/(.*)/%s' out_dir='%s'" %(self.initial_model_directory,self.pdb_style,self.initial_model_directory,self.pdb_style,self.panddas_directory))


class check_if_pandda_can_run:

    # reasons why pandda cannot be run
    # - there is currently a job running in the pandda directory
    # - min datasets available is too low
    # - required input paramters are not complete
    # - map amplitude and phase labels don't exist

    def __init__(self,pandda_params):
        self.data_directory=pandda_params['data_dir']
        self.panddas_directory=pandda_params['out_dir']
        self.min_build_datasets=pandda_params['min_build_datasets']
        self.pdb_style=pandda_params['pdb_style']
        self.mtz_style=pandda_params['mtz_style']

        self.problem_found=False
        self.error_code=-1

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