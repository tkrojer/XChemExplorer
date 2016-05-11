import os, sys, glob
from datetime import datetime
from PyQt4 import QtGui, QtCore
from XChemUtils import mtztools
import XChemDB
import csv

class run_pandda_export(QtCore.QThread):

    def __init__(self,panddas_directory,datasource):
        QtCore.QThread.__init__(self)
        self.panddas_directory=panddas_directory
        self.datasource=datasource
        self.db=XChemDB.data_source(self.datasource)
        self.db.create_missing_columns()

    def run(self):
        print 'fuirguir'
        self.import_samples_into_datasouce()

#    def get_db_dict(self):
#
#        sample_dict={}
#
#        with open(os.path.join(self.panddas_directory,'pandda_inspect.csv'),'rb') as csv_import:
#            csv_dict = csv.DictReader(csv_import)
#            for i,line in enumerate(csv_dict):
#                sampleID=line['dtag']
#                site_index=line['site_idx']
#                if int(site_index) > 15:        # currently data source does not support more than 15 sites
#                    continue
#                if sampleID not in sample_dict:
#                    sample_dict[sampleID]={}
#
#                db_dict=sample_dict[sampleID]
#
#                db_dict['PANDDA_site_'+str(site_index)+'_index']=line['site_idx']
#                db_dict['PANDDA_site_'+str(site_index)+'_comment']=line['Comment']
#                db_dict['PANDDA_site_'+str(site_index)+'_confidence']=line['Ligand Confidence']
#                db_dict['PANDDA_site_'+str(site_index)+'_ligand_placed']=line['Ligand Placed']
#                db_dict['PANDDA_site_'+str(site_index)+'_viewed']=line['Viewed']
#                db_dict['PANDDA_site_'+str(site_index)+'_interesting']=line['Interesting']
#                db_dict['PANDDA_site_'+str(site_index)+'_z_peak']=line['z_peak']
#
#                sample_dict[sampleID]=db_dict

    def import_samples_into_datasouce(self):

        db_dict=self.get_db_dict()
        print 'hallo'
#        for xtal in db_dict:
#            self.db.update_data_source(xtal,db_dict[xtal])


    def get_db_dict(self):


#        db_dict['PANDDA_site_A_index']
#        db_dict['PANDDA_site_A_name']=line['site_idx']
#        db_dict['PANDDA_site_A_comment']=line['Comment']
#        db_dict['PANDDA_site_A_confidence']=line['Ligand Confidence']
#        db_dict['PANDDA_site_A_ligand_placed']=line['Ligand Placed']
#        db_dict['PANDDA_site_A_viewed']=line['Viewed']
#        db_dict['PANDDA_site_A_interesting']=line['Interesting']
#        db_dict['PANDDA_site_A_z_peak']=line['z_peak']


        sample_dict={}

        with open(os.path.join(self.panddas_directory,'analyses','pandda_inspect.csv'),'rb') as csv_import:
            csv_dict = csv.DictReader(csv_import)
            for i,line in enumerate(csv_dict):
                sampleID=line['dtag']
                site_index=line['site_idx']
                if int(site_index) > 15:        # currently data source does not support more than 15 sites
                    continue
                if sampleID not in sample_dict:
                    sample_dict[sampleID]=[]

                db_list=sample_dict[sampleID]

                db_list.append([    line['site_idx'],
                                    line['Comment'],
                                    line['Ligand Confidence'],
                                    line['Ligand Placed'],
                                    line['Viewed'],
                                    line['Interesting'],
                                    line['z_peak']              ])


        for xtal in sample_dict:
            # sort by site index
            sample_dict[xtal].sort()
            for entry in sample_dict[xtal]:
                print xtal,entry

#pandda.export pandda_dir=/dls/labxchem/data/2015/lb13385-2/processing/test/pandda export_dir=./testx export_ligands=False generate_occupancy_groupings=True

        return sample_dict

class run_pandda_analyse(QtCore.QThread):

    def __init__(self,pandda_params):
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

    def run(self):
        if os.path.isfile(os.path.join(self.panddas_directory,'pandda.running')):
            return None
        else:
            if os.getenv('SHELL') == '/bin/tcsh' or os.getenv('SHELL') == '/bin/csh':
                source_file=os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-csh')
            elif os.getenv('SHELL') == '/bin/bash':
                source_file=os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-sh')
            else:
                source_file=''

            os.chdir(self.panddas_directory)
            Cmds = (
                '#!'+os.getenv('SHELL')+'\n'
                '\n'
                'source '+source_file+'\n'
                '\n'
                'cd '+self.panddas_directory+'\n'
                '\n'
                'pandda.analyse '
                ' data_dirs="'+self.data_directory+'"'
                ' out_dir='+self.panddas_directory+
                ' min_build_datasets='+self.min_build_datasets+
                ' maps.ampl_label=FWT maps.phas_label=PHWT'
                ' cpus='+self.nproc+
                ' pdb_style='+self.pdb_style+
                ' mtz_style='+self.mtz_style+'\n'
                )
            print Cmds

            f = open('pandda.sh','w')
            f.write(Cmds)
            f.close()
            if self.submit_mode=='local machine':
                print '==> running PANDDA on local machine'
                os.system('chmod +x pandda.sh')
                os.system('./pandda.sh &')
            else:
                print '==> running PANDDA on cluster, using qsub...'
                os.system('qsub pandda.sh')

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