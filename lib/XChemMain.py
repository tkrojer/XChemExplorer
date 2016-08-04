# last edited: 04/08/2016 - 17:00

import os,glob
import sys
import subprocess
import getpass
from datetime import datetime
from PyQt4 import QtGui, QtCore

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

import XChemDB
import XChemLog


def space_group_list():
    space_group_list =   [      'P1',
                                'P2','P21','C121','P1211','P121',
                                'P222','P2122','P2212','P2221',
                                'P21212','P21221','P22121','P212121',
                                'C222','C2221','F222', 'I222', 'I212121',
                                'P4','P41','P42','P43','I4','I41',
                                'P422','P4212','P4122','P41212','P4222',
                                'P42212','P4322','P43212','I422','I4122' ,
                                'P3','P31','P32',
                                'P312','P321','P3112','P3121','P3212','P3221',
                                'P6','P61','P65','P62','P64','P63',
                                'P622','P6122','P6522','P6222','P6422','P6322',
                                'H3','H32',
                                'P23','F23','I23','P213','I213',
                                'P432','P4232','F432','F4132','I432','P4332','P4132','I4132'        ]
    return space_group_list


def get_target_and_visit_list(beamline_directory):
#    target_list=['*']      # always give the option to read in all targets
    target_list=['=== SELECT TARGET ===']      # always give the option to read in all targets
    visit_list=[]
    # the beamline directory could be a the real directory or
    # a directory where the visits are linked into
    if len(beamline_directory.split('/')) and \
        beamline_directory.split('/')[1]=='dls' and beamline_directory.split('/')[3]=='data' \
        and not 'labxchem' in beamline_directory:
        visit_list.append(beamline_directory)
    elif os.path.islink(beamline_directory):
        visit_list.append(os.path.realpath(beamline_directory))
    else:
        for dir in glob.glob(beamline_directory+'/*'):
            visit_list.append(os.path.realpath(dir))

    for visit in visit_list:
        for target in glob.glob(os.path.join(visit,'processed','*')):
            if target[target.rfind('/')+1:] not in ['results','README-log','edna-latest.html']:
                if target[target.rfind('/')+1:] not in target_list:
                    target_list.append(target[target.rfind('/')+1:])
    return target_list,visit_list

def get_target_and_visit_list_for_Pietro(beamline_directory):
#    target_list=['*']      # always give the option to read in all targets
    target_list=['=== SELECT TARGET ===']      # always give the option to read in all targets
    visit_list=[]
    # the beamline directory could be a the real directory or
    # a directory where the visits are linked into
    for stuff in glob.glob(os.path.join(beamline_directory,'*')):
        visit_list.append(stuff)

    for visit in visit_list:
        for target in glob.glob(os.path.join(visit,'processed','*')):
            if target[target.rfind('/')+1:] not in ['results','README-log','edna-latest.html']:
                if target[target.rfind('/')+1:] not in target_list:
                    target_list.append(target[target.rfind('/')+1:])
    return target_list,visit_list


def get_dewar_configuration(beamline_directory):
    dewar_sample_configuration_dict={}
    # need to do it like this for now until we got a way to query ispyb
    for xml in glob.glob(os.path.join(beamline_directory,'xml','exptTableParams-*.xml')):
        prefix=''
        container_reference=''
        sample_location=''
        for line in open(xml):
            if 'prefix' in line:
                prefix=line[line.find('>')+1:line.rfind('<')]
            if 'container_reference' in line:
                container_reference=line[line.find('>')+1:line.rfind('<')]
            if 'sample_location' in line:
                sample_location=line[line.find('>')+1:line.rfind('<')]
        dewar_sample_configuration_dict[str(container_reference)+'-'+str(sample_location)]=prefix
    return dewar_sample_configuration_dict

def get_jobs_running_on_cluster():
    out_dict={}

    dimple_jobs=[]
    acedrg_jobs=[]
    pandda_jobs=[]
    refmac_jobs=[]
    others_jobs=[]

    # note: each job_details list contains a list with
    # [job_ID, status, run_time]
    out=subprocess.Popen(['qstat'],stdout=subprocess.PIPE)
    for n,line in enumerate(iter(out.stdout.readline,'')):
        if len(line.split()) >= 7:
            if line.split()[3] == getpass.getuser():

                job_id = line.split()[0]
                job_name = line.split()[2]
                job_status = line.split()[4]

                ##########################################################
                # determine run time of each job in minutes
                start_date=''
                start_time=''
                run_time_minutes=''
                start_date=line.split()[5]
                if len(start_date.split('/')) == 3:
                    month_start=start_date.split('/')[0]
                    day_start=start_date.split('/')[1]
                    year_start=start_date.split('/')[2]
                    
                start_time=line.split()[6]
                if len(start_time.split(':')) == 3:
                    hour_start=start_time.split(':')[0]
                    minute_start=start_time.split(':')[1]
                    second_start=start_time.split(':')[2]

                if start_time != '' and start_date != '':
                    start='%s-%s-%s %s:%s:%s' %(year_start,month_start,day_start,hour_start,minute_start,second_start)
                    run_time=datetime.now()-datetime.strptime(start,"%Y-%m-%d %H:%M:%S")
                    run_time_seconds=int(run_time.total_seconds())
                    run_time_minutes=int(run_time.total_seconds() / 60)
                    run_time_hours=int(run_time.total_seconds() / 3600)

                ##########################################################
                # determine run time of each job in minutes
                if 'dimple' in job_name:
                    dimple_jobs.append([job_id,job_status,run_time_minutes])
                elif 'acedrg' in job_name:
                    acedrg_jobs.append([job_id,job_status,run_time_minutes])
                elif 'pandda' in job_name:
                    pandda_jobs.append([job_id,job_status,run_time_minutes])
                elif 'refmac' in job_name:
                    refmac_jobs.append([job_id,job_status,run_time_minutes])
                else:
                    others_jobs.append([job_id,job_status,run_time_minutes])

    out_dict['dimple']=dimple_jobs
    out_dict['pandda']=pandda_jobs
    out_dict['refmac']=refmac_jobs
    out_dict['others']=others_jobs

    return out_dict

def display_queue_status_in_terminal(in_dict):
    # latest_run=max(tmp,key=lambda x: x[1])[0]
    max_dimple_runtime=max(in_dict['dimple'],key=lambda  x:x[2])[2]
    print '----------------------------------------------------------------------------'
    print '| Task                   | Nr. Jobs               | Max. Runtime (minutes) |'
    print '----------------------------------------------------------------------------'
    print '{0:24} {1:24} {2:24} {3:1}'.format('| DIMPLE', '| %s' %len(in_dict['dimple']),'| %s' &max_dimple_runtime,'|')
    print '----------------------------------------------------------------------------'

def get_datasource_summary(db_file):
    db=XChemDB.data_source(db_file)

    out_dict={}

    out_dict['nr_samples']=len(db.execute_statement("select CrystalName from mainTable where CrystalName is not NULL;"))
    out_dict['nr_samples_failed_to_mount']=len(db.execute_statement("select HarvestStatus from mainTable where HarvestStatus is 'fail';"))

    out_dict['nr_smiles_for_samples']=len(db.execute_statement("select compoundSMILES from mainTable where compoundSMILES is not (NULL or '')"))

    out_dict['nr_data_collection_success']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'success';"))
    out_dict['nr_data_collection_centring_fail']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - centring failed';"))
    out_dict['nr_data_collection_no-diffraction']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - no diffraction';"))
    out_dict['nr_data_collection_processing_fail']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - processing';"))
    out_dict['nr_data_collection_loop-empty']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - loop empty';"))
    out_dict['nr_data_collection_loop-broken']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - loop broken';"))
    out_dict['nr_data_collection_low-resolution']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - low resolution';"))
    out_dict['nr_data_collection_no-X-rays']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - no X-rays';"))
    out_dict['nr_data_collection_unknown']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - unknown';"))

    out_dict['nr_data_collection_failed']=  out_dict['nr_data_collection_centring_fail'] + \
                                            out_dict['nr_data_collection_no-diffraction'] + \
                                            out_dict['nr_data_collection_processing_fail'] + \
                                            out_dict['nr_data_collection_loop-empty'] + \
                                            out_dict['nr_data_collection_loop-broken'] + \
                                            out_dict['nr_data_collection_low-resolution'] + \
                                            out_dict['nr_data_collection_no-X-rays'] + \
                                            out_dict['nr_data_collection_unknown']


    out_dict['nr_data_collection_pending']=out_dict['nr_samples'] - out_dict['nr_data_collection_success'] - \
                                                                    out_dict['nr_data_collection_centring_fail'] - \
                                                                    out_dict['nr_data_collection_no-diffraction'] - \
                                                                    out_dict['nr_data_collection_processing_fail'] - \
                                                                    out_dict['nr_data_collection_loop-empty'] - \
                                                                    out_dict['nr_data_collection_loop-broken'] - \
                                                                    out_dict['nr_data_collection_low-resolution'] - \
                                                                    out_dict['nr_data_collection_no-X-rays'] - \
                                                                    out_dict['nr_data_collection_unknown']

    out_dict['nr_initial_maps_available']=len(db.execute_statement("select DimplePathToMTZ from mainTable where DimplePathToMTZ is not '';"))
    out_dict['nr_initial_maps_fail']=len(db.execute_statement("select DataProcessingDimpleSuccessful from mainTable where DataProcessingDimpleSuccessful = 'False';"))
    out_dict['nr_initial_maps_pending']=out_dict['nr_data_collection_success']-out_dict['nr_initial_maps_available']-out_dict['nr_initial_maps_fail']

    out_dict['nr_pandda_hits']=len(db.execute_statement("select DimplePANDDAhit from mainTable where DimplePANDDAhit = 'True';"))
    out_dict['nr_pandda_reject']=len(db.execute_statement("select DimplePANDDAreject from mainTable where DimplePANDDAreject = 'True';"))
    out_dict['nr_pandda_processed']=len(db.execute_statement("select DimplePANDDAwasRun from mainTable where DimplePANDDAwasRun = 'True';"))-out_dict['nr_pandda_hits']-out_dict['nr_pandda_reject']
    out_dict['nr_pandda_pending']=out_dict['nr_initial_maps_available']-out_dict['nr_pandda_hits']-out_dict['nr_pandda_reject']-out_dict['nr_pandda_processed']

    out_dict['nr_cif_files']=len(db.execute_statement("select RefinementCIF from mainTable where RefinementCIF is not (Null or '');"))

    out_dict['nr_analysis-pending']=len(db.execute_statement("select RefinementOutcome from mainTable where RefinementOutcome is '1 - Analysis Pending';"))
    out_dict['nr_pandda-models']=len(db.execute_statement("select RefinementOutcome from mainTable where RefinementOutcome is '2 - PANDDA model';"))
    out_dict['nr_in-refinement']=len(db.execute_statement("select RefinementOutcome from mainTable where RefinementOutcome is '3 - In Refinement';"))
    out_dict['nr_comp-chem-ready']=len(db.execute_statement("select RefinementOutcome from mainTable where RefinementOutcome is '4 - ComChem ready';"))
    out_dict['nr_deposition-ready']=len(db.execute_statement("select RefinementOutcome from mainTable where RefinementOutcome is '5 - Deposition ready';"))

    return out_dict


def remove_all_refmac_jobs_from_cluster_and_reinstate_last_stable_state():
    print 'hallo'


def change_links_to_selected_data_collection_outcome(sample,data_collection_dict,data_collection_column_three_dict,dataset_outcome_dict,initial_model_directory,data_source_file,xce_logfile):
    Logfile=XChemLog.updateLog(xce_logfile)
    # find out which row was selected in respective data collection table
    selected_processing_result='n/a'
    indexes=data_collection_column_three_dict[sample][0].selectionModel().selectedRows()
    if indexes != []:       # i.e. logfile exists
        for index in sorted(indexes):
            selected_processing_result=index.row()

    for n,entry in enumerate(data_collection_dict[sample]):
        if entry[0]=='logfile':
            if entry[7]==selected_processing_result:
                visit=entry[1]
                run=entry[2]
                autoproc=entry[4]
                db_dict=entry[6]
                outcome=dataset_outcome_dict[sample]
                path_to_logfile=db_dict['DataProcessingPathToLogfile']
                path_to_mtzfile=db_dict['DataProcessingPathToMTZfile']
                mtz_filename=db_dict['DataProcessingMTZfileName']
                log_filename=db_dict['DataProcessingLOGfileName']
                relative_path_to_mtzfile='.'+path_to_mtzfile.replace(initial_model_directory,'')
                relative_path_to_logfile='.'+path_to_logfile.replace(initial_model_directory,'')


                # first check if folders and files exist
                # since user might do this before data are actually copied over

                if os.path.isdir(os.path.join(initial_model_directory,sample,'autoprocessing',visit+'-'+run+autoproc)):
                    db_dict['DataProcessingAutoAssigned']='False'
                    Logfile.insert('changing directory to: '+os.path.join(initial_model_directory,sample))
                    os.chdir(os.path.join(initial_model_directory,sample))
                    # first remove old links
                    os.system('/bin/rm '+sample+'.mtz 2> /dev/null')
                    os.system('/bin/rm '+sample+'.log 2> /dev/null')
                    # make new links
#                    Logfile.insert('setting symlink: '+os.path.join(path_to_logfile,log_filename)+' -> '+sample+'.log')
#                    os.symlink(os.path.join(path_to_logfile,log_filename),sample+'.log')
#                    Logfile.insert('setting symlink: '+os.path.join(path_to_mtzfile,mtz_filename)+' -> '+sample+'.mtz')
#                    os.symlink(os.path.join(path_to_mtzfile,mtz_filename),sample+'.mtz')
                    Logfile.insert('setting relative symlink: '+os.path.join(relative_path_to_logfile,log_filename)+' -> '+sample+'.log')
                    os.symlink(os.path.join(relative_path_to_logfile,log_filename),sample+'.log')
                    Logfile.insert('setting relative symlink: '+os.path.join(relative_path_to_mtzfile,mtz_filename)+' -> '+sample+'.mtz')
                    os.symlink(os.path.join(relative_path_to_mtzfile,mtz_filename),sample+'.mtz')

                    # update data source
                    data_source=XChemDB.data_source(data_source_file)
                    data_source.update_insert_data_source(sample,db_dict)

                else:
                    Logfile.insert('please copy data to PROJECT DIRECTORY first!')


#def find_diffraction_image_directory(diffraction_data_directory):
#    data_dict={}
#    diffraction_image_extension = ['.img','.cbf','.mccd','.mar2560','.mar2300']
#    os.chdir(diffraction_data_directory)
#    for xtal in glob.glob('*'):
#        data_dict[xtal]=[]
#        for root,dirs,files in os.walk(xtal):
#            if 'screening' in root:
#                continue
#            for n,image_file in enumerate(glob.glob(os.path.join(root,'*'))):
#                file_extension=image_file[image_file.rfind('.'):]
#                if n > 20 and file_extension in diffraction_image_extension:
#                    data_dict[xtal].append(os.path.join(diffraction_data_directory,root))
#                    break
#        if data_dict[xtal]==[]:
#            del data_dict[xtal]
#    return data_dict

class find_diffraction_image_directory(QtCore.QThread):
    def __init__(self,diffraction_data_directory):
        QtCore.QThread.__init__(self)
        self.diffraction_data_directory=diffraction_data_directory
        self.data_dict={}
        self.diffraction_image_extension = ['.img','.cbf','.mccd','.mar2560','.mar2300']

    def run(self):
        os.chdir(self.diffraction_data_directory)
        if len(glob.glob(os.path.join(self.diffraction_data_directory,'*'))) != 0:
            progress_step=100/float(len(glob.glob(os.path.join(self.diffraction_data_directory,'*'))))
        else:
            progress_step=100
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'searching diffraction data directory')
        for xtal in glob.glob('*'):
            self.data_dict[xtal]=[]
            for root,dirs,files in os.walk(xtal):
                if 'screening' in root:
                    continue
                for n,image_file in enumerate(glob.glob(os.path.join(root,'*'))):
                    file_extension=image_file[image_file.rfind('.'):]
                    if n > 20 and file_extension in self.diffraction_image_extension:
                        self.data_dict[xtal].append(os.path.join(self.diffraction_data_directory,root))
                        break
            if self.data_dict[xtal]==[]:
                del self.data_dict[xtal]

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.emit(QtCore.SIGNAL('update_reprocess_datasets_table'), self.data_dict)
