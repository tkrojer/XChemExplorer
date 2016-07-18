import os,glob
import sys
import subprocess
import getpass
from datetime import datetime

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

import XChemDB
import XChemLog

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
    print 'db file', db_file
    db=XChemDB.data_source(db_file)
    out_dict={}

    out_dict['no_samples']=len(db.execute_statement("select CrystalName from mainTable where CrystalName is not NULL;"))
    out_dict['no_samples_failed_to_mount']=len(db.execute_statement("select HarvestStatus from mainTable where HarvestStatus is 'fail';"))

    out_dict['no_smiles_for_samples']=len(db.execute_statement("select compoundSMILES from mainTable where compoundSMILES is not (NULL or '')"))

    out_dict['no_data_collection_success']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'success';"))
    out_dict['no_data_collection_centring_fail']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - centring failed';"))
    out_dict['no_data_collection_no-diffraction']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - no diffraction';"))
    out_dict['no_data_collection_processing_fail']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - processing';"))
    out_dict['no_data_collection_loop-empty']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - loop empty';"))
    out_dict['no_data_collection_loop-broken']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - loop broken';"))
    out_dict['no_data_collection_low-resolution']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - low resolution';"))
    out_dict['no_data_collection_no-X-rays']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - no X-rays';"))
    out_dict['no_data_collection_unknown']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - unknown';"))

    out_dict['no_cif_files']=len(db.execute_statement("select RefinementCIF from mainTable where RefinementCIF is not (Null or '');"))

    out_dict['no_analysis-pending']=len(db.execute_statement("select RefinementOutcome from mainTable where RefinementOutcome is '1 - Analysis Pending';"))
    out_dict['no_pandda-models']=len(db.execute_statement("select RefinementOutcome from mainTable where RefinementOutcome is '2 - PANDDA model';"))
    out_dict['no_in-refinement']=len(db.execute_statement("select RefinementOutcome from mainTable where RefinementOutcome is '3 - In Refinement';"))
    out_dict['no_comp-chem-ready']=len(db.execute_statement("select RefinementOutcome from mainTable where RefinementOutcome is '4 - ComChem ready';"))
    out_dict['no_deposition-ready']=len(db.execute_statement("select RefinementOutcome from mainTable where RefinementOutcome is '5 - Deposition ready';"))
    print out_dict
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


                # first check if folders and files exist
                # since user might do this before data are actually copied over

                if os.path.isdir(os.path.join(initial_model_directory,sample,'autoprocessing',visit+'-'+run+autoproc)):
                    db_dict['DataProcessingAutoAssigned']='False'
                    os.chdir(os.path.join(initial_model_directory,sample))
                    # first remove old links
                    os.system('/bin/rm '+sample+'.mtz')
                    os.system('/bin/rm '+sample+'.log')
                    # make new links
                    Logfile.insert('setting symlink: '+os.path.join(path_to_logfile,log_filename)+' -> '+sample+'.log')
                    os.symlink(os.path.join(path_to_logfile,log_filename),sample+'.log')
                    Logfile.insert('setting symlink: '+os.path.join(path_to_mtzfile,mtz_filename)+' -> '+sample+'.mtz')
                    os.symlink(os.path.join(path_to_mtzfile,mtz_filename),sample+'.mtz')

                    # update data source
                    data_source=XChemDB.data_source(data_source_file)
                    data_source.update_insert_data_source(sample,db_dict)

                else:
                    Logfile.insert('please copy data to PROJECT DIRECTORY first!')
