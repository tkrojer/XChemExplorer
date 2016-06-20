import os,glob
import sys
import subprocess
import getpass
from datetime import datetime

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

import XChemDB


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

    dimple_jobs=0
    dimple_job_details=[]
    pandda_jobs=0
    pandda_job_details=[]
    refmac_jobs=0
    refmac_job_details=[]
    others_jobs=0
    others_job_details=[]

    # note: each job_details list contains a list with
    # [job_ID, status, run_time]
    out=subprocess.Popen(['qstat'],stdout=subprocess.PIPE)
    for n,line in enumerate(iter(out.stdout.readline,'')):
        print len(line.split()),line.split()
        if len(line.split()) >= 7:
            if line.split()[3] == getpass.getuser():
                start_date=''
                start_time=''
                start_date=line.split()[5]
                print start_date
                print len(start_date.split('/'))
                if len(start_date.split('/')) == 3:
                    month_start=start_date.split('/')[0]
                    day_start=start_date.split('/')[1]
                    year_start=start_date.split('/')[2]
                    
                start_time=line.split()[6]
                print start_time
                print len(start_time.split(':'))
                if len(start_time.split(':')) == 3:
                    hour_start=start_time.split(':')[0]
                    minute_start=start_time.split(':')[1]
                    second_start=start_time.split(':')[2]

                if start_time != '' and start_date != '':
                    start='%s-%s-%s %s:%s:%s' %(year_start,month_start,day_start,hour_start,minute_start,second_start)
                    print start
                    delta_time=datetime.now()-datetime.strptime(start,"%Y-%m-%d %H:%M:%S")
                    print delta_time
                    print int(delta_time.total_seconds() / 60)

#                    print  datetime.strptime(start,"%Y-%m-%d %H:%M:%S")-datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#                    now =  datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#                    delta_time = now - then
#                    print delta_time
                if 'dimple' in line.split()[2]:
                    dimple_jobs+=1
                    dimple_job_details.append(line.split()[0])
                elif 'pandda' in line.split()[2]:
                    pandda_jobs+=1
                    pandda_job_details.append([line.split()[0]],line.split(4)  )
                elif 'refmac' in line.split()[2]:
                    refmac_jobs+=1
                    refmac_job_details.append(line.split()[0])
                else:
                    others_jobs+=1
                    others_job_details.append(line.split()[0])

    out_dict['dimple']=[dimple_jobs,dimple_job_details]
    out_dict['pandda']=[pandda_jobs,pandda_job_details]
    out_dict['refmac']=[refmac_jobs,refmac_job_details]
    out_dict['others']=[others_jobs,others_job_details]

    return out_dict

def get_datasource_summary(db_file):
    print 'db file', db_file
    db=XChemDB.data_source(db_file)
    out_dict={}

    out_dict['no_samples']=len(db.execute_statement("select CrystalName from mainTable where CrystalName is not NULL;"))
    out_dict['no_smiles_for_samples']=len(db.execute_statement("select compoundSMILES from mainTable where compoundSMILES is not (NULL or '')"))
    out_dict['no_data_collection_success']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'success';"))

    out_dict['no_data_collection_success']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - centring failed';"))
    out_dict['no_data_collection_fail-no-diffraction']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - no diffraction';"))
    out_dict['no_data_collection_fail-processing']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - processing';"))
    out_dict['no_data_collection_fail-loop-empty']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - loop empty';"))
    out_dict['no_data_collection_fail-loop-broken']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - loop broken';"))
    out_dict['no_data_collection_fail-low-resolution']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - low resolution';"))
    out_dict['no_data_collection_fail-no-X-rays']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - no X-rays';"))
    out_dict['no_data_collection_fail-unknown']=len(db.execute_statement("select DataCollectionOutcome from mainTable where DataCollectionOutcome is 'Failed - unknown';"))

    out_dict['no_cif_files']=len(db.execute_statement("select RefinementCIF from mainTable where RefinementCIF is not (Null or '');"))

    out_dict['no_analysis-pending']=len(db.execute_statement("select RefinementOutcome from mainTable where RefinementOutcome is '1 - Analysis Pending';"))
    out_dict['no_pandda-models']=len(db.execute_statement("select RefinementOutcome from mainTable where RefinementOutcome is '2 - PANDDA model';"))
    out_dict['no_in-refinement']=len(db.execute_statement("select RefinementOutcome from mainTable where RefinementOutcome is '3 - In Refinement';"))
    out_dict['no_comp-chem-ready']=len(db.execute_statement("select RefinementOutcome from mainTable where RefinementOutcome is '4 - ComChem ready';"))
    out_dict['no_deposition-ready']=len(db.execute_statement("select RefinementOutcome from mainTable where RefinementOutcome is '5 - Deposition ready';"))

    return out_dict

def remove_all_refmac_jobs_from_cluster_and_reinstate_last_stable_state():
    print 'hallo'

