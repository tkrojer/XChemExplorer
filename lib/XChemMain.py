import os,glob
import subprocess
import getpass

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
    dimple_job_ID=[]
    pandda_jobs=0
    pandda_job_ID=[]
    refmac_jobs=0
    refmac_job_ID=[]
    others_jobs=0
    others_job_ID=[]

    out=subprocess.Popen(['qstat'],stdout=subprocess.PIPE)
    for n,line in enumerate(iter(out.stdout.readline,'')):
        if len(line.split()) >= 4:
            if line.split()[3] == getpass.getuser():
                if 'dimple' in line.split()[2]:
                    dimple_jobs+=1
                    dimple_job_ID.append(line.split()[0])
                elif 'pandda' in line.split()[2]:
                    pandda_jobs+=1
                    pandda_job_ID.append(line.split()[0])
                elif 'refmac' in line.split()[2]:
                    refmac_jobs+=1
                    refmac_job_ID.append(line.split()[0])
                else:
                    others_jobs+=1
                    others_job_ID.append(line.split()[0])

    out_dict['dimple']=[dimple_jobs,dimple_job_ID]
    out_dict['pandda']=[pandda_jobs,pandda_job_ID]
    out_dict['refmac']=[refmac_jobs,refmac_job_ID]
    out_dict['others']=[others_jobs,others_job_ID]

    return out_dict