import os
import sqlite3
from functools import partial


class Proasis:
    def __init__(self):
        pass

    def proasis_menu(self, xce_object):

        xce_object.proasis_directory = '/dls/science/groups/proasis/'

        if os.path.isfile(os.path.join(xce_object.database_directory, xce_object.data_source_file)):
            conn = sqlite3.connect(os.path.join(xce_object.database_directory, xce_object.data_source_file))
            c = conn.cursor()

        # Project details or add project in menu
        counter = 0
        try:
            # get protein name from soakDB - this will be the proasis project name
            for row in c.execute('SELECT Protein FROM soakDB;'):
                counter += 1
                # if there is only one protein name in soakDB - all is good - happy days
                if counter == 1:
                    xce_object.proasis_name = str(row[0])
                # otherwise - give a warning
                # TODO: If this is actually ever encountered - deal with it. Should be fine.
                if counter > 1:
                    print('WARNING: More than one protein name found (proasis)')
                # If the project directory already exists in proasis dir, project should exist in proasis

                if os.path.isdir(os.path.join('/dls/science/groups/proasis/LabXChem/', xce_object.proasis_name)):
                    # show project name in menu (no action when clicked)
                    xce_object.proasis_project_label = str('Project Name: ' + xce_object.proasis_name)
                    xce_object.proasis_project_function = lambda:Proasis().blank()

                else:
                    xce_object.proasis_project_label = str('Create Project for ' + xce_object.proasis_name + '...')
                    xce_object.proasis_project_function = partial(Proasis().create_project, name=xce_object.proasis_name, xce_object=xce_object)
                    xce_object.needtoaddproj=True


        # should catch if project doesnt exitst
        except (AttributeError, UnboundLocalError):
            print('cannot find %s' % os.path.join(xce_object.database_directory, xce_object.data_source_file))

        # Lead details or add lead in menu
        counter = 0
        try:


            # check if there is a lead (from soakDB)
            for row in c.execute('SELECT proasisID FROM proasisLead'):
                counter += 1
            if counter == 1:
                # If so, display id of lead in menu, no action if clicked
                xce_object.proasis_lead_label = str('Lead ID: ' + str(row[0]))
                xce_object.proasis_lead_function = lambda:Proasis().blank()



        # If no lead or sites file, error message. No action on click
        except:

            # otherwise, if you can find the pandda_analyse_sites.csv file, allow lead to be added
            if hasattr(xce_object, 'needtoaddproj'):
                xce_object.proasis_lead_label = 'To add a lead: create project, and then reload XCE'
                xce_object.proasis_lead_function = lambda:Proasis().blank()
                xce_object.needtoaddleads = True
            elif os.path.isfile(os.path.join(xce_object.panddas_directory, 'analyses/pandda_analyse_sites.csv')):
                xce_object.proasis_lead_label = str('Create lead from pandda sites...')
                xce_object.proasis_lead_function = partial(Proasis().add_lead, xce_object=xce_object)
                xce_object.needtoaddleads = True
            else:

                xce_object.proasis_lead_label = str('Site info not found... please run pandda analyse before adding lead')
                xce_object.proasis_lead_function = lambda:Proasis().blank()
                xce_object.needtoaddleads = True

        # Hit details or add hits (refined) in menu
        counter = 0
        try:
            if hasattr(xce_object,'needtoaddleads'):
                xce_object.proasis_hits_label = 'To add hits: add a lead, and then reload XCE'
                xce_object.proasis_hits_function = lambda:Proasis().blank()
            else:
                # count the number of hits in proasis if they exist (from soakDB)
                for row in c.execute('SELECT proasisID FROM proasis'):
                    counter += 1
                no_hits = counter
                # display no of hits (proasis) in menu, no action if clicked
                xce_object.proasis_hits_label = str('Hits in proasis: ' + str(no_hits))
                xce_object.proasis_hits_function = lambda:Proasis().blank()
        # otherwise, try to add hits to proasis (if there are no hits, the job will still run and hits will be added as
        # they are refined - i.e. when refine.bound.pdb file is detected for a refinement detailed in soakDB)
        except:
            xce_object.proasis_hits_label = str('Attempt to add refined hits to proasis...')
            xce_object.proasis_hits_function = lambda: Proasis().add_hits(xce_object)


    def create_project(self, name, xce_object):
        # make relevant project directory in proasis LabXChem folder
        print(str('Making Proasis project directory: ' + str(
            'mkdir ' + os.path.join(xce_object.proasis_directory, 'LabXChem', name))))
        os.system(str('mkdir ' + os.path.join(xce_object.proasis_directory, 'LabXChem', name)))
        perm_string = str('chmod u=rwx,g=rwx,o=r ' + os.path.join(xce_object.proasis_directory, 'LabXChem', name))
        os.system(perm_string)
        # make reference file directory in project directory
        os.system(str('mkdir ' + os.path.join(xce_object.proasis_directory, 'LabXChem', name, 'reference')))
        perm_string = str('chmod u=rwx,g=rwx,o=r ' + os.path.join(xce_object.proasis_directory, 'LabXChem',
                                                                  name, 'reference'))
        os.system(perm_string)
        # create a temporary job to add the project in proasis schedule
        temp_job = open(os.path.join(xce_object.proasis_directory, 'Scripts/scheduled_jobs/temp_jobs',
                                     str(name + '.sh')), 'w')
        perm_string = str(
            'chmod u=rwx,g=rwx,o=r ' + os.path.join(xce_object.proasis_directory, 'Scripts/scheduled_jobs/temp_jobs',
                                                    str(name
                                                        + '.sh')))
        os.system(perm_string)
        job_string = str('/usr/local/Proasis2/utils/addnewproject.py -q OtherClasses -p ' + name)
        temp_job.write(str(job_string))
        temp_job.close()
        print('==> XCE: Added job to queue... please restart XCE...')
        print('==> XCE: Exiting!')
        xce_object.quit_xce()

    def add_lead(self, xce_object):
        if hasattr(xce_object, 'leadadded'):
            print('The lead has already been added,, please reload XCE to check status')
        # in case directories don't exist...
        print(str('Making Proasis project directory: ' + str(
            'mkdir ' + os.path.join(xce_object.proasis_directory, 'LabXChem', xce_object.proasis_name))))
        os.system(str('mkdir ' + os.path.join(xce_object.proasis_directory, 'LabXChem', xce_object.proasis_name)))
        perm_string = str(
            'chmod u=rwx,g=rwx,o=r ' + os.path.join(xce_object.proasis_directory, 'LabXChem', xce_object.proasis_name))
        os.system(perm_string)
        # make reference file directory in project directory
        os.system(str('mkdir ' + os.path.join(xce_object.proasis_directory, 'LabXChem', xce_object.proasis_name, 'reference')))
        perm_string = str('chmod u=rwx,g=rwx,o=r ' + os.path.join(xce_object.proasis_directory, 'LabXChem',
                                                                  xce_object.proasis_name, 'reference'))

        # copy pandda_analyse_sites.csv to proasis directory for lead build
        os.system(str('cp ' + str(os.path.join(xce_object.panddas_directory, 'analyses/pandda_analyse_sites.csv')) + ' ' +
                      str(os.path.join(xce_object.proasis_directory, 'LabXChem', xce_object.proasis_name, 'reference'))))
        # copy reference pdb (from pandda directory to make sure same as in sites file)
        os.system(str('cp ' + str(os.path.join(xce_object.panddas_directory, 'reference/reference.pdb')) + ' ' +
                      str(os.path.join(xce_object.proasis_directory, 'LabXChem', xce_object.proasis_name, 'reference'))))
        # open a temporary job file to write to for proasis scheduling
        temp_job = open(
            os.path.join(xce_object.proasis_directory, 'Scripts/scheduled_jobs/temp_jobs', str(xce_object.proasis_name + '.sh')),
            'w')
        # change file permissions of job
        perm_string = str(
            'chmod u=rwx,g=rwx,o=r ' + os.path.join(xce_object.proasis_directory, 'Scripts/scheduled_jobs/temp_jobs',
                                                    str(xce_object.proasis_name + '.sh')))
        os.system(perm_string)
        # string to add leads in temp job file
        job_string = str('/dls/science/groups/proasis/Scripts/generate_leads.py -n ' + xce_object.proasis_name
                         + ' -r '
                         + str(
            os.path.join(xce_object.proasis_directory, 'LabXChem', xce_object.proasis_name, 'reference/reference.pdb'))
                         + ' -p '
                         + str(os.path.join(xce_object.proasis_directory, 'LabXChem', xce_object.proasis_name,
                                            'reference/pandda_analyse_sites.csv'))
                         + ' -d '
                         + str(xce_object.current_directory))

        temp_job.write(str(job_string))
        temp_job.close()
        # remove option from menu so lead can't be added multiple times
        xce_object.leadadded=True
        print('==> XCE: Added job to queue... please restart XCE...')
        print('==> XCE: Exiting!')
        xce_object.quit_xce()

    def add_hits(self, xce_object):
        # open the list of pernament jobs to append
        perm_job = open(os.path.join(xce_object.proasis_directory, 'Scripts/scheduled_jobs/test.sh'), 'a')
        # string for job to add and update hits in proasis
        job_string = (str(os.path.join(xce_object.proasis_directory, 'Scripts/populate_hits.py') + ' -d ' +
                          xce_object.current_directory + ' > ' + os.path.join(xce_object.proasis_directory,
                                                                        'Scripts/scheduled_jobs_logs',
                                                                        str(xce_object.proasis_name + '_proasis.out'))))
        perm_job.write(job_string)
        perm_job.close()
        print('==> XCE: Added job to queue... please restart XCE...')
        print('==> XCE: Exiting!')
        xce_object.quit_xce()


    def blank(self):
        print('==> XCE: This does nothing...')
