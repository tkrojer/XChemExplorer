### RACHAEL'S PROASIS STUFF ###

        self.proasis_directory = '/dls/science/groups/proasis/'

        # function for adding a new project
        def create_project(name):
            # make relevant project directory in proasis LabXChem folder
            print(str('Making Proasis project directory: ' + str(
                'mkdir ' + os.path.join(self.proasis_directory, 'LabXChem', name))))
            os.system(str('mkdir ' + os.path.join(self.proasis_directory, 'LabXChem', name)))
            perm_string = str('chmod u=rwx,g=rwx,o=r ' + os.path.join(self.proasis_directory, 'LabXChem', name))
            os.system(perm_string)
            # make reference file directory in project directory
            os.system(str('mkdir ' + os.path.join(self.proasis_directory, 'LabXChem', name, 'reference')))
            perm_string = str('chmod u=rwx,g=rwx,o=r ' + os.path.join(self.proasis_directory, 'LabXChem',
                                                                      name, 'reference'))
            os.system(perm_string)
            # create a temporary job to add the project in proasis schedule
            temp_job = open(os.path.join(self.proasis_directory, 'Scripts/scheduled_jobs/temp_jobs',
                                         str(name + '.sh')), 'w')
            perm_string = str(
                'chmod u=rwx,g=rwx,o=r ' + os.path.join(self.proasis_directory, 'Scripts/scheduled_jobs/temp_jobs',
                                                        str(name
                                                            + '.sh')))
            os.system(perm_string)
            job_string = str('/usr/local/Proasis2/utils/addnewproject.py -q OtherClasses -p ' + name)
            temp_job.write(str(job_string))
            temp_job.close()

        # add proasis menu to main menu
        self.proasis_menu = menu_bar.addMenu('Proasis')
        # connect to soakDB to get proasis info
        if os.path.isfile(os.path.join(self.database_directory, self.data_source_file)):
            conn = sqlite3.connect(os.path.join(self.database_directory, self.data_source_file))
            c = conn.cursor()

        # Project details or add project in menu
        counter = 0
        try:
            # get protein name from soakDB - this will be the proasis project name
            for row in c.execute('SELECT Protein FROM soakDB;'):
                counter += 1
                # if there is only one protein name in soakDB - all is good - happy days
                if counter == 1:
                    self.proasis_name = str(row[0])
                # otherwise - give a warning
                # TODO: If this is actually ever encountered - deal with it. Should be fine.
                if counter > 1:
                    print('WARNING: More than one protein name found (proasis)')
                # If the project directory already exists in proasis dir, project should exist in proasis
                if os.path.isdir(os.path.join('/dls/science/groups/proasis/LabXChem/', self.proasis_name)):
                    # show project name in menu (no action when clicked)
                    self.proasis_project = QtGui.QAction(str('Project Name: ' + self.proasis_name), self.window)
                    self.proasis_menu.addAction(self.proasis_project)
        # should catch if project doesnt exitst
        except (AttributeError, UnboundLocalError):
            self.update_log.insert('cannot find %s' % os.path.join(self.database_directory, self.data_source_file))
        # except UnboundLocalError:
        #     self.update_log.insert('cannot find %s' % os.path.join(self.database_directory, self.data_source_file))
        except:
            # option to create project, action = create_project()
            self.proasis_project = QtGui.QAction(str('Create Project for ' + self.proasis_name + '...'),
                                                 self.window)
            self.proasis_project.triggered.connect(lambda: create_project(self.proasis_name))
            self.proasis_menu.addAction(self.proasis_project)

        # Lead details or add lead in menu
        counter = 0
        try:
            # check if there is a lead (from soakDB)
            for row in c.execute('SELECT proasisID FROM proasisLead'):
                counter += 1
                if counter == 1:
                    # If so, display id of lead in menu, no action if clicked
                    self.proasis_lead = QtGui.QAction(str('Lead ID: ' + str(row[0])), self.window)
                    self.proasis_menu.addAction(self.proasis_lead)
                # otherwise, if you can find the pandda_analyse_sites.csv file, allow lead to be added
                elif os.path.isfile(os.path.join(self.panddas_directory, 'analyses/pandda_analyse_sites.csv')):
                    self.proasis_lead = QtGui.QAction(str('Create lead from pandda sites...'), self.window)
                    self.proasis_menu.addAction(self.proasis_lead)
        # If no lead or sites file, error message. No action on click
        except:
            self.proasis_lead = QtGui.QAction(str('Site info not found... '
                                                  'please run pandda analyse before adding lead'), self.window)
            self.proasis_lead.triggered.connect(lambda: self.add_lead())
            self.proasis_menu.addAction(self.proasis_lead)

        # Hit details or add hits (refined) in menu
        counter = 0
        try:
            # count the number of hits in proasis if they exist (from soakDB)
            for row in c.execute('SELECT proasisID FROM proasis'):
                counter += 1
            no_hits = counter
            # display no of hits (proasis) in menu, no action if clicked
            self.proasis_hits = QtGui.QAction(str('Hits in proasis: ' + str(no_hits)), self.window)
            self.proasis_menu.addAction(self.proasis_hits)
        # otherwise, try to add hits to proasis (if there are no hits, the job will still run and hits will be added as
        # they are refined - i.e. when refine.bound.pdb file is detected for a refinement detailed in soakDB)
        except:
            self.proasis_hits = QtGui.QAction(str('Attempt to add refined hits to proasis...'), self.window)
            self.proasis_hits.triggered.connect(lambda: self.add_hits())
            self.proasis_menu.addAction(self.proasis_hits)

        ### end of proasis shit ###


def add_lead(self):
    # in case directories don't exist...
    print(str('Making Proasis project directory: ' + str(
        'mkdir ' + os.path.join(self.proasis_directory, 'LabXChem', self.proasis_name))))
    os.system(str('mkdir ' + os.path.join(self.proasis_directory, 'LabXChem', self.proasis_name)))
    perm_string = str(
        'chmod u=rwx,g=rwx,o=r ' + os.path.join(self.proasis_directory, 'LabXChem', self.proasis_name))
    os.system(perm_string)
    # make reference file directory in project directory
    os.system(str('mkdir ' + os.path.join(self.proasis_directory, 'LabXChem', self.proasis_name, 'reference')))
    perm_string = str('chmod u=rwx,g=rwx,o=r ' + os.path.join(self.proasis_directory, 'LabXChem',
                                                              self.proasis_name, 'reference'))

    # copy pandda_analyse_sites.csv to proasis directory for lead build
    os.system(str('cp ' + str(os.path.join(self.panddas_directory, 'analyses/pandda_analyse_sites.csv')) + ' ' +
                  str(os.path.join(self.proasis_directory, 'LabXChem', self.proasis_name, 'reference'))))
    # copy reference pdb (from pandda directory to make sure same as in sites file)
    os.system(str('cp ' + str(os.path.join(self.panddas_directory, 'reference/reference.pdb')) + ' ' +
                  str(os.path.join(self.proasis_directory, 'LabXChem', self.proasis_name, 'reference'))))
    # open a temporary job file to write to for proasis scheduling
    temp_job = open(
        os.path.join(self.proasis_directory, 'Scripts/scheduled_jobs/temp_jobs', str(self.proasis_name + '.sh')),
        'w')
    # change file permissions of job
    perm_string = str(
        'chmod u=rwx,g=rwx,o=r ' + os.path.join(self.proasis_directory, 'Scripts/scheduled_jobs/temp_jobs',
                                                str(self.proasis_name + '.sh')))
    os.system(perm_string)
    # string to add leads in temp job file
    job_string = str('/dls/science/groups/proasis/Scripts/generate_leads.py -n ' + self.proasis_name
                     + ' -r '
                     + str(
        os.path.join(self.proasis_directory, 'LabXChem', self.proasis_name, 'reference/reference.pdb'))
                     + ' -p '
                     + str(os.path.join(self.proasis_directory, 'LabXChem', self.proasis_name,
                                        'reference/pandda_analyse_sites.csv'))
                     + ' -d '
                     + str(self.current_directory))

    temp_job.write(str(job_string))
    temp_job.close()
    # remove option from menu so lead can't be added multiple times
    self.proasis_lead.setVisible(False)


def add_hits(self):
    # open the list of pernament jobs to append
    perm_job = open(os.path.join(self.proasis_directory, 'Scripts/scheduled_jobs/test.sh'), 'a')
    # string for job to add and update hits in proasis
    job_string = (str(os.path.join(self.proasis_directory, 'Scripts/populate_hits.py') + ' -d ' +
                      self.current_directory + ' > ' + os.path.join(self.proasis_directory,
                                                                    'Scripts/scheduled_jobs_logs',
                                                                    str(self.proasis_name + '_proasis.out'))))
    perm_job.write(job_string)
    perm_job.close()
    # remove option from menu so job is not added multiple times
    self.proasis_hits.setVisible(False)