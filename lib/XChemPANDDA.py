import os

class PANDDAs(object):
    def __init__(self,pandda_params):
        self.data_directory=pandda_params['data_dir']
        self.panddas_directory=pandda_params['out_dir']
        self.nproc=pandda_params['nproc']
        self.submit_mode=pandda_params['submit_mode']
        self.xtalform=pandda_params['xtalform']

    def run_pandda_analyse(self):
        if os.path.isfile(os.path.join(self.panddas_directory,'PANDDA_RUN_IN_PROGRESS')):
            return None
        else:
            os.chdir(self.panddas_directory)
            os.system('touch PANDDA_RUN_IN_PROGRESS')

            Cmds = (
                '#!'+os.getenv('SHELL')+'\n'
                '\n'
                'cd '+self.panddas_directory+'\n'
                '\n'
                'pandda.analyse data_dirs="'+self.data_directory+'"'
                ' pdb_style=final.pdb out_dir='+self.panddas_directory+'\n'
                '\n'
                '/bin/rm PANDDA_RUN_IN_PROGRESS\n'
            )
            print Cmds
#            f = open('pandda.sh','w')
#            f.write(Cmds)
#            f.close()
#            os.system('chmod +x pandda.sh')
#            os.system('./pandda.sh &')

    def launch_pandda_inspect(self):
        print 'hallo'