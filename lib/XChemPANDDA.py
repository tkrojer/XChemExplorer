import os

class PANDDAs(object):
    def __init__(self,initial_model_directory,panddas_directory):
        self.initial_model_directory=initial_model_directory
        self.panddas_directory=panddas_directory

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
                'pandda.analyse data_dirs="'+os.path.join(self.initial_model_directory,'*','Dimple','dimple')+'"'
                ' pdb_style=final.pdb out_dir='+self.panddas_directory+'\n'
                '\n'
                '/bin/rm PANDDA_RUN_IN_PROGRESS\n'
            )

            f = open('pandda.sh','w')
            f.write(Cmds)
            f.close()
#            os.system('chmod +x pandda.sh')
#            os.system('./pandda.sh &')


