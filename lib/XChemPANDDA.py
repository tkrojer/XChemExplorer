import os, sys, glob
from datetime import datetime
from PyQt4 import QtGui, QtCore


class run_pandda_analyse(QtCore.QThread):

    def __init__(self,pandda_params):
        QtCore.QThread.__init__(self)
        self.data_directory=pandda_params['data_dir']
        self.panddas_directory=pandda_params['out_dir']
        self.nproc=pandda_params['nproc']
        self.submit_mode=pandda_params['submit_mode']
        self.min_build_datasets=pandda_params['min_build_datasets']

    def run(self):
        if os.path.isfile(os.path.join(self.panddas_directory,'pandda.running')):
            return None
        else:
            os.chdir(self.panddas_directory)
            Cmds = (
                '#!'+os.getenv('SHELL')+'\n'
                '\n'
                'cd '+self.panddas_directory+'\n'
                '\n'
                'pandda.analyse data_dirs="'+self.data_directory+'"'
                ' pdb_style=final.pdb out_dir='+self.panddas_directory+
                ' min_build_datasets='+self.min_build_datasets+
                ' cpus='+self.nproc+'\n'
                )
            print Cmds
            f = open('pandda.sh','w')
            f.write(Cmds)
            f.close()
            if self.submit_mode=='local machine':
                os.system('chmod +x pandda.sh')
                os.system('./pandda.sh &')
            else:
                os.system('qsub pandda.sh')

