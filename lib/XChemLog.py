import os,glob
import sys
from datetime import datetime

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

class startLog:

    def __init__(self,logfile):
        self.logfile=logfile

    def create_logfile(self):
        if not os.path.isfile(self.logfile):
            os.system('touch '+self.logfile)
