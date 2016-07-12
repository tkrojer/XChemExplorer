from __future__ import print_function
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

class updateLog:

    def __init__(self,logfile):
        self.logfile=open(logfile, "w")

    def insert(self,message):
        print(message, file = self.logfile)