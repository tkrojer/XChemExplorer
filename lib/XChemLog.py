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
            message='==> XCE: creating new logfile for the current XChemExplorer session: '+self.logfile
        else:
            message='==> XCE: writing into existing logfile for current XChemExplorer session: '+self.logfile
        updateLog(self.logfile).insert(message)

class updateLog:

    def __init__(self,logfile):
        self.logfile=open(logfile, "w")

    def insert(self,message):
        present_time=datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')
        print( str(present_time)+'  ==> XCE: '+message, file = self.logfile)
        print(message)
