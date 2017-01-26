# last edited: 26/01/2017, 15:00

from __future__ import print_function
import os,glob
import sys
from datetime import datetime

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

class startLog:

    def __init__(self,logfile):
        self.logfile=logfile

    def create_logfile(self,version):

        pasteVersion=version
        for i in range(0,20-len(version)):
            pasteVersion+=' '

        message = (
            '\n\n'
            '     ######################################################################\n'
            '     #                                                                    #\n'
            '     # XCHEMEXPLORER - multi dataset analysis                             #\n'
            '     #                                                                    #\n'
            '     # Version: %s                                      #\n' %pasteVersion+
            '     #                                                                    #\n'
            '     # Date:                                                              #\n'
            '     #                                                                    #\n'
            '     # Author: Tobias Krojer, Structural Genomics Consortium, Oxford, UK  #\n'
            '     #         tobias.krojer@sgc.ox.ac.uk                                 #\n'
            '     #                                                                    #\n'
            '     ######################################################################\n'
            '\n'
        )

        if not os.path.isfile(self.logfile):
            os.system('touch '+self.logfile)
            message+='creating new logfile for the current XChemExplorer ('+version+') session:\n'+self.logfile+'\n'

        else:
            message+='writing into existing logfile for current XChemExplorer ('+version+') session:\n'+self.logfile+'\n'
        updateLog(self.logfile).insert(message)

class updateLog:

    def __init__(self,logfile):
        self.logfile=open(logfile, "a")

    def insert(self,message):
        present_time=datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S')
        print( str(present_time)+' ==> XCE: '+message, file = self.logfile)
        print('==> XCE: '+message)
