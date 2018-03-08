import os, sys

# set XCE environment variables
os.environ['XChemExplorer_DIR']=os.path.join(os.getcwd(), 'XChemExplorer/')

print('XCE path = ' + os.environ['XChemExplorer_DIR'])

# Append all directories to path
sys.path.append(os.environ['XChemExplorer_DIR'])

# print PATH to check
print('Current path:- ')
for p in sys.path:
    print(p)

# try to import all of XCE stuff
from XChemExplorer import *
print('Successfully imported XCE... compile was succesful')
