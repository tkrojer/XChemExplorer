import os,sys
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

#from XChemUtils import parse
#from XChemUtils import external_software
#from XChemUtils import helpers
#import XChemThread
import XChemDB
#import XChemDialogs
#import XChemPANDDA
#import XChemToolTips
#import XChemMain

if __name__=='__main__':
    db_file=sys.argv[1]
    xtal=sys.argv[2]
    inital_model_directory=sys.argv[3]

    print 'XCE path',os.path.join(os.getenv('XChemExplorer_DIR'),'lib')


#    db=XChemDB.data_source(db_file)
    if os.path.isfile(os.path.join(inital_model_directory,xtal,'dimple.pdb')):
        print os.path.join(inital_model_directory,xtal,'dimple.pdb')
