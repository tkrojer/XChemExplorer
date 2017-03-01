# last edited: 01/03/2017, 15:00

import os,glob
import sys

def relink_interesting_datasets():
    dirList=[]
    for dirs in glob.glob('*'):
        if os.path.isdir(dirs):
            dirList.append(dirs)

    for dirs in dirList:
        print '/bin/rm -fr %s' %dirs
        os.system('/bin/rm -fr %s' %dirs)
        print 'ln -s ../processed_datasets/%s .' %dirs
        os.system('ln -s ../processed_datasets/%s .' %dirs)

def check_interesting_datasets(panddaDir):

    os.chdir(os.path.join(panddaDir,'interesting_datasets'))
    for xtal in glob.glob('*'):
        if not os.path.islink(xtal):
#            print xtal
#            print os.path.join(panddaDir,'interesting_datasets',xtal,'modelled_structures',xtal+'-pandda-model.pdb')
            if os.path.isfile(os.path.join(panddaDir,'interesting_datasets',xtal,'modelled_structures',xtal+'-pandda-model.pdb')):
                print 'mv %s %s' %(os.path.join(panddaDir,'interesting_datasets',xtal,'modelled_structures'),os.path.join(panddaDir,'processed_datasets',xtal))
                print 'cp %s %s' %(os.path.join(panddaDir,'interesting_datasets',xtal,'*params'),os.path.join(panddaDir,'processed_datasets',xtal))




if __name__=='__main__':
    panddaDir=sys.argv[1]
    check_interesting_datasets(panddaDir)
