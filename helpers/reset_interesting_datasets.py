#!/usr/bin/python
# last edited: 01/03/2017, 15:00

import os,sys

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
    for folder in glob.glob('*'):
        if not os.path.islink(folder):
            print 'not symlink',folder


if __name__=='__main__':
    panddaDir=sys.argv[1]
    check_interesting_datasets(panddaDir)
