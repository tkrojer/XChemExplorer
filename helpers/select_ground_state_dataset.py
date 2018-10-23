# last edited: 03/08/2017, 15:00

import os,sys
import glob
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

from XChemUtils import parse


# - select datasets with highest resolution
# - select only those without an event map
# - take the one with the lowest Rfree


def find_highest_resolution_datasets(panddaDir):
    found=False
    datasetList=[]

    datasetList = ['TCRUFPPS - x0045', 'TCRUFPPS - x0053', 'TCRUFPPS - x0063', 'TCRUFPPS - x0064', 'TCRUFPPS - x0065',
                   'TCRUFPPS - x0073', 'TCRUFPPS - x0074', 'TCRUFPPS - x0076', 'TCRUFPPS - x0079', 'TCRUFPPS - x0082',
                   'TCRUFPPS - x0087', 'TCRUFPPS - x0088', 'TCRUFPPS - x0092', 'TCRUFPPS - x0096', 'TCRUFPPS - x0097',
                   'TCRUFPPS - x0103', 'TCRUFPPS - x0108', 'TCRUFPPS - x0111', 'TCRUFPPS - x0112', 'TCRUFPPS - x0124',
                   'TCRUFPPS - x0125', 'TCRUFPPS - x0126', 'TCRUFPPS - x0128', 'TCRUFPPS - x0129', 'TCRUFPPS - x0134',
                   'TCRUFPPS - x0135', 'TCRUFPPS - x0136', 'TCRUFPPS - x0138', 'TCRUFPPS - x0141', 'TCRUFPPS - x0142',
                   'TCRUFPPS - x0144', 'TCRUFPPS - x0145', 'TCRUFPPS - x0147', 'TCRUFPPS - x0152', 'TCRUFPPS - x0154',
                   'TCRUFPPS - x0155', 'TCRUFPPS - x0159', 'TCRUFPPS - x0160', 'TCRUFPPS - x0161', 'TCRUFPPS - x0162',
                   'TCRUFPPS - x0163', 'TCRUFPPS - x0164', 'TCRUFPPS - x0165', 'TCRUFPPS - x0169', 'TCRUFPPS - x0176',
                   'TCRUFPPS - x0180', 'TCRUFPPS - x0182', 'TCRUFPPS - x0183', 'TCRUFPPS - x0184', 'TCRUFPPS - x0187',
                   'TCRUFPPS - x0192', 'TCRUFPPS - x0194', 'TCRUFPPS - x0195', 'TCRUFPPS - x0196', 'TCRUFPPS - x0200',
                   'TCRUFPPS - x0203', 'TCRUFPPS - x0204', 'TCRUFPPS - x0208', 'TCRUFPPS - x0209', 'TCRUFPPS - x0216']

    # Issue is likely to be in belwo code. Assess on recent run of pandda prerun

    # for logFile in glob.glob(os.path.join(panddaDir,'logs','*.log')):
    #     for n,line in enumerate(open(logFile)):
    #
    #         if line.startswith('Statistical Electron Density Characterisation') and len(line.split()) == 6:
    #             resolution=line.split()[5]
    #             found=True
    #             foundLine=n
    #         if found and n>=foundLine+3:
    #             if line.startswith('---'):
    #                 break
    #             else:
    #                 tmpLine=line.replace(' ','').replace('\t','').replace('\n','').replace('\r','')
    #                 for item in tmpLine.split(','):
    #                     if item != '':
    #                         datasetList.append(item)
    print datasetList
    return datasetList

def get_datasets_without_event_map(panddaDir,datasetList):
    datasetListwithoutEvent=[]
    for dataset in datasetList:
        noEvent=True
        for files in glob.glob(os.path.join(panddaDir,'processed_datasets',dataset,'*')):
            if 'event' in files:
                noEvent=False
                break
        if noEvent:
            datasetListwithoutEvent.append(dataset)
    print datasetListwithoutEvent
    return datasetListwithoutEvent

def select_dataset_with_lowest_Rfree(panddaDir,datasetListwithoutEvent):
    datasetList=[]
    lowestRfree=''
    for dataset in datasetListwithoutEvent:
        if os.path.isfile(os.path.join(panddaDir,'processed_datasets',dataset,dataset+'-pandda-input.pdb')):
            stats=parse().PDBheader(os.path.join(panddaDir,'processed_datasets',dataset,dataset+'-pandda-input.pdb'))
            Rfree=stats['Rfree']
            try:
                print dataset,Rfree,stats['ResolutionHigh']
                datasetList.append([dataset,float(Rfree)])
            except ValueError:
                pass
    if datasetList != []:
        lowestRfree=min(datasetList,key=lambda x: x[1])[0]
    return lowestRfree

def link_pdb_mtz_files(panddaDir,lowestRfree):
    targetDir='/'.join(panddaDir.split('/')[:len(panddaDir.split('/'))-1])
    panddaFolder=panddaDir.split('/')[len(panddaDir.split('/'))-1]
    print targetDir
    print panddaFolder
    os.chdir(targetDir)
    if os.path.isfile(os.path.join(panddaDir,'processed_datasets',lowestRfree,lowestRfree+'-pandda-input.pdb')):
        os.system('/bin/rm %s-ground-state.pdb 2> /dev/null' %lowestRfree)
        os.symlink(os.path.join(panddaFolder,'processed_datasets',lowestRfree,lowestRfree+'-pandda-input.pdb'),lowestRfree+'-ground-state.pdb')
    if os.path.isfile(os.path.join(panddaDir,'processed_datasets',lowestRfree,lowestRfree+'-pandda-input.mtz')):
        os.system('/bin/rm %s-ground-state.free.mtz 2> /dev/null' %lowestRfree)
        os.symlink(os.path.join(panddaFolder,'processed_datasets',lowestRfree,lowestRfree+'-pandda-input.mtz'),lowestRfree+'-ground-state.free.mtz')
    if os.path.isfile(os.path.join(panddaDir,'processed_datasets',lowestRfree,lowestRfree+'-ground-state-mean-map.native.ccp4')):
        os.system('/bin/rm %s-ground-state-mean-map.native.ccp4 2> /dev/null' %lowestRfree)
        os.symlink(os.path.join(panddaFolder,'processed_datasets',lowestRfree,lowestRfree+'-ground-state-mean-map.native.ccp4'),lowestRfree+'-ground-state-mean-map.native.ccp4')



if __name__=='__main__':
    panddaDir=sys.argv[1]
    datasetList=find_highest_resolution_datasets(panddaDir)
    datasetListwithoutEvent=get_datasets_without_event_map(panddaDir,datasetList)
    lowestRfree=select_dataset_with_lowest_Rfree(panddaDir,datasetListwithoutEvent)
    link_pdb_mtz_files(panddaDir,lowestRfree)
