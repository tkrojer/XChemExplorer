# last edited: 10/01/2017, 15:00

import os,sys
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))


def copy_files(htmlDir):
        if not os.path.isdir('zenodo'):
            os.mkdir('zenondo')
        os.chdir(os.path.join(htmlDir,'zenodo'))

        print 'copying compoundImages'
        os.system('/bin/cp ../compoundImages/* .')

        print 'copying pdbs'
        os.system('/bin/cp ../pdbs/* .')

        print 'copying maps'
        os.system('/bin/cp ../maps/* .')

        print 'copying residueplots'
        os.system('/bin/cp ../residueplots/* .')

        print 'copying mapImages'
        os.system('/bin/cp ../mapImages/* .')

        print 'copying icbs'
        os.system('/bin/cp ../icbs/* .')

def edit_index_html(htmlDir):
    os.chdir(os.path.join(htmlDir,'zenodo'))

    out=''
    for line in open('../index.html'):
        if 'compoundImages/' in line:
            line.replace('compoundImages/','')
        if 'pdbs/' in line:
            line.replace('pdbs/','')
        if 'maps/' in line:
            line.replace('maps/','')
        if 'residueplots/' in line:
            line.replace('residueplots/','')
        if 'mapImages/' in line:
            line.replace('mapImages/','')
        if 'icbs/' in line:
            line.replace('icbs/','')
        out+=line

    print 'writing index.html to',os.path.join(htmlDir,'zenodo')
    f=open('index.html','w')
    f.write(out)
    f.close()


if __name__=='__main__':
    htmlDir=sys.argv[1]

    if os.path.isfile(htmlDir):
        os.chdir(htmlDir)
        copy_files(htmlDir)
        edit_index_html(htmlDir)
    else:
        print 'sorry, no html export directory'



