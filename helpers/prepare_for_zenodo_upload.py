# last edited: 11/01/2017, 15:00

import os,sys
import glob
#sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))


def copy_files(htmlDir):
        os.chdir(htmlDir)
        if not os.path.isdir('zenodo'):
            os.mkdir('zenodo')
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

        print 'copying css'
        os.system('/bin/cp ../css/* .')

        print 'copying js'
        os.system('/bin/cp ../js/* .')

def edit_index_html(htmlDir,uploadID):
    os.chdir(os.path.join(htmlDir,'zenodo'))

    out=''
    for line in open('../index.html'):
        if 'compoundImages/' in line:
            line=line.replace('compoundImages/','https://zenodo.org/record/'+uploadID+'/files/')
        if 'pdbs/' in line:
            line=line.replace('pdbs/','https://zenodo.org/record/'+uploadID+'/files/')
        if 'maps/' in line:
            line=line.replace('maps/','https://zenodo.org/record/'+uploadID+'/files/')
        if 'residueplots/' in line:
            line=line.replace('residueplots/','https://zenodo.org/record/'+uploadID+'/files/')
        if 'mapImages/' in line:
            line=line.replace('mapImages/','https://zenodo.org/record/'+uploadID+'/files/')
        if 'icbs/' in line:
            line=line.replace('icbs/','https://zenodo.org/record/'+uploadID+'/files/')
        if 'css/' in line:
            line=line.replace('css/','https://zenodo.org/record/'+uploadID+'/files/')
        if 'js/' in line:
            line=line.replace('js/','https://zenodo.org/record/'+uploadID+'/files/')
        out+=line

    print 'writing index.html to',os.path.join(htmlDir,'zenodo')
    f=open('0_index.html','w')
    f.write(out)
    f.close()

def edit_icb_html_files(htmlDir,uploadID):
    # https://zenodo.org/record/48768/files/README.txt
    os.chdir(os.path.join(htmlDir,'zenodo'))
    for file in glob.glob('*event*html'):
        out=''
        for line in open(file):
            if 'src="http://molsoft.com/lib/acticm.js">' in line:
                line=line.replace('src="http://molsoft.com/lib/acticm.js">','src="https://molsoft.com/lib/acticm.js">')
            if '../mapImages/' in line:
                line=line.replace('../mapImages/','https://zenodo.org/record/'+uploadID+'/files/')
            if '../compoundImages/' in line:
                line=line.replace('../compoundImages/','https://zenodo.org/record/'+uploadID+'/files/')
            if ' act.projectFile = "' in line:
                line=line.replace(' act.projectFile = "',' act.projectFile = "https://zenodo.org/record/'+uploadID+'/files/')
            out+=line
        print 'updating',file
        f=open(file,'w')
        f.write(out)
        f.close()


if __name__=='__main__':
    htmlDir=sys.argv[1]
    uploadID=sys.argv[2]

    if os.path.isdir(htmlDir):
#        copy_files(htmlDir)

        # MISSING: copy .cif files!!!

        edit_index_html(htmlDir,uploadID)
        edit_icb_html_files(htmlDir,uploadID)
    else:
        print 'sorry, no html export directory'



