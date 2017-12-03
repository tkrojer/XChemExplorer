import os

def run_script_locally():
    os.system(cmd)

def run_script_remotely(cmd):
    print 'hallo'

def pdb_extract(outRoot,refSoft,pdb,log,data_integration_software,phasing_software):
    if os.path.isdir('/dls'):
        pdb_extract_init='source /dls/science/groups/i04-1/software/pdb-extract-prod/setup.sh\n'
        pdb_extract_init+='/dls/science/groups/i04-1/software/pdb-extract-prod/bin/pdb_extract'
    else:
        pdb_extract_init='source '+os.path.join(os.getenv('XChemExplorer_DIR'),'pdb_extract/pdb-extract-prod/setup.sh')+'\n'
        pdb_extract_init+=+os.path.join(os.getenv('XChemExplorer_DIR'),'pdb_extract/pdb-extract-prod/bin/pdb_extract')

    Cmd = ( pdb_extract_init+
        ' -r {0!s}'.format(refSoft)+
        ' -iPDB {0!s}'.format(pdb)+
        ' -i {0!s}'.format(data_integration_software) +
        ' -p {0!s}'.format(phasing_software) +
        ' -e MR'
        ' -s AIMLESS'
        ' -iLOG %s'         %log+
        ' -iENT data_template.cif '
        ' -o {0!s}.mmcif > {1!s}.mmcif.log'.format(outRoot, outRoot)       )
    return Cmd


def sf_convert(outRoot):
        if os.path.isdir('/dls'):
            pdb_extract_init='source /dls/science/groups/i04-1/software/pdb-extract-prod/setup.sh\n'
            pdb_extract_init+='/dls/science/groups/i04-1/software/pdb-extract-prod/bin/sf_convert'
        else:
            pdb_extract_init='source '+os.path.join(os.getenv('XChemExplorer_DIR'),'pdb_extract/pdb-extract-prod/setup.sh')+'\n'
            pdb_extract_init+=+os.path.join(os.getenv('XChemExplorer_DIR'),'pdb_extract/pdb-extract-prod/bin/sf_convert')

        Cmd = ( pdb_extract_init+
                ' -o mmcif'
                ' -sf %s' %mtzin+
                ' -out {0!s}_sf.mmcif  > {1!s}.sf_mmcif.log'.format(outRoot, outRoot) )
        return Cmd