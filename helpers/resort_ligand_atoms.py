# last edited: 07/06/2017, 17:00

import os,sys
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

def get_atom_order_of_ensemble_model(ensembleModel):

    ligandResnames = ['LIG','FRS','DRG','UNL']

    ensembleLIGdir = {}

    for line in open(ensembleModel):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atomName_line=  str(line[11:16])
            resname_line=   str(line[17:20])
            chainID_line=   str(line[20:23])
            resseq_line=    str(line[23:26])
            residueID=resname_line.replace(' ','')+'-'+chainID_line.replace(' ','')+'-'+resseq_line.replace(' ','')
            if resname_line in ligandResnames:
                if residueID not in ensembleLIGdir:
                    ensembleLIGdir[residueID]=[]
                ensembleLIGdir[residueID].append(atomName_line)

    return ensembleLIGdir

def resort_ligand_atoms_in_refined_model(refinedModel,ensembleLIGdir):
    ligandResnames = ['LIG','FRS','DRG','UNL']
    out=''
    refineLIGDir={}
    LIGlineAdditional=''
    LIGlist=[]
    for line in open(refinedModel):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            resname_line=   str(line[17:20])
            chainID_line=   str(line[20:23])
            resseq_line=    str(line[23:26])
            residueID=resname_line.replace(' ','')+'-'+chainID_line.replace(' ','')+'-'+resseq_line.replace(' ','')
            if residueID in ensembleLIGdir:
                if residueID not in refineLIGDir:
                    refineLIGDir[residueID]=[]
                else:
                    # if a ligand was for modelled after pandda.inspect
                    LIGlineAdditional+=line
                refineLIGDir[residueID].append(line)


    LIGline=''
    for ligand in ensembleLIGdir:
        for atom in ensembleLIGdir[ligand]:
            for line in refineLIGDir[ligand]:
                atomName_line=str(line[11:16])
                if atomName_line==atom:
                    LIGline+=line
                    break

    for line in open(refinedModel):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            resname_line=   str(line[17:20])
            if resname_line in ligandResnames:
                continue
            else:
                out+=line

        elif line.startswith('END'):
            continue

        elif line.startswith('TER'):
            continue

        else:
            out+=line

    out+=LIGline
    out+='END\n'

    os.system('/bin/mv %s %s' %(refinedModel,refinedModel.replace('.pdb','_original_atom_order.pdb')))
    f=open(refinedModel,'w')
    f.write(out)
    f.close()



if __name__=='__main__':
    ensembleModel=sys.argv[1]
    refinedModel=sys.argv[2]
    ensembleLIGdir=get_atom_order_of_ensemble_model(ensembleModel)
    resort_ligand_atoms_in_refined_model(refinedModel,ensembleLIGdir)