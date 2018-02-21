# last edited: 12/12/2016, 15:00

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem,Draw
from rdkit.Chem.Draw import IPythonConsole



#class generate_all_stereoisomers:

    # https://github.com/rdkit/rdkit/issues/626


mol = Chem.MolFromSmiles('FC(Cl)C=CC(C(O)N)C=CC(Cl)Br')

def spam(n):
    out=[]
    for perm in getPerms(n):
        elem = [ int(i) for i in list(perm) ]
        out.append(elem)
    return out

def getPerms(n):
    from itertools import permutations
    for i in getCandidates(n):
        for perm in set(permutations(i)):
            yield ''.join(perm)

def getCandidates(n):
    for i in range(0, n+1):
        res = "1" * i + "0" * (n - i)
        yield res

def GetStereoIsomers(mol):
    from rdkit import Chem
    from copy import copy
    out = []

    chiralCentres = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

    #return the molecule object when no chiral centres where identified
    if chiralCentres == []:
        return [mol]

    #All bit permutations with number of bits equals number of chiralCentres
    elements = spam(len(chiralCentres))

    for isoId,element in enumerate(elements):
        for centreId,i in enumerate(element):
            atomId = chiralCentres[centreId][0]
            if i == 0:
                mol.GetAtomWithIdx(atomId).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
            elif i == 1:
                mol.GetAtomWithIdx(atomId).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
        outmol = copy(mol        )
        out.append(outmol)
        print Chem.MolToSmiles(mol,isomericSmiles=True)
    return out

Draw.MolsToGridImage(GetStereoIsomers(mol))
