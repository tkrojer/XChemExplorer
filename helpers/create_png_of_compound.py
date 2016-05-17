import os,sys
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


if __name__=='__main__':

    smiles =                    sys.argv[1]
    compoundID =                sys.argv[2]
    xtal =                      sys.argv[3]
    inital_model_directory =    sys.argv[4]

    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    # Draw to a file

    print 'path',os.path.join(inital_model_directory,xtal,'compound')
    os.chdir(os.path.join(inital_model_directory,xtal,'compound'))
    print 'compoundID',compoundID
    print 'mol',mol
    Draw.MolToFile(mol, "%s.png" %compoundID.replace(' ',''))
