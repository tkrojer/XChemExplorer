import os,sys
#sys.path.append(os.path.join(os.getenv('CCP4'),'lib','python2.7','site-packages'))
sys.path.append(os.path.join(os.getenv('PHENIX'),'base','Python.framework','Versions','2.7','lib','python2.7','site-packages'))

import PIL
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

def MakeImage(CompoundID,smiles):
	mol = Chem.MolFromSmiles(smiles)
	AllChem.Compute2DCoords(mol)
	# Draw to a file
	Draw.MolToFile(mol, "%s.png" %CompoundID)
	# Draw to image to be used in Python
	py_image = Draw.MolToImage(mol)
	print 'making %s.png' %CompoundID


if __name__=='__main__':
#	compoundID=sys.argv[1]
#	smiles=sys.argv[2]
#	folder=sys.argv[3]

	MakeImage('XXX','OCCO')
