#/bin/csh -f

export XChemExplorer_DIR='/usr/local/scripts/tobias/XChemExplorer'

source $XChemExplorer_DIR/setup-scripts/xce.setup-sh

ccp4-python $XChemExplorer_DIR/XChemExplorer.py 

