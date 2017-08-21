#/bin/csh -f

export XChemExplorer_DIR='/dls/science/groups/i04-1/software/xce-elliot-dev/xce/XChemExplorer'

source $XChemExplorer_DIR/setup-scripts/xce.setup-sh

ccp4-python $XChemExplorer_DIR/XChemExplorer.py 

