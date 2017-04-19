#/bin/bash -f

export XChemExplorer_DIR="/dls/science/groups/i04-1/software/XChemExplorer"
#source $XChemExplorer_DIR/setup-scripts/pandda.setup-sh 
source $XChemExplorer_DIR/setup-scripts/xce.setup-sh
#$CCP4/libexec/python $XChemExplorer_DIR/XChemExplorer.py
ccp4-python $XChemExplorer_DIR/XChemExplorer.py 
#/dls_sw/apps/ccp4/64/7.0/update21/ccp4-7.0/bin/ccp4-python $XChemExplorer_DIR/XChemExplorer.py


