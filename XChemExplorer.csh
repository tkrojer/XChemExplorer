#/bin/csh -f

setenv XChemExplorer_DIR /usr/local/scripts/tobias/XChemExplorer
source $XChemExplorer_DIR/setup-scripts/xce.setup-csh 
$CCP4/libexec/python $XChemExplorer_DIR/XChemExplorer.py $1

