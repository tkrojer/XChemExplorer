#!/bin/bash
#module purge
export XChemExplorer_DIR="/home/uzw12877/XCE_currentdev/XChemExplorer"
source $XChemExplorer_DIR/setup-scripts/xce.setup-sh
module unload ccp4
source /dls/science/groups/i04-1/software/pandda-update/ccp4-linux64-2016-10-22-0115/bin/ccp4.setup-sh
ccp4-python $XChemExplorer_DIR/XChemExplorer.py 



