#!/bin/bash

export XChemExplorer_DIR="/dls/science/groups/i04-1/elliot-dev/XChemExplorer"
source $XChemExplorer_DIR/setup-scripts/xce.setup-sh
module unload ccp4
source /dls/science/groups/i04-1/software/pandda-update/ccp4/ccp4-7.0/bin/ccp4.setup-sh

ccp4-python $XChemExplorer_DIR/XChemExplorer.py
