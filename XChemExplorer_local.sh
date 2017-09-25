#!/bin/bash

export XChemExplorer_DIR="~/mounting_dir"
source $XChemExplorer_DIR/setup-scripts/xce.setup-sh
#module unload ccp4
source /Applications/ccp4-7.0/setup-scripts/ccp4.setup-sh

ccp4-python $XChemExplorer_DIR/XChemExplorer.py
