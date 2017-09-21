#!/bin/bash

export XChemExplorer_DIR="."
source $XChemExplorer_DIR/setup-scripts/xce.setup-sh
#module unload ccp4
source ccp4/bin/ccp4.setup-sh

ccp4-python $XChemExplorer_DIR/XChemExplorer.py
