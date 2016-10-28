#! /bin/sh -e

XChemExplorer_DIR=/tmp/XChemExplorer

if [ ! -x $CCP4/libexec/python ];then
  echo "You must first load your CCP4 (at least version 7.0) environment."
  exit
fi

. $XChemExplorer_DIR/setup-scripts/xce.setup-sh

if [ -r ${XChemExplorer_DIR}/XChemExplorer.py ];then
  export XChemExplorer_DIR 
  $CCP4/libexec/python $XChemExplorer_DIR/XChemExplorer.py
fi
