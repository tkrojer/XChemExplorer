#!/bin/bash

if [[ "$OSTYPE" == "linux-gnu" ]]; then
        wget http://devtools.fg.oisin.rc-harwell.ac.uk/nightly/ccp4-linux64-latest.tar.bz2
	bunzip2 ccp4-linux64-latest.tar.bz2
	mkdir ./ccp4
	tar -xf ccp4-linux64-latest.tar -C ./ccp4 --strip-components=1
elif [[ "$OSTYPE" == "darwin"* ]]; then
        wget http://devtools.fg.oisin.rc-harwell.ac.uk/nightly/ccp4-osx-clang-latest.tar.gz
	gunzip ccp4-osx-clang-latest.tar.gz
	mkdir ./ccp4
	tar -xf ccp4-osx-clang-latest.tar -C ./ccp4 --strip-components=1
fi
cd ccp4
yes y| ./BINARY.setup > /dev/null 2>&1
source bin/ccp4.setup-sh
yes y | ccp4-python -m pip uninstall panddas
ccp4-python -m pip install panddas
#git clone https://www.github.com/xchem/XChemExplorer
