# [](#header-1)Introduction

XChemExplorer (XCE) is a data management and workflow tool for the parallel determination of protein-ligand structures. It was initially written to support crystallographic fragment screening at [beamline I04-1 at the Diamond Light Source](http://www.diamond.ac.uk/Beamlines/Mx/Fragment-Screening.html), but it is a generic program that can be used to facilitate any structure-based drug design project.


# [](#header-1)Download

[v1.3](https://github.com/xchem/XChemExplorer/archive/v1.3.tar.gz)


# [](#header-1)Mailing List

For questions, suggestions, hints, bug reports, please sign up to XCHEMBB@jiscmail.ac.uk


# [](#header-1)Operating System
* LINUX
* Mac OS X


# [](#header-1)Prerequisites:
* CCP4 version 7.0 (or higher)
* PHENIX (optional, but recommended)


# [](#header-1)Installation

Put the gzipped tar archive to wherever you want XCE to be installed. In case you have no root privileges, put it somewhere into your home directory, e.g.:

```
/home/tkrojer/software
```

Then change to the respective directory and unpack the archive, e.g.:

```shell
cd /home/tkrojer/software
gunzip XChemExplorer-1.2.tar.gz
tar –xvf XChemExplorer-1.2.tar
```

This will create a new directory, i.e. from now on your XChemExplorer directory. Change into this directory, e.g.:

```shell
cd XChemExplorer-1.2
```

The contents of the directory should look something like this when you type ‘ls –l’:
```
-rw-r--r--@  1 tobiaskrojer  staff     182 26 Jan 10:03 Dockerfile
-rw-r--r--@  1 tobiaskrojer  staff    2832 26 Jan 10:03 README.md
lrwxr-xr-x@  1 tobiaskrojer  staff      20 26 Jan 10:03 XChemExplorer -> XChemExplorer_dmd.sh
-rwxr-xr-x@  1 tobiaskrojer  staff  222991 26 Jan 10:03 XChemExplorer.py
-rwxr-xr-x@  1 tobiaskrojer  staff     316 26 Jan 10:03 XChemExplorer_dmd.sh
-rwxr-xr-x@  1 tobiaskrojer  staff     269 26 Jan 10:03 XChemExplorer_local.sh
-rw-r--r--@  1 tobiaskrojer  staff     465 26 Jan 10:03 compile_test.py
drwxr-xr-x@ 13 tobiaskrojer  staff     442 26 Jan 10:03 gui_scripts
drwxr-xr-x@ 14 tobiaskrojer  staff     476 26 Jan 10:03 helpers
drwxr-xr-x@ 10 tobiaskrojer  staff     340 26 Jan 10:03 icons
drwxr-xr-x@ 11 tobiaskrojer  staff     374 26 Jan 10:03 image
drwxr-xr-x@ 21 tobiaskrojer  staff     714 26 Jan 10:03 lib
-rwxr-xr-x@  1 tobiaskrojer  staff      43 26 Jan 10:03 run_tests
drwxr-xr-x@  5 tobiaskrojer  staff     170 26 Jan 10:03 setup-scripts
-rwxr-xr-x@  1 tobiaskrojer  staff     553 26 Jan 10:03 setupssh.sh
-rwxr-xr-x@  1 tobiaskrojer  staff     809 26 Jan 10:03 test_build.sh
drwxr-xr-x@  6 tobiaskrojer  staff     204 26 Jan 10:03 web
```

The only thing left to do is to edit the XChemExplorer_dmd.sh file. After you open XChemExplorer_dmd.sh with your editor of choice, the file will look like this:

```shell
#!/bin/bash

export XChemExplorer_DIR="/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer"
source $XChemExplorer_DIR/setup-scripts/xce.setup-sh
module unload ccp4
source /dls/science/groups/i04-1/software/pandda-update/ccp4/ccp4-7.0/bin/ccp4.setup-sh

ccp4-python $XChemExplorer_DIR/XChemExplorer.py
```

Change the first line (export XChemExplorer_Dir=...) to wherever you have the program installed (in the example, XChemExplorer is installed in /Users/tobiaskrojer/Downloads/XChemExplorer-1.2). And remove lines 3 and 4! After these edits, the file should look something like this:

```shell
#!/bin/bash

export XChemExplorer_DIR="/Users/tobiaskrojer/Downloads/XChemExplorer-1.1"
source $XChemExplorer_DIR/setup-scripts/xce.setup-sh

ccp4-python $XChemExplorer_DIR/XChemExplorer.py
```

That’s it! 




# Getting Started

You can now run XCE by typing

```shell
/Users/tobiaskrojer/Downloads/XChemExplorer-1.2/XChemExplorer
```

It may however be easier if you insert an alias into your .bashrc file:

```shell
alias XChemExplorer='/Users/tobiaskrojer/Downloads/XChemExplorer-1.1/XChemExplorer'
```

# Manual - under development -

* [XChemExplorer v1.3](https://github.com/tkrojer/XChemExplorer/blob/gh-pages/XCE_manual_2019-01-11.pdf)



# Reference

Krojer, T., Talon, R., Pearce, N., Collins, P., Douangamath, A., Brandao-Neto, J., Dias, A., Marsden, B., and von Delft, F. (2017). The XChemExplorer graphical workflow tool for routine or large-scale protein–ligand structure determination. _Acta Cryst D_ 73, 267–278.


# Contact

Tobias Krojer
tobias.krojer@sgc.ox.ac.uk

Rachael Skyner
rachael.skyner@diamond.ac.uk

# Feature Requests

*** under construction ***


# Version History

* v1.3 - 20/12/2018
* v1.2.1 - 20/06/2018
* v1.2 - 18/06/2018
* v1.1 - 31/01/2018
* v1.0 - 10/08/2017


# Links

* [SGC Oxford](http://www.thesgc.org/scientists/groups/oxford)

* [XChem](http://www.diamond.ac.uk/Beamlines/Mx/Fragment-Screening.html)

* [PanDDA](https://pandda.bitbucket.io)


# Resources

* <a href="https://github.com/tkrojer/XChemExplorer/blob/gh-pages/example.deposit">example.deposit file</a>
* <a href="https://github.com/tkrojer/XChemExplorer/blob/gh-pages/data_template.cif">data_template.cif file</a>
* <a href="https://github.com/tkrojer/XChemExplorer/blob/gh-pages/index.txt">index.txt file</a>
* <a href="https://github.com/tkrojer/XChemExplorer/blob/gh-pages/NUDT7A.tar.bz2">NUDT7A.tar.bz2 file</a>

