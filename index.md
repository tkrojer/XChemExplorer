# [](#header-1)Introduction

XChemExplorer (XCE) is a data management and workflow tool for the parallel determination of protein-ligand structures. It was initially written to support crystallographic fragment screening at [beamline I04-1 at the Diamond Light Source](http://www.diamond.ac.uk/Beamlines/Mx/Fragment-Screening.html), but it is a generic program that can be used to facilitate any structure-based drug design project.


# [](#header-1)Download
A pre-release is currently available for download:

[v1.0-beta.3.4](https://github.com/tkrojer/XChemExplorer/archive/v1.0-beta.3.4.tar.gz)



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
gunzip XChemExplorer-1.0-beta.3.4.tar.gz
tar –xvf XChemExplorer-1.0-beta.3.4.tar
```

This will create a new directory, i.e. from now on your XChemExplorer directory. Change into this directory, e.g.:

```shell
cd XChemExplorer-1.0-beta.3.4
```

The contents of the directory should look something like this when you type ‘ls –l’:
```
-rwxr-xr-x   1 tobiaskrojer  staff   185B  2 Mar 09:53 XChemExplorer.csh
-rwxr-xr-x   1 tobiaskrojer  staff   309K  2 Mar 09:53 XChemExplorer.py
-rwxr-xr-x   1 tobiaskrojer  staff   186B  2 Mar 09:53 XChemExplorer.sh
drwxr-xr-x  11 tobiaskrojer  staff   374B  2 Mar 09:53 helpers
drwxr-xr-x  12 tobiaskrojer  staff   408B 26 Jan 10:53 image
drwxr-xr-x  34 tobiaskrojer  staff   1.1K  9 Mar 14:41 lib
drwxr-xr-x   4 tobiaskrojer  staff   136B  2 Mar 09:53 setup-scripts
drwxr-xr-x   7 tobiaskrojer  staff   238B  9 Mar 14:41 web
```

The only thing left to do is to edit the XChemExplorer.sh or XChemExplorer.csh file, depending on which shell you are using. Open XChemExplorer.sh for bash shells or XChemExplorer.csh for C-shells with your editor of choice and edit the line

```shell
export XChemExplorer_DIR='/usr/local/scripts/tobias/XChemExplorer'
```

to where you XChemExplorer is installed. In our example this would be 

```shell
export XChemExplorer_DIR='/home/tkrojer/software/XChemExplorer-1.0-beta.3.4’
```

That’s it! 




# Getting Started

You can now run XCE by typing

```shell
/home/tkrojer/software/XChemExplorer-1.0-beta.3.4/XChemExplorer.sh
```

It may however be easier if you insert an alias into your .bashrc or .cshrc file:

```shell
alias XChemExplorer='/home/tkrojer/software/XChemExplorer-1.0-beta.3.4/XChemExplorer.sh
```

# Manual

[Download](ftp://ftp.sgc.ox.ac.uk/pub/tkrojer/XChemExplorer/Export_HTML_summary.pdf)

# Reference

Krojer, T., Talon, R., Pearce, N., Collins, P., Douangamath, A., Brandao-Neto, J., Dias, A., Marsden, B., and von Delft, F. (2017). The XChemExplorer graphical workflow tool for routine or large-scale protein–ligand structure determination. _Acta Cryst D_ 73, 267–278.


# Contact

Tobias Krojer
tobias.krojer@sgc.ox.ac.uk


# Links

* [SGC Oxford](http://www.thesgc.org/scientists/groups/oxford)

* [XChem](http://www.diamond.ac.uk/Beamlines/Mx/Fragment-Screening.html)

* [PanDDA](https://pandda.bitbucket.io)
