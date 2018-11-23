[![Build Status](https://travis-ci.org/xchem/XChemExplorer.svg?branch=master)](https://travis-ci.org/xchem/XChemExplorer)
<a href="https://codeclimate.com/github/xchem/XChemExplorer/"><img src="https://codeclimate.com/github/xchem/XChemExplorer/badges/gpa.svg" /></a>
<a href="https://codeclimate.com/github/xchem/XChemExplorer/"><img src="https://codeclimate.com/github/xchem/XChemExplorer/badges/issue_count.svg" /></a>
[![HitCount](http://hits.dwyl.io/xchem/XChemExplorer.svg)](http://hits.dwyl.io/xchem/XChemExplorer)

# XChemExplorer (XCE)
<i> "The XChemExplorer graphical workflow tool for routine or large-scale protein-ligand structure determination." </i>
<p>Acta Crystallogr D Struct Biol. 2017 Mar 1;73(Pt 3):267-278. (https://doi.org/10.1107/S2059798316020234)</p> 

## Scope 

XChemExplorer (XCE) is a data-management and workflow tool to support large-scale simultaneous analysis of protein-ligand complexes during structure-based ligand discovery (SBLD). 

The user interfaces of established crystallographic software packages such as CCP4 [Winn et al. (2011), Acta Cryst. D67, 235-242] or PHENIX [Adams et al. (2010), Acta Cryst. D66, 213-221] have entrenched the paradigm that a `project' is concerned with solving one structure. This does not hold for SBLD, where many almost identical structures need to be solved and analysed quickly in one batch of work. Functionality to track progress and annotate structures is essential. 

XCE provides an intuitive graphical user interface which guides the user from data processing, initial map calculation, ligand identification and refinement up until data dissemination. It provides multiple entry points depending on the need of each project, enables batch processing of multiple data sets and records metadata, progress and annotations in an SQLite database. 

## Requirements
Operating Systems:
- Linux
- Mac OSX

<b>Windows users: </b>

Potential solutions:

a) Partition your hard drive and install a light-weight version of linux, such as Ubuntu (https://www.ubuntu.com/download/desktop)

b) Install Ubuntu (or other) on a USB drive, and boot from that: https://tutorials.ubuntu.com/tutorial/tutorial-create-a-usb-stick-on-windows#0

c) VirtualBox - emulate a linux environment on your Windows desktop (https://www.virtualbox.org)

Prerequisites:
- CCP4 version 7.0 (or higher)
- PHENIX (optional, but recommended)

<b>Please note:</b> The recommended installation process, described below, includes an install of CCP4. This is so that we know that our code works with the version of CCP4, rather than your system version, making user support much easier

## Installation
1. Clone the github repository onto your machine with:
```
git clone https://github.com/xchem/XChemExplorer
```

2. Change directory into the repository:
```
cd XChemExplorer/
```

3. Run the test_build.sh script (this is currently only included for bash, but you can view the steps within the script and modify it as you please for other shells):
```
./test_build.sh
```

4. To execute, run the XChemExplorer_local.sh script
```
./XChemExplorer_local.sh
```

(<i>We recommend you add an alias to your bash profile to do this</i>)
```
alias xce="<full_path_to_local_git_repository>/XChemExplorer_local.sh"
```

## For Diamond users: Accessing your data on your local machine

For XChemExplorer to work with your data from Diamond, you will need to mirror diamond's filesystem on your local machine. This can be done in a number of ways, but here is our reccommended route. (NB: sudo access - i.e. your system administrators - may be required)

1. Install FUSE (Filesystem in Userspace) - more info: https://github.com/libfuse/libfuse

Linux distributions:
```
sudo yum install fuse-utils sshfs
```

MacOS: https://osxfuse.github.io

2. Create a mountpoint for the diamond filesystem in your root directory (may require root/sudo access):
```
sudo mkdir /dls
```

3. Change ownership of /dls from sudo/root to yourself:
```
sudo chown <user>:<group> /dls
```

4. Now, you should be able to mount the Diamond filesytem with:
```
sshfs <fed_id>@nx.diamond.ac.uk:/dls /dls
```

5. To run XCE: change directory to your project directory, and then launch your local version (described in installation section):
```
# change to project directory
cd /dls/labxchem/data/2016/lb13385-10/processing

# launch XCE with local version script
<path-to-local-install>/XChemExplorer_local.sh

```

6. And finally, when you are done unmount Diamond's drives:
```
fusermount -u /dls
```
