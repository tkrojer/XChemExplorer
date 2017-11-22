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


