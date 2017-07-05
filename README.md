# pySE - Python Structural Entropy

A collection of python tools used to explore structural entropy, a novel approach to identifying epitope vaccine targets.

This repository contains everything required for data preparation, parallel protein folding orchestration, and result analysis.

Questions/inquiries please visit: https://www.immunityproject.org/


## Install

Quickstart (Ubuntu):

We recommend using [virtualenv](https://github.com/pypa/virtualenv) to help separate project dependencies.


```
sudo apt-get install python python-pip
git clone https://github.com/immunityproject/pySE
pip install -r pySE/requirements.txt
```



## Directory Structure


* **/data** - processed data, such as the output from FoldX PDB analysis via *SequenceOnly* command
*  **/pdbs** - Raw PDB files.  These are the PDBs referenced in our papers.  It is usually required to perform some data cleaning and repair for each PDB if working with new material. 
*  **/src** - Python tasks for working with PDBs, generating FoldX jobs, and calculating Structural Entropy
