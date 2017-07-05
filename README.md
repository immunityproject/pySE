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


## /src scripts
* **best_match.py** - Given an epitope amino sequence or file containing many amino sequences, search all known proteins for close matches for this sequence.  Output chain specific offsets for epitope start and finish and a confidence value.
* **clean_lockfiles.py** - Search all FoldX task directories and clear lockfiles, resetting these tasks for work.
* **compute_shannon_entropy.py** -  Computes Structural Entropy for a given epitope and protein, once the FoldX substitution tasks are finished.
* **convert_foldx_output_directories_to_csv.py** - Collects energy change results from all FoldX directories into a CSV file.
* **generate_buildmodel.py** - Given a protein name, generate a set of FoldX compatible run directories which replace the current WT with each other possible amino for all relevant chains.
* **generate_mutant_string.py** - Deprecated script, replaced with buildmodel.  Performed replacements with FoldX's positionscan command, and did direct parsing of PDBs. 
* **proteins.py** - Settings and configuration variables for other scripts.
* **scan_foldx_jobs.py** - Examine running FoldX job directories to determine status.
* **scan_to_attributes.py** -


