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

If you encounter problems compiling the `bigfloat` package, make sure
you have installed your platform's MPFR development package.
```
sudo apt-get install libmpfr-dev
```

** todo ** how to download and install data files


## Usage

```
 compute_shannon_entropy.py [-h] [--debug] [-q] [-p PROTEIN] [-s SITES]
                            [-e EPITOPE] [-d] [-b {absolute,raw}]
                            [--exclude-wt] [-c CSV_FILENAME]
                            [--no-header] [--tab-delimited]
                            [--scan SCAN] [--energy-map]
                            jobs_dir

positional arguments:
  jobs_dir

optional arguments:
  -h, --help            show this help message and exit
  --debug
  -q                    warnings are not printed
  -p PROTEIN, --protein PROTEIN
                        must be used in conjunction with --sites
  -s SITES, --sites SITES
                        the site range to analyze, using a dash range:
                        "115-121". must also specify --protein
  -e EPITOPE, --epitope EPITOPE
                        the case-insensitive epitope code to compute "KF11"
                        Automatically sets --protein and --site
  -d, --displacement    use displacement calculation instead of energies.
                        requires a flag which sets the protein.
  -b {absolute,raw}, --baseline {absolute,raw}
                        how to preprocess energies before entropy calculation
  --exclude-wt          if set, the Boltzmann distribution does not include
                        the WT energy difference (always 0.0)
  -c CSV_FILENAME, --csv CSV_FILENAME
                        base pathname (without extension) for a csv dump of
                        the data
  --no-header           do not put a header row in the CSV files
  --tab-delimited       use tab as a field delimiter in the CSV file
  --scan SCAN           if set, the entropies of all ranges of the specified
                        width are calculated. If --protein is not specified
                        then --all is implied
  --energy-map          compute the matrix of average energy deltas for all
                        possible amino acid substitutions
```

The `jobs_dir` is the top-level directory containing all the FoldX jobs
which have been previously processed.  See
[HIV-SE](https://github.com/immunityproject/HIV-SE) for HIV data.

The normal mode of operation computes the Structural Entropies for all
epitopes on all proteins.  Using `--displacement` is an experimental
mode which computes the sum of squared displacements due to amino acid
substitutions, allowing analysis of the physical change in the protein
due to mutations from WT.  Using `--scan` computes the SE for all
segments of a given length in a protein, allowing for analysis of areas
of the protein which have low SE's regardless of any assocation with an
epitope.  The output of `--scan 1` may be further processed using
`scan_to_attributes.py` to create a visualization of a protein's SE with
Swiss PDB Viewer.

Output is normally sent to stdout, but a csv filename may be supplied
with `--csv <basename>`.

#### Other Options

The `compute_shannon_entropy.py` contains a few other settings which
were used during exploration but did not get command-line switches.  See
the source file for adjusting energies based on the rotamer symmetry,
for changing units from shannons to nats, for adjusting energies with
the Boltzmann Constant and ambient temperature, etc.


## Directory Structure


* **data/** - processed data, such as the output from FoldX PDB analysis via *SequenceOnly* command
*  **pdbs/** - Raw PDB files.  These are the PDBs referenced in our papers.  It is usually required to perform some data cleaning and repair for each PDB if working with new material.
*  **src/** - Python tasks for working with PDBs, generating FoldX jobs, and calculating Structural Entropy


## src/ scripts
* **best_match.py** - Given an epitope amino sequence or file containing many amino sequences, search all known proteins for close matches for this sequence.  Output chain specific offsets for epitope start and finish and a confidence value.
* **clean_lockfiles.py** - Search all FoldX task directories and clear lockfiles, resetting these tasks for work.
* **compute_shannon_entropy.py** -  Computes Structural Entropy for a given epitope and protein, once the FoldX substitution tasks are finished.
* **convert_foldx_output_directories_to_csv.py** - Collects energy change results from all FoldX directories into a CSV file.
* **generate_buildmodel.py** - Given a protein name, generate a set of FoldX compatible run directories which replace the current WT with each other possible amino for all relevant chains.
* **generate_mutant_string.py** - Deprecated script, replaced with buildmodel.  Performed replacements with FoldX's positionscan command, and did direct parsing of PDBs. 
* **proteins.py** - Settings and configuration variables for other scripts.
* **scan_foldx_jobs.py** - Examine running FoldX job directories to determine status.
* **scan_to_attributes.py** -

