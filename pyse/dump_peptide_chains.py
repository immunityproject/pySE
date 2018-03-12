#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Dump Peptide Chains - Read a pdb file and construct the peptide chains
"""
from __future__ import print_function

import click
import os

from pyse.pdb import parse_pdb,get_peptide_chains

@click.command()
@click.argument('pdbfile')
def main(pdbfile):
    pdbfn = os.path.basename(pdbfile)
    with open(pdbfile, 'r') as pdb:
        for peptide in get_peptide_chains(parse_pdb(pdb)):
            print('{},{},{},{},{}'.format(pdbfn,
                                          peptide['startsite'],
                                          peptide['endsite'],
                                          peptide['chain'],
                                          peptide['peptide']))

if __name__ == '__main__':
    main()
