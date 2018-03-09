#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Dump Peptide Chains - Read a pdb file and construct the peptide chains
"""
from __future__ import print_function

import click
import os

from pyse.pdb import parse_pdb

@click.command()
@click.argument('pdbfile')
def main(pdbfile):
    pdbfn = os.path.basename(pdbfile)
    with open(pdbfile, 'r') as pdbdata:
        startsite = None
        chain = None
        peptide = ''
        endsite = None
        prevsite = None
        for pdbentry in parse_pdb(pdbdata):
            if pdbentry['atom'] != 'N':
                continue

            # Reset Case: The chain changes
            if chain != pdbentry['chain']:
                if chain:
                    print('{},{},{},{},{}'.format(pdbfn, startsite, endsite,
                                                  chain, peptide))
                startsite = pdbentry['position']
                chain = pdbentry['chain']
                peptide = ''
                endsite = None
                prevsite = None

            peptide += pdbentry['remnant']
            endsite = int(pdbentry['position'])
            if prevsite != None and endsite - 1 != prevsite:
                print('Missing sites between {} and {}'.format(endsite,
                                                               prevsite))
            prevsite = int(pdbentry['position'])

if __name__ == '__main__':
    main()
