#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Create Foldx Jobs

This tool creates foldx job directories to be parsed and processed by
scan_foldx_jobs. It takes as input a source pdb file, which is parsed
to generate all possible combinations of mutations based on the ATOM
lines read from the file.

The jobs created by this tool are used to run the foldx BuildModel
command.  The job directories contain 2 data elements:
  1. list.txt - the name of the source pdb file. This is usually a repaired
                pdb file.
  2. individual_list.txt - wild type amino (starting amino letter), the chain
                           (A, B, C, ...), the residue number (site position),
                           and the mutant residue (target mutation).

See the BuildModel command in the foldx documentation for more information.

We generate a separate job directory for every site position and
target mutation, grouped by the chains with matching wild types in
that position. For example, if site 2 of a given pdb has the following
wild types and chains:
  A. LYS (K)
  B. LYS (K)
  C. GLU (E)
  D. HIS (H)
  E. GLU (E)

This will produce 3 sets: {A, B}, {C, E}, and {D}. Each set will produce a
mutation against each target, including itself.

The original version of this tool can be found in generate_buildmodel.py

"""
from __future__ import print_function

import click
import os

from collections import defaultdict

from pyse.pdb import parse_pdb
from pyse.proteins import codes


def makeJobDir(pdbfile, jobid, line, basedir='.'):
    """Generate a job directory for a work item.
    :param pdbfile: the name of the pdbfile that this work item is used with.
    :param line: the mutation(s) that this work item should perform.
    :returns: None.  Generates a folder and FoldX work files on disk. """
    dir_name = "foldxbm-{}".format(jobid)
    full_loc = os.path.join(basedir, dir_name)
    os.makedirs(full_loc)
    with open(os.path.join(full_loc, 'list.txt'), 'w') as il:
        il.write(os.path.basename(pdbfile) + '\n')
    with open(os.path.join(full_loc, 'individual_list.txt'), 'w') as l:
        l.write(line + '\n')

def generate_chaingroups(pdbfn):
    """This parses the pdb then generates the following structure for
    generating build job folders:
      site: {wildtype: [groups]}
    """
    chaingroups = defaultdict(dict)
    with open(pdbfn, 'r') as pdbfile:
        for pdbentry in parse_pdb(pdbfile):
            wildtype = pdbentry['remnant']
            chain = pdbentry['chain']
            site = pdbentry['position']
            chains = chaingroups[site].get(wildtype, set())
            chains.add(chain)
            chaingroups[site][wildtype] = chains
    return chaingroups

@click.command()
@click.option('--outdir', '-o', default='.',
              help='The output job directory')
@click.option('--protein', '-p', default=None,
              help=('The name of the protein (jobs are put in '
                    'jobdir/proteinname/jobname)'))
@click.argument('pdbfile')
def main(outdir, protein, pdbfile):
    outdir = os.path.join(outdir, protein)
    print('Reading {} and putting job directories in {}'.format(pdbfile,
                                                                outdir))
    chaingroups = generate_chaingroups(pdbfile)
    for site, chaingroup in chaingroups.items():
        for wildtype, chains in chaingroup.items():
            wtamino = codes[wildtype]
            for chain in chains:
                for mut in codes.values():
                    mutlist = '{}{}{}{};'.format(wtamino, chain, site, mut)
                    jobid = '{}-{}-{}-{}'.format(site, chain, wtamino, mut)
                    makeJobDir(pdbfile, jobid, mutlist, outdir)


if __name__ == '__main__':
    main()
