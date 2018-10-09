#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=no-value-for-parameter

"""
Create Foldx Jobs

This tool creates foldx job directories to be parsed and processed by
scan_foldx_jobs. It takes as input a source pdb file, which is parsed
to generate all possible combinations of mutations based on the ATOM
lines read from the file.

The jobs created by this tool are used to run the foldx BuildModel
command.  The job directories contain 2 data elements:
  1. run.cfg - the run configuration, with pdb file included
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

import os
from collections import defaultdict

import click

from pyse.pdb import parse_pdb
from pyse.proteins import codes


def make_job_dir(pdbfile, jobid, line, cbeta, water_crystal, basedir='.'):
    """Generate a job directory for a work item.
    :param pdbfile: the name of the pdbfile that this work item is used with.
    :param line: the mutation(s) that this work item should perform.
    :returns: None.  Generates a folder and FoldX work files on disk. """
    dir_name = "foldxbm-{}".format(jobid)
    full_loc = os.path.join(basedir, dir_name)
    os.makedirs(full_loc)
    with open(os.path.join(full_loc, 'run.cfg'), 'w') as list_file:
        list_file.write('command=BuildModel\n')
        list_file.write('pdb=' + os.path.basename(pdbfile) + '\n')
        list_file.write('mutant-file=individual_list.txt\n')
        list_file.write('numberOfRuns=5\n')
        if cbeta:
          list_file.write('fixedCbeta=false\n')
        if water_crystal:
          list_file.write('water=-CRYSTAL\n')
    with open(os.path.join(full_loc, 'individual_list.txt'), 'w') as ilist_file:
        ilist_file.write(line + '\n')

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
@click.option('--cbeta/--no_cbeta', default=False,
              help=('The cbeta paramter is set to false'))
@click.option('--water_crystal/--no_water_crystal', default=False,
              help=('The water parameter is set to crystal'))
@click.argument('pdbfile')
def main(outdir, protein, pdbfile,cbeta,water_crystal):
    """CLI for creating jobs"""
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
                    make_job_dir(pdbfile, jobid, mutlist, cbeta, water_crystal, outdir)


if __name__ == '__main__':
    main()
