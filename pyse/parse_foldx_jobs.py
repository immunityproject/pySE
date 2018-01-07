#!/usr/bin/env python
# # -*- coding: utf-8 -*-

"""
Parse Foldx Jobs

This tool looks for foldx jobs and generates a json file database with
the desired raw foldx data to perform further calculations.

A progress bar is output to stderr. Stdout will produce a log of any
problems found for directories that look like foldx job directories
but which are missing files or produced errors while processing.
"""
import click
import json
import os
import sys

from collections import defaultdict
from multiprocessing import Pool, cpu_count

proteins = [ 'RT', 'TAT', 'P24', 'INT', 'ZIKA_E', 'PRO', 'P17', 'REV',
             'GP120', 'NEF' ]
pdb2protein = {
    'RT_1DLO_N.pdb': 'RT',
    'p24_Hexa_3H4E_N.pdb': 'P24',
    'p17_2GOL_N.pdb': 'P17',
    'gp120_3JWD_N.pdb': 'GP120',
    'Integrase_1BIS_N.pdb': 'INT',
    'Nef_1EFN_N.pdb': 'NEF',
    'protease_Dimer_3IXO_N.pdb': 'PRO',
    'Rev_3LPH_N.pdb': 'REV',
    'tat_3MI9_N.pdb': 'TAT'
}
protein2pdb = {v: k for k,v in pdb2protein.items()}
errstream = sys.stderr

def find_foldx_jobs(directory):
    """
    A directory is considered a foldx job when it contains:
      - list.txt
      - A valid pdb name in the list.txt file per pdb2proteins

    Return a list of (protein, jobdir)
    """
    foldx_jobs = set()
    for root,dirs,files in os.walk(directory):
        if 'list.txt' not in files:
            continue
        with open(os.path.join(root, 'list.txt'), 'r') as listfile:
            pdb = listfile.read().strip()
            foldx_jobs.add((pdb2protein[pdb], root))
    return foldx_jobs

def parse_raw_buildmodel(protein, bmdir):
    pdb = protein2pdb[protein][:-len('.pdb')]
    rawmodel_fn = os.path.join(bmdir, 'Raw_BuildModel_{}.fxout'.format(pdb))
    energies = list()
    wt_energies = list()
    with open(rawmodel_fn, 'r') as rawmodel:
        for line in rawmodel.readlines()[9:]:
            if line.startswith(pdb):
                energies.append([float(e) for e in line.split('\t')[1:]])
            elif line.startswith('WT_'):
                wt_energies.append([float(e) for e in line.split('\t')[1:]])

    assert len(energies) == len(wt_energies)
    return energies, wt_energies

def mean(list):
    num = len(list)
    if num == 0:
        return float('NaN')
    return sum(list) / num

def load_energy_deltas(protein, directory):
    energies, wt_energies = parse_raw_buildmodel(protein, directory)
    return [mean([e - w for e, w in zip(es, ws)])
            for es, ws in zip(zip(*(energies)),
                              zip(*(wt_energies)))]


def read_pdb_positions(pdbfile):
    """ used by analyze_displacement.  reads x,y,z positions for each
    atom in a PDB file """
    results = []

    def stripstr(v):
        return str(v).strip()
    with open(pdbfile, 'r') as pdb:
        for line in pdb:
            if not line.startswith('ATOM'):
                continue
            # these are specified as inclusive ranges of 1-based
            # indexes to match the PDB spec
            fieldspec = [
                ('atom', 13, 16, stripstr),
                ('remnant', 18, 20, stripstr),
                ('chain', 22, 22, stripstr),
                ('position', 23, 26, int),
                ('x', 31, 38, float),
                ('y', 39, 46, float),
                ('z', 47, 54, float),
            ]
            results.append({f[0]: f[3](line[f[1] - 1:f[2]]) for f in fieldspec})
    return results

def load_displacements(protein, directory):
    pdb = protein2pdb[protein][:-len('.pdb')]
    i = 0

    displacements = defaultdict(dict)
    while True:
        pdb_basefn = '{}_1_{}.pdb'.format(pdb, i)
        out_pdb_fn = os.path.join(directory, pdb_basefn)
        wt_out_pdb_fn = os.path.join(directory, 'WT_{}'.format(pdb_basefn))

        if not os.path.exists(out_pdb_fn):
            break
        i += 1

        positions = read_pdb_positions(out_pdb_fn)
        wt_positions = read_pdb_positions(wt_out_pdb_fn)
        # Load position data into the displacements
        kfmt = '{}-{}-{}-{}'
        for p in positions:
            k = kfmt.format(p['chain'], p['position'], p['remnant'], p['atom'])
            displacements[protein][k] = dict()
            displacements[protein][k]['position'] = p
        # Generate the delta position data from the wt displacements
        for p in wt_positions:
            displacements[protein][k]['wt_position'] = p
    return displacements


def load_foldx_jobdir(foldx_jobs):
    """
    Parse the files in the provided foldx job dir, return the
    available json data.
    """
    protein, jobdir = foldx_jobs
    energies = list()
    try:
        energies = load_energy_deltas(protein, jobdir)
    except FileNotFoundError as e:
        print('Could not load energy deltas: {}'.format(e), file=errstream)

    displacements = load_displacements(protein, jobdir)

    data = defaultdict(dict)
    data[protein]['energies'] = energies
    data[protein]['displacements'] = displacements
    return data

@click.command()
@click.option('--outfile', '-o', default='EpitopeData.json',
              type=click.File('w', encoding='utf-8'),
              help='The database output file')
@click.option('--error', default='-',
              type=click.File('w', encoding='utf-8'),
              help='The database output file')
@click.argument('jobs_dir')
def main(outfile, error, jobs_dir):
    global errstream
    errstream = error
    foldx_jobs = find_foldx_jobs(jobs_dir)

    workers = Pool(cpu_count())
    results = list()
    with click.progressbar(workers.imap_unordered(load_foldx_jobdir,
                                                  foldx_jobs),
                           length=len(foldx_jobs), label='Parsing',
                           file=sys.stdout) as progbar:
        for j in progbar:
            results.append(j)

    json.dump(results, outfile)
