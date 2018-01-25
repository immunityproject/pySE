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
from __future__ import print_function

import click
import csv
import functools
import gzip
import gc
import json
import os
import re
import sys

from collections import defaultdict
from multiprocessing import Pool, cpu_count

# Define and error printing function in lieu of logging api
eprint = functools.partial(print, file=sys.stderr, flush=True)

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

def find_foldx_jobs(directory):
    """
    A directory is considered a foldx job when it contains:
      - list.txt
      - individual_list.txt
      - A valid pdb name in the list.txt file per pdb2proteins
      - A parseable wt, site, mutation in individual_list.txt

    Return a list of (protein, wt, site, mutation, jobdir). Folders
    missing list.txt or individual_list.txt files are skipped. We make
    no guarantees about validity of the jobs. Checking and reporting
    on that should go elsewhere.
    """
    individual_list_parser = re.compile(r'(\w)(\w)(\d+)(\w)')
    foldx_jobs = set()
    for root,dirs,files in os.walk(directory):
        if 'list.txt' not in files or 'individual_list.txt' not in files:
            continue

        protein = None
        wt = None
        site = None
        mutation = None
        jobdir = None
        with open(os.path.join(root, 'list.txt'), 'r') as listfile:
            pdb = listfile.read().strip()
            protein = pdb2protein[pdb]
            jobdir = root

        with open(os.path.join(root, 'individual_list.txt')) as indfile:
            # This will take a line like: EA28E,EB28E,EC28E,ED28E,EE28E,EF28E;
            # It reads the first element: EA28E
            # Then EA28E -> wt = E, site = 28, and mutation = E.
            # NOTE: A is the chain, which shows up again in the .pdb outputs
            individual_list = indfile.read().strip().split(',')
            parsed_ind_list = individual_list_parser.match(individual_list[0])
            wt, site, mutation = parsed_ind_list.group(1, 3, 4)

        foldx_job = (protein, wt, site, mutation, jobdir)
        foldx_jobs.add(foldx_job)
    return foldx_jobs

def parse_raw_buildmodel(pdb, rawmodel):
    """Parse a raw buildmodel file, returning a tuple of the energies
    and wt_energies. We assert that both lists are non-empty and that
    they are the same length."""
    def parse_energy_line(el):
        energy_keys = [ "total energy", "Backbone Hbond", "Sidechain Hbond",
                        "Van der Waals", "Electrostatics",
                        "Solvation Polar", "Solvation Hydrophobic",
                        "Van der Waals clashes", "entropy sidechain",
                        "entropy mainchain", "sloop_entropy",
                        "mloop_entropy", "cis_bond", "torsional clash",
                        "backbone clash", "DNA stack clash",
                        "helix dipole", "water bridge", "disulfide",
                        "electrostatic kon", "partial covalent bonds",
                        "energy Ionisation"]
        # This number is typically missing, so we don't include it:
        #   , "Entropy Complex" ]
        energies = [float(e) for e in line.split('\t')[1:]]
        return { energy_keys[i]: energies[i] for i in range(len(energy_keys)) }

    energies = list()
    wt_energies = list()
    for line in rawmodel.readlines()[9:]:
        # if line.startswith('Pdb') and energy_keys is None:
        #     energy_keys = [k for k in line.split('\t')[1:]]
        if line.startswith(pdb):
            energies.append(parse_energy_line(line))
        elif line.startswith('WT_'):
            wt_energies.append(parse_energy_line(line))

    assert len(energies) != 0
    assert len(wt_energies) != 0
    assert len(energies) == len(wt_energies)

    return energies, wt_energies

def parse_pdb(pdb):
    """ Read the x,y,z positions for each atom in a PDB file """
    positions = []

    def stripstr(v):
        return str(v).strip()
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
        positions.append({f[0]: f[3](line[f[1] - 1:f[2]]) for f in fieldspec})
    return positions

def calculate_energy_deltas(energies, wt_energies):
    """Sum the columns in the energies and wt_energies lists and
    return a collapsed list of the averages

    This code replaces the following zip magic that does the same thing:

      return [mean([e - w for e, w in zip(es, ws)])
              for es, ws in zip(zip(*(energies)),
                                zip(*(wt_energies)))]
    """
    energy_deltas = { k: 0.0 for k in energies[0].keys() }
    for k in energy_deltas.keys():
        for i in range(len(energies)):
            energy_deltas[k] += (energies[i][k] - wt_energies[i][k])
        energy_deltas[k] = energy_deltas[k]/len(energies)

    return energy_deltas

def get_displacement_files(pdb, directory):
    """ Create the full path filename pairs for the 5 Wild Type (WT)
    and output pdb files """
    disp_files = list()
    for i in range(5):
        pdb_basefn = '{}_1_{}.pdb'.format(pdb, i)
        out_pdb_fn = os.path.join(directory, pdb_basefn)
        wt_out_pdb_fn = os.path.join(directory, 'WT_{}'.format(pdb_basefn))
        disp_files.append((out_pdb_fn, wt_out_pdb_fn))
    return disp_files

def load_displacements(positions, displacements, disp_type):
    """load the displacements into the given dictionary"""
    kfmt = '{}-{}-{}-{}'
    for p in positions:
        k = kfmt.format(p['chain'], p['position'], p['remnant'], p['atom'])
        if k not in displacements:
            displacements[k] = {
                'atom': p['atom'],
                'remnant': p['remnant'],
                'chain': p['chain'],
                'position': p['position']
            }
        displacements[k][disp_type] = {
            'x': p['x'],
            'y': p['y'],
            'z': p['z']
        }

def check_buildmodel(buildmodel):
    """Reads the buildmodel file and raises exceptions if there are
    PROBLEM lines"""
    problems = list()
    for line in buildmodel:
        if "PROBLEM" in line:
            problems.append(line.rstrip())

    if len(problems) != 0:
        raise Exception('{}'.format(problems))

def calculate_displacement_deltas(displacements):
    """If position and wt_position exist, add a 'delta' field with
    the differences in the x,y,z values"""
    for k, v in displacements.items():
        if 'positions' not in v or 'wt_positions' not in v:
            continue
        try:
            displacements[k]['deltas'] = {
                'x': (v['positions']['x'] - v['wt_positions']['x']),
                'y': (v['positions']['y'] - v['wt_positions']['y']),
                'z': (v['positions']['z'] - v['wt_positions']['z'])
            }
        except Exception as e:
            eprint(v)
            raise e

def load_foldx_job(foldx_job):
    """ Parse the files in the provided foldx job, return the
    available json data.  """
    protein, wt, site, mutation, jobdir = foldx_job
    job_bn = os.path.basename(jobdir) # Use basename for reporting

    jobid = "{},{},{},{},{}".format(protein, wt, site, mutation, job_bn)
    pdb = protein2pdb[protein][:-len('.pdb')]

    # First check the buildmodel file. If there are errors present,
    # skip processing and return an empty result. Print the error to
    # the screen.
    try:
        buildmodel_fn = os.path.join(jobdir,
                                     'BuildModel_{}.fxout'.format(pdb))
        with open(buildmodel_fn, 'r') as buildmodel:
            check_buildmodel(buildmodel)
    except Exception as e:
        eprint('{},Detected FoldX Errors,{}'.format(jobid, e))
        return defaultdict(dict)

    # Calculat energy deltas
    energies = list()
    wt_energies = list()
    try:
        rawmodel_fn = os.path.join(jobdir,
                                   'Raw_BuildModel_{}.fxout'.format(pdb))
        with open(rawmodel_fn) as rawmodel:
            energies, wt_energies = parse_raw_buildmodel(pdb, rawmodel)
            energy_deltas = calculate_energy_deltas(energies, wt_energies)
    except FileNotFoundError as e:
        eprint('{},Could not load energy deltas,{}'.format(jobid, e))

    # Calculate displacements
    disp_files = get_displacement_files(pdb, jobdir)
    displacements = defaultdict(dict)
    for out_pdb_fn, wt_out_pdb_fn in disp_files:
        try:
            with open(out_pdb_fn, 'r') as out_pdb:
                load_displacements(parse_pdb(out_pdb), displacements,
                                   'positions')
            with open(wt_out_pdb_fn, 'r') as out_pdb:
                load_displacements(parse_pdb(out_pdb), displacements,
                                   'wt_positions')
            calculate_displacement_deltas(displacements)
        except Exception as e:
            eprint('{},Could not load displacements,{}'.format(jobid, e))

    data = dict()
    data[jobid] = {
        'protein': protein,
        'wt': wt,
        'site': site,
        'mutation': mutation,
        'jobdir': job_bn,
        'energy_deltas': energy_deltas,
        'displacements': displacements
    }
    return data

def parse_foldx_jobs(outfile, epitopes, jobs_dir):
    foldx_jobs = list(find_foldx_jobs(jobs_dir))
    emptyepitope = { k: '' for k in [ "peptide", "protein", "subprotein",
                                      "start", "end", "epitope", "subtype" ]}

    workers = Pool(cpu_count())
    with click.progressbar(workers.imap_unordered(load_foldx_job,
                                                  foldx_jobs),
                           length=len(foldx_jobs), label='Parsing',
                           file=sys.stdout) as progbar:
        for j in progbar:
            if j:
                for key in j.keys():
                    kps = key.split(',')[0:4]
                    epitopekey = '{},{},{}'.format(kps[0], kps[1], kps[2])
                    for epitope in epitopes.get(epitopekey, [emptyepitope]):
                        for ek,ev in epitope.items():
                            if ek == 'protein':
                                ek = 'polyprotein'
                            j[key][ek] = ev
                        outfile.write((json.dumps(j) + '\n').encode())
                outfile.flush()
                j.clear()
                gc.collect()

def load_epitopes(epitopefile):
    """Parse a csv of epitopes of the form:
        peptide,polyprotein,subprotein,start,end,epitope,...
    and return a dict representation
    """
    epitopes = []
    header = []
    epreader = csv.reader(epitopefile)
    for row in epreader:
        if row[0] == 'peptide':
            header = row
            continue
        epitope = { header[i]: row[i] for i in range(len(row)) }
        epitopes.append(epitope)

    return epitopes

def resolve_subprotein_site(site, pidx, length, subprotein):
    """subproteins come in the format subprotein(start-end) and
    subprotein1(start)-subprotein2(end), the latter requires us to
    determine if we are in the first or second subprotein space."""
    parts = subprotein.split('-')
    eprint(parts)
    if ')' not in parts[0]:
        p,s = parts[0].split('(')[:2]
        return p,(int(s) + pidx)
    else:
        p1, start = parts[0].split('(')[:2]
        p2, end = parts[1].split('(')[:2]
        start = int(start[:-1])
        end = int(end[:-1])
        second_site = pidx - (length - end)
        if second_site > 0: # Second subprotein
            return p2,second_site
        else:
            return p1,(pidx+start)

def add_epmap_dict(epmap, k, e):
    for ep in epmap[k]:
        if e['peptide'] == ep['peptide']:
            return
    epmap[k].append(e)

def generate_epitope_map(epitopedb):
    """For each epitope in the epitopes list, map it to a job key

    Limit the key outputs to the proteins and subproteins in the global
    proteins list.
    """
    global proteins

    # Map matched epitopes to job prefixes
    kfmt = '{},{},{}'
    epmap = defaultdict(list)
    for e in epitopedb:
        start = int(e['start'])
        end = int(e['end'])

        for site in range(start, end):
            pidx = site-start
            wt = e['peptide'][pidx]
            protein = e['protein']
            subprotein = e['subprotein']
            subprotein.replace('POL-TF', 'POL_TF')
            for p in proteins:
                if p in protein:
                    k = kfmt.format(p,wt,site)
                    add_epmap_dict(epmap, k, e)
                if p in subprotein:
                    # NOTE: this adds keys outside the protein space,
                    # which is why we keep the full peptide and
                    # subproteins when adding proteins
                    length = start-end
                    sp,subsite = resolve_subprotein_site(site, pidx,
                                                         length,
                                                         subprotein)
                    eprint(sp,subsite)
                    if sp.startswith(p):
                        k = kfmt.format(p,wt,subsite)
                        add_epmap_dict(epmap, k, e)

    return epmap

@click.command()
@click.option('--outfile', '-o', default='FoldXData.json',
              type=click.File('wb'),
              help='The database output file')
@click.option('--epitopedb', '-e', default='epitope-info.csv',
              type=click.File('r'),
              help=('A csv with epitope information '
                    '(pySE/data/epitopes/epitope-info.csv)'))
@click.argument('jobs_dir')
def main(outfile, epitopedb, jobs_dir):
    epitopes = generate_epitope_map(load_epitopes(epitopedb))
    parse_foldx_jobs(outfile, epitopes, jobs_dir)
