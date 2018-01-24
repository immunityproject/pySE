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

def parse_foldx_jobs(outfile, jobs_dir):
    foldx_jobs = list(find_foldx_jobs(jobs_dir))
    jobkeys = set()

    workers = Pool(cpu_count())
    with click.progressbar(workers.imap_unordered(load_foldx_job,
                                                  foldx_jobs),
                           length=len(foldx_jobs), label='Parsing',
                           file=sys.stdout) as progbar:
        for j in progbar:
            if j:
                for key in j.keys():
                    jobkeys.add(key)
                outfile.write((json.dumps(j) + '\n').encode())
                outfile.flush()
                j.clear()
                gc.collect()

    return jobkeys

def extract_peptide_mapping(jobkeys):
    """For each job key, extract the mapping of protein to the list of
    sites to their wild types.

    This mapping can be use to reconstruct peptide chains using the site order.
    """
    pm = defaultdict(dict)
    for jk in jobkeys:
        protein, wildtype, site = jk.split(',')[0:3]
        site = int(site)
        wtlist = pm[protein].get(site, set())
        wtlist.add(wildtype)
        pm[protein][site] = wtlist

    return pm

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

def find_protein_name(protein, subprotein, proteinslist):
    """Given a protein and subprotein name, determine if any of the
    proteins in the proteins list match"""
    if protein in proteinslist:
        return protein
    if subprotein in proteinslist:
        return subprotein
    for p in proteinslist:
        if subprotein.startswith(p):
            return p
    return None

def find_jobkey_epitopes(epitopes, jobkeys):
    """For each epitope in the epitopes list, map it to it's
    corresponding job keys"""
    peptidemap = extract_peptide_mapping(jobkeys)
    proteins = peptidemap.keys()

    # Map matched epitopes to job prefixes
    jkprefixes = defaultdict(list)
    for e in epitopes:
        # Find the protein name in peptidemap
        protein = find_protein_name(e['protein'], e['subprotein'], proteins)
        if not protein:
            continue

        # Now find the start and stop
        pm = peptidemap[protein]
        start = int(e['start'])
        end = int(e['end'])
        if start not in pm or end not in pm:
            continue

        # Determine all peptide sequence candidates in this site range
        candidate_peptides = [ '' ] # Start with 1 empty candidate
        for site in range(start, end + 1):
            new_candidates = []
            # Default case: do not add candidates
            # If more than 1, replicate each existing candidate
            for cp in candidate_peptides:
                for wt in pm[site]:
                    new_candidates.append(cp + wt)
            candidate_peptides = new_candidates

        # If the peptide for this epitope is present in these sites,
        # add it to each site
        peptide = e['peptide']
        if peptide not in candidate_peptides:
            continue
        for site in range(start, end + 1):
            # Job keys start with protein,wildtype,site
            for wt in pm[site]:
                key = '{},{},{}'.format(protein, wt, site)
                jkprefixes[key].append(e)

    # For each job key, set the epitope based on the matches found in the prefix
    jobkey_epitopes = dict()
    for jk in jobkeys:
        p,s,wt = jk.split(',')[0:3]
        jkprefix = '{},{},{}'.format(p, s, wt)
        jobkey_epitopes[jk] = jkprefixes[jkprefix]
    return jobkey_epitopes

@click.command()
@click.option('--outfile', '-o', default='FoldXData.json',
              help='The database output file')
@click.option('--epitopedb', '-e', default='epitope-info.csv',
              type=click.File('r'),
              help=('A csv with epitope information '
                    '(pySE/data/epitopes/epitope-info.csv)'))
@click.argument('jobs_dir')
def main(outfile, epitopedb, jobs_dir):
    epitopes = load_epitopes(epitopedb)

    jobkeys = set()
    with open(outfile + '.tmp', 'wb') as out:
        jobkeys = parse_foldx_jobs(out, jobs_dir)

    jobkey_epitopes = find_jobkey_epitopes(epitopes, jobkeys)

    with open(outfile + '.tmp', 'rb') as tmp:
        with open(outfile, 'wb') as out:
            for line in tmp:
                entry = json.loads(line.rstrip())
                for k in entry.keys():
                    matched_eps = jobkey_epitopes.get(k, ['Unknown'])
                    for matched_ep in matched_eps:
                        entry[k]['epitope'] = matched_ep['epitope']
                        entry[k]['peptide'] = matched_ep['peptide']
                        entry[k]['subprotein'] = matched_ep['subprotein']
                        entry[k]['polyprotein'] = matched_ep['protein']
                        entry[k]['subtype'] = matched_ep['subtype']
                        out.write((json.dumps(entry) + '\n').encode())
                out.flush()

    os.remove(outfile + '.tmp')
