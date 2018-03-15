#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Map Epitopes to a foldx parsed jobs stream

This tool reads in an epitope mapping database, the original pdb, and
the parsed foldx jobs stream. It starts by generating the peptide
chains from the pdb, then it maps the epitopes from those pdbs against
the provided epitope database. The mapping is done by:

1. majority - picking the chains with the most matches for a given epitope.
2. threshhold - picking the chains with matches of at least the threshhold

After mapping, the stream is read and the mapped epitopes
are added to each site. If multiple epitopes match the same site, the
entry is duplicated at this time. For error tracking, a
'peptide_status' key is added, this indicates what peptide fragments
are missing ('-') or mismatched ('x').
"""

import click
import csv
import json
import logging

from collections import defaultdict

from pyse.pdb import parse_pdb,get_peptide_chains

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

### NOTE: The following 3 functions are unused at this time, but kept
### as they will be used when we add subprotein handling to this code.

def resolve_subprotein_site(pidx, length, subprotein):
    """subproteins come in the format subprotein(start-end) and
    subprotein1(start)-subprotein2(end), the latter requires us to
    determine if we are in the first or second subprotein space."""
    if '(' not in subprotein and ')' not in subprotein:
        return None, None
    parts = subprotein.split('-')
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
    kfmt = '{},{}'
    epmap = defaultdict(list)
    for e in epitopedb:
        start = int(e['start'])
        end = int(e['end'])

        for site in range(start, end + 1):
            pidx = site-start
            # eprint(e['peptide'],pidx, len(e['peptide']), start, end)
            wt = '-'
            if pidx < len(e['peptide']):
                wt = e['peptide'][pidx]
            # protein = e['protein']
            subprotein = e['subprotein']
            subprotein.replace('POL-TF', 'POL_TF')
            k = kfmt.format(wt,site)
            add_epmap_dict(epmap, k, e)
            for p in proteins:
                if p in subprotein:
                    # NOTE: this adds keys outside the protein space,
                    # which is why we keep the full peptide and
                    # subproteins when adding proteins
                    length = end-start
                    sp,subsite = resolve_subprotein_site(pidx,
                                                         length,
                                                         subprotein)
                    if sp.startswith(p):
                        k = kfmt.format(wt,subsite)
                        add_epmap_dict(epmap, k, e)
    # for k,v in epmap.items():
    #     eprint('{}: {}'.format(k, len(v)))
    return epmap

def get_peptide_status(epitope, peptidechain, threshhold=5):
    """Return status if epitope is plausibly in this chain using the
    threshhold method. Since epitopes are usually 9-11 in length,
    anything matching more than 5 (>50%) sequence characters is considered a
    match.

    :param epitope_peptide: A string with the complete epitope peptide sequence
    :param peptide_frag: A string with the epitopes corresponding
                         peptide subsequence
    :return: Status if epitope contains peptide chain, None otherwise
    """
    e_peptide = epitope['peptide']
    # Start is the index into the peptide string for this chain, which
    # is the starting location for the epitope - the starting location
    # of the peptide chain
    start = int(epitope['start']) - int(peptidechain['startsite'])
    end = int(epitope['end']) - int(peptidechain['startsite'])
    if start < 0 or end < 0:
        return None
    if (end - start + 1) != len(e_peptide):
        logging.warn('Skipped! len({}) = {} does not match declared start, end:'
                     ' {}, {} -> {}, {}'.format(e_peptide, len(e_peptide),
                                                epitope['start'],
                                                epitope['end'],
                                                start, end))
        return None

    c_peptide = peptidechain['peptide'][start:(end + 1)]

    logging.debug("{} v {} ({}, {})".format(e_peptide, c_peptide, start, end))
    peptide_status = list()
    hits = 0
    # NOTE: c_peptide can be shorter when epitope describes end of chain
    for i in range(len(c_peptide)):
        if c_peptide[i] == '-':
            # pass through missing peptides
            peptide_status.append('-')
        elif c_peptide[i] == e_peptide[i]:
            hits = hits + 1
            peptide_status.append(e_peptide[i])
        else:
            peptide_status.append('x')

    if hits > threshhold:
        return ''.join(peptide_status)

    return None

def find_epitope_chains(epitopes, peptidechains):
    """For each epitope, find all the peptide chains which may contain
    this epitope. Return a mapping of sites and chains to an epitope
    list that includes a peptide_status key indicating fidelity of the
    match."""
    epitopemap = defaultdict(list)
    for epitope in epitopes:
        for peptidechain in peptidechains:
            peptide_status = get_peptide_status(epitope, peptidechain)
            startsite = int(epitope['start'])
            endsite = int(epitope['end'])
            length = endsite - startsite
            if not peptide_status:
                subprotein = epitope['subprotein']
                sp,subsite = resolve_subprotein_site(0,
                                                     length,
                                                     subprotein)
                if sp and subsite:
                    # NOTE: Updating startsite is important to get the
                    # key right map key, wihch has to be relative to
                    # the PDB, not the epitope info
                    startsite = subsite
                    epitope['start'] = subsite
                    epitope['end'] = subsite + length
                    peptide_status = get_peptide_status(epitope, peptidechain)

            if not peptide_status:
                continue

            if epitope['peptide'] == 'MHEDIISLW':
                print('LOVE: {}'.format(peptide_status))

            for i in range(length):
                mapkey = '{}{}'.format(peptidechain['chain'],i + startsite)
                if epitope['peptide'] == 'MHEDIISLW':
                    print('LOVE: {}'.format(mapkey))
                new_e = epitope.copy()
                new_e['peptide_status'] = peptide_status
                epitopemap[mapkey].append(new_e)
    return epitopemap

def get_epitopes(epitopemap, site, chains):
    """Search the epitopemap for epitopes matching this site and chain"""
    epitopes = list()
    ekeys = set()
    for c in chains.split(','):
        mapkey = '{}{}'.format(c,site)
        for e in epitopemap[mapkey]:
            ekey = '{}{}{}'.format(e['peptide'], e['epitope'], e['start'])
            if ekey not in ekeys:
                epitopes.append(e)
                ekeys.add(ekey)
    return epitopes

@click.command()
@click.option('--epitopecsv', '-e', default='epitope-info.csv',
              type=click.File('r'),
              help=('A csv with epitope information '
                    '(pySE/data/epitopes/epitope-info.csv)'))
@click.option('--pdb', '-p', default='protein.pdb',
              type=click.File('r'),
              help=('A pdb file that was used for (at least part of) '
                    'the input data stream.'))
@click.option('--infile', '-i', default='-',
              type=click.File('r'),
              help=('Input file in jsonl format'))
@click.option('--outfile', '-o', default='-',
              type=click.File('w'),
              help=('Output file in jsonl format'))
@click.option('--protein', '-f', default=None,
              help=('Limit to given protein name'))
def map_epitopes(epitopecsv, pdb, infile, outfile, protein):
    peptidechains = get_peptide_chains(parse_pdb(pdb))
    epitopes = load_epitopes(epitopecsv)

    if protein:
        new_eps = list()
        for e_v in epitopes:
            # skip this entry if we fail to get a protein match on the
            # protein or the subprotein. NOTE the startswith on
            # subprotein, because they often contain the site
            # information.
            if (protein != e_v['protein']
                and not e_v['subprotein'].startswith(protein)):
                continue
            new_eps.append(e_v)
        epitopes = new_eps

    matched_epitopes = find_epitope_chains(epitopes, peptidechains)
    for line in infile:
        foldx_job = json.loads(line.rstrip())
        site = ''
        chains = ''
        fkey = foldx_job.keys()[0]
        site = foldx_job[fkey]['site']
        chains = foldx_job[fkey]['chains']
        job_epitopes = get_epitopes(matched_epitopes, site, chains)
        if not job_epitopes:
            job_epitopes.append({ k: '' for k in [ "peptide", "protein",
                                                   "subprotein",
                                                   "start", "end", "epitope",
                                                   "subtype" ]})
        for e in job_epitopes:
            for k,v in e.items():
                if k not in ['protein']:
                    foldx_job[fkey][k] = v
            outfile.write((json.dumps(foldx_job) + '\n').encode())

if __name__ == '__main__':
    map_epitopes()
