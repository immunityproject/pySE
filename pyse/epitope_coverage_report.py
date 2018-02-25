#!/usr/bin/env python
# # -*- coding: utf-8 -*-

"""
Epitope Coverage Report

This tool scans output from parse_foldx_jobs.py and generates a report of
each epitope found in the data set.

The report includes the full peptide chain of the epitope with each part
of the chain that is found, the start and stop sites of the subprotein,
and the polyprotein that the subprotein belongs to, if any. It also
includes other available metadata.
"""
from __future__ import print_function

import click
import json
import csv

from collections import defaultdict
from operator import itemgetter

def load_parsed_foldx_data(parsed_foldx_jobs):
    for j in parsed_foldx_jobs:
        with open(j, 'rb') as jobfile:
            for line in jobfile:
                fxjob = json.loads(line.rstrip())
                for fxk, v in fxjob.items():
                    yield v


def get_subprotein_idx(subprotein, length, site):
    parts = subprotein.split('-')
    if ')' not in parts[0]:
        p,s = parts[0].split('(')[:2]
        return site-int(s)
    else:
        p1, start = parts[0].split('(')[:2]
        p2, end = parts[1].split('(')[:2]
        start = int(start[:-1])
        end = int(end[:-1])
        second_site = site - end
        if second_site <= 0: # Second subprotein
            return (site-1) + (length-end)
        else:
            return site - start

@click.command()
@click.argument('parsed_foldx_jobs', nargs=-1)
def main(parsed_foldx_jobs):
    # The field names to copy into the entry and dump to csv
    copyfields = [ 'subprotein', 'peptide', 'protein', 'start', 'end',
                   'epitope', 'subtype' ]
    # (protein,peptide,start) -> info
    ei = defaultdict(dict)

    # The foldx data contains epitope info at each site.
    #
    # The first time we see a given peptide, we add the current site
    # relative to the start position to the peptide chain and add
    # the rest of the peptide info we know about.
    #
    # The second time, we scan the list and update the peptide entry,
    # taking care to ensure the rest of the data remains the same.
    #
    # Note that there is an implied "empty" peptide, with no peptide
    # and a chain with a start at 0 and an end at the last site.
    for foldx_data in load_parsed_foldx_data(parsed_foldx_jobs):
        k = itemgetter('subprotein','peptide','start')(foldx_data)
        # Remap the key to just the protein name if it's empty
        if not k or k == ('','',''):
            k = foldx_data['protein']
        if not ei[k]:
            entry = { k: foldx_data[k] for k in copyfields }
            entry['peptide_fragment'] = ''
            ei[k] = entry

        site, start, end, wt, sub = itemgetter('site', 'start', 'end',
                                               'wt', 'subprotein')(foldx_data)
        frag = ei[k]['peptide_fragment']
        if not site:
            site = 0
        if not start:
            start = 0
        if not end:
            end = 0
        site = int(site)
        start = int(start)
        end = int(end)
        pidx = site - start
        if sub: # subprotein boundary
            pidx = get_subprotein_idx(sub, len(ei[k]['peptide']), site)

        # Note that start and end can both be 0 for entries w/ no epitope
        for i in range(max(len(ei[k]['peptide'])-len(frag),
                           (pidx - len(frag) + 1))):
            frag += '-'

        tmp = list(frag)
        tmp[pidx] = wt
        frag = ''.join(tmp)
        ei[k]['peptide_fragment'] = frag

    with open('out.csv', 'w') as outfile:
        copyfields.insert(1, 'peptide_fragment')
        copyfields.append('mutations')
        writer = csv.DictWriter(outfile, fieldnames=copyfields)
        writer.writeheader()
        for k in ei.keys():
            writer.writerow(ei[k])

if __name__ == '__main__':
    main()
