#!/usr/bin/python
# coding=utf-8

# This is the final analysis script.  It reads data from the FoldX
# jobs directory and can compute various things depending on
# command-line switches:
# - Boltzmann Distribution for all mutations,
#   the corresponding Site Shannon Entropies, and the epitope Structural
#   Entropies
# - Structural Entropy for all sites (--scan)
# - 3-D displacements due to mutations (experimental, --displacement)

from __future__ import print_function
import functools

import math
import os
import json
import re
import sys
import csv
import click
import datetime
from traceback import print_exc
from multiprocessing import Pool, Manager, cpu_count

import bigfloat
from bigfloat import BigFloat

from proteins import proteins, pdbs, codes

from collections import defaultdict

# Define and error printing function in lieu of logging api
eprint = functools.partial(print, file=sys.stderr)


# In trying to replicate the results of Pereyra et al, we explored
# dividing the change in free energy by the symmetry number of a
# rotamer.  For example, a setting of rotamerSymmetriesToTry = [1,6]
# will calculate results using both the total free energy deltas and
# the total energy deltas divided by 6
rotamerSymmetriesToTry = [1]

# Trying to figure out units in the original paper, we explored
# dividing the energies by the Boltzmann Constant times temperature.
# This affects the magnitude of the elements in the Boltzmann
# Distribution, but it does not impact the final Shannon Entropy
# calculations.  Setting divideByKTToTry = [True,False] will produce
# results with and without the division by the kT constant
divideByKTToTry = [False]

site_summary = {}
comparison = []


# Constants
boltzmannConstant = 0.0019872041  # kcal/(mole.Kelvin)
temperature = 298  # Kelvin
# We are dealing with small values that _should_ fit in an IEEE double
# precision float, but we use BigFloat throughout for certainty
bigfloatPrecision = 100

# The basic approach of this script is to recursively read all
# directories in the jobs_dir and mark whether they belong to the
# protein/site range/epitope selected for analysis.  Each selected job
# directory is then represented by a MutationEvaluator object which
# performs lazy evaluation of various values for the amino acid
# substitution represented by that job directory.  The Structural
# Entropy (or displacement) for each epitope or site range is accessed
# by invoking the relevant MutationEvaluators, which use cached
# computations if the site has already been processed by another
# epitope or site range.  Since many epitopes overlap, this caching
# structure speeds up computation considerably.
@click.command()
@click.option('--debug/--no-debug', default=False,
              help='Show script debugging information')
@click.option('--quiet/--no-quiet', '-q', default=False,
              help='Suppress warnings about missing data')
@click.option('--protein', '-p', type=click.Choice(proteins.keys()),
              help=('Limit analysis to a specific protein (used in conjunction'
                    ' with --sites)'))
@click.option('--sites', '-s',
              help=('Site range to analyze, using a dash range: "115-121".'
                    ' Ranges are inclusive. Must also specify --protein'))
@click.option('--epitope', '-e',
              help=('Limit computation to a specific epitope. Shortcut for'
                    ' setting --protein and --site'))
# experimental mode to compute the sum of squared displacements for each epitope,
@click.option('--displacement/--no-displacement', '-d', default=False,
              help=('Use displacement calculation instead of energies.'
                    ' Requires a flag which sets the protein.'))
@click.option('--baseline', '-b', type=click.Choice(['absolute', 'raw']),
              default='absolute',
              help='How to preprocess energies before entropy calculation')
@click.option('--exclude-wt/--no-exclude-wt', default=False,
              help=('If set, the Boltzmann distribution does not include the'
                    ' WT energy difference (which is always 0.0)'))
@click.option('--csv-file', '-c', default=None,
              help=('Base pathname (without extension) for a CSV dump of the'
                    ' results'))
@click.option('--json/--no-json', default=False,
              help='Generate json results in addition to csv')
@click.option('--header/--no-header', default=True,
              help='Include or exclude the header row from the CSV files')
@click.option('--tab-delimited', default=False,
              help='Use tab delimiters instead of commas in the CSV files')
@click.option('--scan', type=int,
              help=('Computes SE\'s for a sliding window across all sites.'
                    ' This is is the number of sites in the sliding window.'))
@click.option('--energy-map', default=False,
              help=('Compute average energy deltas of all possible amino acid'
                    ' substitutions vs WT'))
@click.option('--energies', default=False,
              help='outputs the raw energy delta data')
@click.argument('jobs_dir')
def main(debug, quiet, protein, sites, epitope, displacement, baseline,
         exclude_wt, csv_file, json, header, tab_delimited, scan, energy_map,
         energies, jobs_dir):
    include_wt = not exclude_wt
    print("Jobs dir...{}".format(datetime.datetime.now()))
    sys.stdout.flush()
    initJobsDirs(jobs_dir, include_wt, displacement, debug)
    print("Done...{}".format(datetime.datetime.now()))
    sys.stdout.flush()

    # energy map is a separate mode
    if energy_map:
        energy_map(protein, epitope, baseline, csv_file)
        exit(0)

    # the results array holds tuples of result data for output either
    # as CSV or to console
    results = []
    has_epitopes = False
    has_ranges = False

    # this sets up the output files
    if csv_file:
        setup_entropy_csv(csv_file, tab_delimited)

    # scan is a separate mode
    if scan:
        results.extend(scan(protein, scan, debug, baseline, displacement,
                            include_wt))
        has_ranges = True

    elif epitope is not None:
        if protein is None:
            eprint('If you specify --epitope, you must also specify --protein')
            click.get_current_context().get_help()
            return
        # EpitopeResults handles the cached calculations
        results.append(EpitopeResults(protein, epitope, debug, baseline,
                                      displacement, include_wt))
        has_epitopes = True

    elif protein:
        # ProteinResults handles the cached calculations
        results.append(ProteinResults(protein, debug, baseline, displacement,
                                      include_wt))
        has_epitopes = True

    elif sites:
        if not protein:
            eprint('If you specify --sites, you must also specify --protein')
            click.get_current_context().get_help()
            return
        first, last = (int(i) for i in sites.split('-'))
        # RangeResults handles the cached calculations
        results.append(RangeResults(protein, first, last, debug, baseline,
                                    displacement, include_wt))

    else:  # all proteins & epitopes
        has_epitopes = True

        print("Running all proteins...{}".format(datetime.datetime.now()))

        for proteinName in proteins.keys():
            results.append(ProteinResults(proteinName, debug, baseline,
                                          displacement, include_wt))

    if energies:  # output raw energy data
        if csv_file:
            output_energies_csv(results, tab_delimited, csv_file)
        else:
            outputEnergies(results)
        return

    # output results
    if csv_file:
        output_csv(results, has_epitopes, has_ranges, tab_delimited, csv_file,
                   header)
    else:
        output(results)

    if json:
        output_json(results, csv_file)

def dump_energies(range_results):
    rr = range_results
    energies = defaultdict(dict)
    for site in rr.siteResults:
        for mutation in site.evaluators:
            mutations = energies[mutation.protein].get(mutation.site, [])
            mutations.append({
                "wt": str(mutation.wt),
                "mutation": str(mutation.mutation),
                "energyDelta": str(mutation.energyDelta)
            })
            energies[mutation.protein][mutation.site] = mutations
    return energies

def dump_entropies(site_results):
    return sorted([str(sr.entropy) for sr in site_results
                   if sr.entropy is not None])

def dump_displacements(site_results):
    return sorted([str(sr.averageDisplacement) for sr in site_results
                   if sr.averageDisplacement is not None])

def dump_SiteResult(site_result):
    sr = site_result
    mutations = defaultdict(dict)
    mutations['entropy'] = str(sr.entropy)
    for mutation in sr.evaluators:
        idx = sr.evaluators.index(mutation)
        mutations[mutation.wt][mutation.mutation] = {
            "energy_delta": str(mutation.energyDelta),
            "displacement": str(mutation.displacement),
            "boltzman_probability": str(sr.boltzmann_probabilities[idx])
        }

    return mutations

def dump_EpitopeResult(epitope_result, output):
    er = epitope_result

    output[er.proteinName][er.epitopeName] = { sr.site: dump_SiteResult(sr)
                                               for sr in er.siteResults }
    # {
    #     "avgE": str(er.averageEntropy),
    #     "entropies": dump_entropies(er.siteResults),
    #     "avgD": str(er.averageDisplacement),
    #     "displacements": dump_displacements(er.siteResults),
    #     "energies": dump_energies(epitope_result)
    # }
    return output

def dump_ProteinResult(protein_result, output):
    """Update and return the given output dict with the given protein_result"""
    for er in protein_result.epitopeResults:
        output = dump_EpitopeResult(er, output)
    return output

def dump_RangeResult(range_result, output):
    rr = range_result
    output[rr.proteinName]['{}-{}'.format(rr.first, rr.last)] = {
        "avgE": str(rr.averageEntropy),
        "entropies": dump_entropies(rr.siteResults),
        "avgD": str(rr.averageDisplacement),
        "displacements": dump_displacements(rr.siteResults),
        "energies": dump_energies(range_result)
    }
    return output


def output_json(results, filename):
    """Output all results as machine-parseable json files

    :param results: the results array
    :param filename: the output filename prefix

    FIXME: The following is very naive and unlikely to work.
    """
    outputs = defaultdict(dict)
    # If the results are one of the various results classes, parse into a dict
    # FIXME: We could improve this by aggregating results into a single object
    #        or by passing an output buffer to print results if they
    #        don't fit in memory
    for r in results:
        if isinstance(results[0], SiteResults):
            raise IOError('Unimplemented!')
        elif isinstance(results[0], EpitopeResults):
            outputs = dump_EpitopeResult(r, outputs)
        elif isinstance(results[0], RangeResults):
            outputs = dump_RangeResult(r, outputs)
        elif isinstance(results[0], ProteinResults):
            outputs = dump_ProteinResult(r, outputs)
        else:
            outputs.append(r)

    with open('{}.json'.format(filename), 'wb') as outfile:
        json.dump(outputs, outfile)

csv_entropies_file = None  # file object
csv_site_results = None  # csv writer


def setup_entropy_csv(csv_file, tab_delimited):
    global csv_site_results, csv_entropies_file
    if tab_delimited:
        dialect = 'excel-tab'
    else:
        dialect = 'excel'
    csv_entropies_file = open(csv_file + '-site_entropies.csv', 'wb')
    csv_site_results = csv.writer(csv_entropies_file, dialect=dialect)
    csv_site_results.writerow(('protein', 'site', 'entropy'))


def scan(protein, scan, debug, baseline, displacement, include_wt):
    """runs through all sites within the ranges specified by
    command-line args and accumulates the calculation results"""
    if protein:
        return scan_protein(protein, scan, debug, baseline, displacement,
                            include_wt)

    results = []
    for protein in proteins.keys():
        results.extend(scan_protein(protein, scan, debug, baseline,
                                    displacement, include_wt))
    return results


def energy_map(protein, epitope, baseline, csv_file):
    bannedList = ['X', 'B', 'Z']  # these aminos are ignored
    aaCodes = sorted([v for v in codes.values() if v not in bannedList])
    energies = {wt: {mutation: [] for mutation in aaCodes} for wt in aaCodes}
    for evaluator in evaluators.values():
        if evaluator.wt in bannedList or evaluator.mutation in bannedList:
            continue
        if protein and evaluator.protein != protein:
            continue
        if epitope:
            epitopeRecord = proteins[evaluator.protein][epitope]
            if (evaluator.site < epitopeRecord.first
                or evaluator.site > epitopeRecord.last):
                continue
        try:
            energy = evaluator.energyDelta
        except Exception:
            print_exc()
            print('energy_map is ignoring the bad data.')
            continue
        if baseline == 'absolute':
            energy = abs(energy)
        energies[evaluator.wt][evaluator.mutation].append(energy)
    if csv_file:
        with open(csv_file, 'wb') as f:
            out = csv.writer(f)
            out.writerow(['WT \ Mutation'] + aaCodes)
            for wt in aaCodes:
                out.writerow([wt] + [mean(energies[wt][mutation])
                                     for mutation in aaCodes])
    else:
        print('WT \ Mutation \t' + '\t'.join(aaCodes))
        for wt in aaCodes:
            # print wt, < -- not sure if below is equivalent -- FIXME
            print(wt, end=',')
            for mutation in aaCodes:
                # print '\t' + str(),
                # FIXME: Not sure if the following is equivalent to above
                print('\t', end=',')
            print


def scan_protein(proteinName, scan, debug, baseline, displacement, include_wt):
    proteinKey = proteinName.upper()
    minIndex, maxIndex = sys.maxint, 0
    for protein, site, mutation in evaluators.keys():
        if protein == proteinKey:
            minIndex = min(minIndex, site)
            maxIndex = max(maxIndex, site)

    results = []
    window_size = scan
    for first in range(minIndex, maxIndex - window_size + 2):
        results.append(RangeResults(
            proteinName, first, first + window_size - 1, debug, baseline,
            displacement, include_wt))
    return results


def output_csv(results, has_epitopes, has_ranges, tab_delimited, csv_file,
               csv_header):
    csv_epitopes = None
    global csv_site_results
    if tab_delimited:
        dialect = 'excel-tab'
    else:
        dialect = 'excel'
    epitope_file = None
    try:
        if has_epitopes:
            epitope_file = open(csv_file +
                                '-epitope_summaries.csv', 'wb')
            csv_epitopes = csv.writer(epitope_file, dialect=dialect)
        elif has_ranges:
            epitope_file = open(csv_file +
                                '-range_summaries.csv', 'wb')
            csv_epitopes = csv.writer(epitope_file, dialect=dialect)
        if csv_header:
            if has_epitopes:
                csv_epitopes.writerow(('protein', 'epitope', 'avg_entropy',
                                       'max_entropy', 'max_2_entropy',
                                       'max_3_entropy',
                                       'avg_displacement', 'min_displacement',
                                       'min_2_displacement',
                                       'min_3_displacement'))
            elif has_ranges:
                csv_epitopes.writerow(('protein', 'first', 'last',
                                       'avg_entropy', 'max_entropy',
                                       'max_2_entropy',
                                       'max_3_entropy', 'avg_displacement',
                                       'min_displacement', 'min_2_displacement',
                                       'min_3_displacement'))
        for result in results:
            if isinstance(result, ProteinResults):
                for er in result.epitopeResults:
                    output_epitope_result_csv(er, csv_epitopes)
            elif isinstance(result, EpitopeResults):
                output_epitope_result_csv(result, csv_epitopes)
            elif isinstance(result, RangeResults):
                output_range_result_csv(result, csv_epitopes)
    finally:
        csv_entropies_file.close()
        if epitope_file is not None:
            epitope_file.close()


def output(results):
    for result in results:
        if isinstance(result, ProteinResults):
            for er in result.epitopeResults:
                output_epitope_result(er)
            if result.r2:
                print('%s r2 = %.2f' % (result.proteinName, result.r2))
        elif isinstance(result, EpitopeResults):
            output_epitope_result(result)
        elif isinstance(result, RangeResults):
            output_range_result(result)


def output_epitope_result(er):
    if er.averageEntropy is not None:
        print('%s %s entropy = %.2f' % (er.proteinName, er.epitopeName,
                                        er.averageEntropy), end=',')
        if er.expectedAverageEntropy is not None:
            print('(%.2f)' % er.expectedAverageEntropy, end=',')
    if er.averageDisplacement is not None:
        print('with average displacement of %.4f' % er.averageDisplacement,
              end=',')
    print


def output_range_result(rr):
    if rr.averageEntropy is not None:
        print('%s %d-%d entropy = %.2f' % (rr.proteinName, rr.first, rr.last,
                                           rr.averageEntropy), end=',')
    if rr.averageDisplacement is not None:
        print('with average displacement of %.4f' % rr.averageDisplacement)
    print('')


def output_epitope_result_csv(er, out):
    avgE = er.averageEntropy or ''
    entropies = sorted(
        [sr.entropy for sr in er.siteResults if sr.entropy is not None],
        reverse=True)
    while len(entropies) < 3:
        entropies += [None]
    avgD = er.averageDisplacement or ''
    displacements = sorted(
        [sr.averageDisplacement for sr in er.siteResults
         if sr.averageDisplacement is not None])
    while len(displacements) < 3:
        displacements += [None]
    row = (er.proteinName, er.epitopeName,
           avgE, entropies[0], entropies[1], entropies[2],
           avgD, displacements[0], displacements[1], displacements[2])
    out.writerow(row)


def output_range_result_csv(rr, out):
    avgE = rr.averageEntropy or ''
    entropies = sorted(
        [sr.entropy for sr in rr.siteResults if sr.entropy is not None],
        reverse=True)
    while len(entropies) < 3:
        entropies += [None]
    avgD = rr.averageDisplacement or ''
    displacements = sorted(
        [sr.averageDisplacement for sr in rr.siteResults
         if sr.averageDisplacement is not None])
    while len(displacements) < 3:
        displacements += [None]
    row = (rr.proteinName, rr.first, rr.last,
           avgE, entropies[0], entropies[1], entropies[2],
           avgD, displacements[0], displacements[1], displacements[2])
    out.writerow(row)


def output_energies_csv(results, tab_delimited, csv_file):
    if tab_delimited:
        dialect = 'excel-tab'
    else:
        dialect = 'excel'
    with open(csv_file + '-energies.csv', 'wb') as outfile:
        csvfile = csv.writer(outfile, dialect=dialect)
        csvfile.writerow(('protein', 'site', 'wt', 'mutation', 'deltaE'))

        def writeline(items):
            csvfile.writerow(items)
        outputEnergies(results, writeline)


def writeStdout(items):
    print('\t'.join(items))


def outputEnergies(results, writeline=writeStdout):
    writeline_queue = []

    def handleRangeResults(result):
        for site in result.siteResults:
            for mutation in site.evaluators:
                #items =
                # (mutation.protein,mutation.site,mutation.wt,mutation.mutation,
                # mutation.energyDelta)
                # writeline(items)
                items = "{0},{1},{2},{3},{4}".format(
                    mutation.protein, mutation.site, mutation.wt,
                    mutation.mutation, mutation.energyDelta)
                writeline_queue.append(items)
    for result in results:
        if isinstance(result, ProteinResults):
            for eresult in result.epitopeResults:
                handleRangeResults(eresult)
        else:
            handleRangeResults(result)

    # sort results
    print(writeline_queue)
    sorted_wl_queue = sorted(writeline_queue, key=str.lower)
    # writeline(sorted_wl_queue)
    for entry in sorted_wl_queue:
        writeline(tuple(entry.split(",")))


proteinToPdb = {protein: pdbPath.split(
    '/')[-1] for protein, pdbPath in pdbs.iteritems()}


class MutationEvaluator(object):
    """ evaluates a single job directory which represents a single
    site mutation """

    def __init__(self, directory, protein, site, mutation, wt, displacement,
                 debug):
        self.initialized = False
        self.directory = directory
        self.protein = protein
        self.site = site
        self.mutation = mutation
        self.wt = wt
        self.displacement = displacement
        self.debug = debug
        self.initialize()

    def initialize(self):
        if self.initialized:
            return
        pdb = proteinToPdb[self.protein]
        self.pdb = pdb
        pdbStub = pdb[:-len('.pdb')]
        # each entry is a list whose first element is the total energy and
        # other elements are the component energies
        self.energies = []
        # each entry is a list whose first element is the total energy
        # and other elements are the component energies
        self.wtEnergies = []
        self.displacements = []

        # load energies
        raw_buildmodel_filename = os.path.join(
            self.directory, 'Raw_BuildModel_%s.fxout' % pdbStub)
        try:
            with open(raw_buildmodel_filename) as buildmodel_file:
                for line in buildmodel_file.readlines()[9:]:
                    if line.startswith(pdbStub):
                        self.energies.append([float(e)
                                              for e in line.split('\t')[1:]])
                    elif line.startswith('WT_'):
                        self.wtEnergies.append([float(e)
                                                for e in line.split('\t')[1:]])
        except IOError as err:
            eprint('Could not read file %s' % raw_buildmodel_filename)
            raise err
        assert len(self.energies) == len(self.wtEnergies)
        # average deltas for each energy component
        self.energyDeltas = [mean([e - w for e, w in zip(es, ws)])
                             for es, ws in zip(zip(*(self.energies)),
                                               zip(*(self.wtEnergies)))]
        self.energyDelta = self.energyDeltas[0]

        if self.displacement:
            # calculate displacements
            i = 0
            while True:
                outputPdb = os.path.join(
                    self.directory, '%s_1_%d.pdb' % (pdbStub, i))
                if not os.path.exists(outputPdb):
                    break
                wtOutputPdb = os.path.join(
                    self.directory, 'WT_%s_1_%d.pdb' % (pdbStub, i))
                with open(outputPdb) as a, open(wtOutputPdb) as b:
                    self.displacements.append(analyze_displacement(a, b,
                                                                   self.debug))
                i += 1
            if len(self.displacements) == 0:
                eprint('Warning: no displacement data found for'
                       ' %s %s %s' % (self.protein, self.site, self.mutation))
                self.displacement = None
            else:
                self.displacement = mean(self.displacements)
        else:
            self.displacements = []
            self.displacement = None

        self.initialized = True


manager = Manager()
evaluators = manager.dict()
mutations = codes.values()


def initJobsDirs(jobs_dir, include_wt, displacement, debug):
    """ initializes the jobDirs variable """
    pdbToProtein = { pdbPath.split('/')[-1]: protein
                     for protein, pdbPath in pdbs.iteritems() }

    roots = list()
    filelists = list()
    for root, dirs, files in os.walk(jobs_dir, followlinks=True):
        roots.append(root)
        filelists.append(files)
        # checkJobsDir(pdbToProtein, root, files, include_wt, displacement, debug)

    doers = Pool(cpu_count())
    jobs = zip([pdbToProtein]*len(roots),
               roots,
               filelists,
               [include_wt]*len(roots),
               [displacement]*len(roots),
               [debug]*len(roots))

    print("Loading {} evaluators...{}".format(len(roots),
                                              datetime.datetime.now()))
    sys.stdout.flush()
    with click.progressbar(doers.imap_unordered(checkJobsDir, jobs),
                           length=len(roots), label='Running',
                           file=sys.stderr) as progbar:
        for j in progbar:
            pass
    print("Done...{}".format(datetime.datetime.now()))




def checkJobsDir(args):
    """Looks for FoldX job files in the directory and creates a
    MutationEvaluator for the directory"""
    pdbToProtein, dirname, filenames, include_wt, displacement, debug = args
    list_file_path = os.path.join(dirname, 'list.txt')
    if not os.path.exists(list_file_path):
        # no list.txt file, so this is not a job directory.
        return
    with open(list_file_path) as list_file:
        pdbFilename = list_file.read().strip()
        protein = pdbToProtein.get(pdbFilename)
        if protein is None:
            raise Exception('%s was not found in the pdbs variable in'
                            ' proteins.py' % pdbFilename)
    ind_list_file_path = os.path.join(dirname, 'individual_list.txt')
    with open(ind_list_file_path) as ind_list_file:
        m = re.match(r'(\w)(\w)(\d+)(\w)',
                     ind_list_file.readline().split(',')[0])
        wt = m.group(1)
        # chain = m.group(2)
        site = int(m.group(3))
        mutation = m.group(4)
    if mutation != wt or include_wt:
        evaluatorKey = (protein, site, mutation)
        if evaluatorKey in evaluators:
            eprint('Warning: duplicate data found for %s %d %s in directories '
                   '%s and %s' % (protein, site, mutation, dirname,
                                  evaluators[evaluatorKey].directory))
        else:
            try:
                evaluators[evaluatorKey] = MutationEvaluator(
                    dirname, protein, site, mutation, wt, displacement, debug)
            except Exception as e:
                eprint("Exception: {}".format(e))


def analyze_displacement(original_pdb, mutated_pdb, debug):
    """ This mode computes the sum of squared displacements across all
    sites in the protein due to the mutation.

    :param original_pdb: file-like of the original PDB file
    :param mutated_pdb: file-like of the mutated PDB file
    :return: a single non-negative float which represents the sum of
             squared displacements caused by the mutation
    """

    orig_pos = read_positions(original_pdb, debug)
    mutated_pos = read_positions(mutated_pdb, debug)

    total = BigFloat(0)
    matches = 0
    for k, v in orig_pos.iteritems():
        if k in mutated_pos:
            matches += 1
            ox, oy, oz = v
            mx, my, mz = mutated_pos[k]
            dx, dy, dz = BigFloat(
                mx - ox), BigFloat(my - oy), BigFloat(mz - oz)
            total += dx * dx + dy * dy + dz * dz
    if debug:
        print('found %d matching ATOM lines between PDB files' % matches)
    return float(total)


def read_positions(file, debug):
    """ used by analyze_displacement.  reads x,y,z positions for each
    atom in a PDB file """
    result = {}

    def stripstr(v):
        return str(v).strip()
    for line in file:
        if line.startswith('ATOM'):
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
            info = {f[0]: f[3](line[f[1] - 1:f[2]]) for f in fieldspec}
            key = '%s-%s-%s-%s' % (info['chain'], info['position'],
                                   info['remnant'], info['atom'])
            result[key] = (info['x'], info['y'], info['z'])
    if len(result) == 0:
        print('WARNING: displacement file had no ATOM lines')
    if debug:
        print('displacement file had ' + str(len(result)) + ' ATOM lines')
    return result


def mean(list):
    num = len(list)
    if num == 0:
        return float('NaN')
    return sum(list) / num


class SiteResults:
    """This uses the MutationEvaluators for the current site to
    calculate site-wide metrics like its Boltzmann Distribution"""

    instances = {}

    @classmethod
    def instance(cls, proteinName, site, debug, baseline, displacement,
                 include_wt):
        proteinKey = proteinName.upper()
        key = (proteinKey, site)
        found = cls.instances.get(key)
        if found is None:
            found = SiteResults(proteinName, site, debug, baseline,
                                displacement, include_wt)
            cls.instances[key] = found
        return found

    def __init__(self, proteinName, site, debug, baseline, displacement,
                 include_wt):
        self.proteinName = proteinName
        self.site = site
        self.evaluators = []  # MutationEvaluators
        proteinKey = proteinName.upper()
        if debug:
            print('Site %s %d Results' % (proteinName, site))
        # Collect all MutationEvaluators for this site
        for mutation in mutations:
            if mutation not in ['B', 'Z', 'X']:
                evaluator = evaluators.get((proteinKey, site, mutation))
                if evaluator is not None:
                    self.evaluators.append(evaluator)
        if debug:
            for ev in self.evaluators:
                print('%s %s %s energy: %.3f' % (ev.protein, ev.site,
                                                 ev.mutation, ev.energyDelta))
        # Collect energy deltas for the Boltzmann Distribution
        # todo different baselines.  currently only ddE
        self.energies = [ev.energyDelta for ev in self.evaluators]
        if baseline == 'absolute':
            self.energies = [abs(e) for e in self.energies]
        elif baseline == 'raw':
            pass
        else:
            eprint('Unknown baseline method: %s' % baseline)
            click.get_current_context().get_help()
            raise IOError('Unknown baseline method: %s' % baseline)
        # Displacements
        self.displacements = [ ev.displacement for ev in self.evaluators
                               if ev.displacement is not None ]
        if len(self.displacements):
            self.averageDisplacement = mean(self.displacements)
        else:
            if displacement:
                eprint('Warning: no displacement data for %s %d'
                       % (self.proteinName, self.site))
            self.averageDisplacement = None
        # Compute boltzmann distribution from energy deltas
        if include_wt:
            hasWt = False
            for ev in self.evaluators:
                if ev.mutation == ev.wt:
                    hasWt = True
                    break
            if len(self.energies) and not hasWt:
                self.boltzmann_probabilities = compute_boltzmann(
                    self.energies + [0.0], debug)
            else:
                self.boltzmann_probabilities = compute_boltzmann(self.energies,
                                                                 debug)
        else:
            self.boltzmann_probabilities = compute_boltzmann(self.energies,
                                                             debug)
        if len(self.boltzmann_probabilities):
            self.entropy = compute_entropy(self.boltzmann_probabilities)
        else:
            self.entropy = None
        if debug and self.entropy is not None:
            print('entropy %.3f' % self.entropy)
        # Validate that we have one evaluator/job directory for every
        # amino acid mutation
        num_evaluators = len(self.evaluators)
        if (num_evaluators != 19 and not include_wt or num_evaluators != 20
            and include_wt):
            eprint('Warning: %s site %d only had %d mutations: %s'
                   % (proteinName, site, num_evaluators,
                      [ ev.mutation for ev in self.evaluators ]))
        # If CSV is selected, output per-site data to the site results file
        if csv_site_results is not None:
            csv_site_results.writerow(
                (proteinName, site, self.entropy, self.averageDisplacement))


class RangeResults:
    """Computes Structural Entropy for a range of sites"""

    def __init__(self, proteinName, first, last,
                 debug, baseline, displacement, include_wt):
        self.proteinName = proteinName
        self.siteResults = []
        self.first = first
        self.last = last
        # site ranges are inclusive by convention
        for site in range(first, last + 1):
            try:
                self.siteResults.append(
                    SiteResults.instance(proteinName, site, debug, baseline,
                                         displacement, include_wt))
            except Exception:
                eprint('Warning: ignoring %s %d due to errors'
                       % (proteinName, site))
        # Collect entropies from every site in the range
        entropies = [ ev.entropy for ev in self.siteResults
                      if ev.entropy is not None ]
        if len(entropies) > 0:
            self.averageEntropy = mean(entropies)
        else:
            self.averageEntropy = None

        # Average displacement for a range
        self.displacements = [
            ev.averageDisplacement for ev in self.siteResults
            if ev.averageDisplacement is not None]
        if len(self.displacements) > 0:
            self.averageDisplacement = mean(self.displacements)
        else:
            self.averageDisplacement = None


class EpitopeResults(RangeResults):
    """Uses RangeResults to compute Structural Entropy for a given epitope"""

    def __init__(self, proteinName, epitopeName,  debug, baseline, displacement,
                 include_wt):
        self.proteinName = proteinName
        self.epitopeName = epitopeName
        self.epitopeKey = epitopeName.upper()
        self.proteinKey = proteinName.upper()
        protein = get_protein(proteinName)
        epitope = get_epitope(protein, epitopeName)
        first = epitope['first']
        last = epitope['last']
        RangeResults.__init__(self, proteinName, first, last,  debug, baseline,
                              displacement, include_wt)
        self.expectedAverageEntropy = epitope['reference_value']
        if len(self.siteResults) != last - first + 1:
            eprint('Warning: only %d sites were found for %s %s'
                   % (len(self.siteResults), proteinName, epitopeName))


class ProteinResults:
    """Compares SE results for an entire protein against the reference
    data from Florencia et al"""

    def __init__(self, proteinName, debug, baseline, displacement, include_wt):
        self.proteinName = proteinName
        protein = get_protein(proteinName)
        self.epitopeResults = [EpitopeResults(
            proteinName, epitopeName,  debug, baseline, displacement,
            include_wt) for epitopeName in protein]
        comparisons = [(er.averageEntropy, er.expectedAverageEntropy)
                       for er in self.epitopeResults
                       if (er.expectedAverageEntropy is not None
                           and er.averageEntropy is not None)]
        if len(comparisons) > 1:
            self.r2 = r_squared(comparisons)
        else:
            self.r2 = None


def get_protein(proteinName):
    try:
        return proteins[proteinName.upper()]
    except KeyError:
        eprint("Invalid protein name: " + proteinName)
        eprint("Valid protein names: " +
               ', '.join(proteins.keys()))
        exit(1)


def get_epitope(protein, epitopeName):
    try:
        return protein[epitopeName.upper()]
    except KeyError:
        eprint("Invalid epitope name: " + epitopeName)
        eprint("Valid epitope names in protein: " +
               ', '.join(protein.keys()))
        exit(1)


def r_squared(pairs):  # pass a list of tuples [ (x1,y1), (x2,y2), … ]
    """r-squared computation used to measure correlation of our
    results with Pereyra et al"""
    num = len(pairs)
    avgX = sum([x for x, y in pairs]) / num
    avgY = sum([y for x, y in pairs]) / num
    residuals = [(x - avgX, y - avgY) for x, y in pairs]
    numerator = sum(rx * ry for rx, ry in residuals)
    denominator = math.sqrt(sum(
        [x * x for x, y in residuals])) * math.sqrt(sum([y * y
                                                         for x, y
                                                         in residuals]))
    if denominator == 0:
        return float('NaN')
    return numerator / denominator


def compute_entropy(probabilities):
    # Compute Shannon Entropy of the given distribution
    # https://en.wikipedia.org/wiki/Entropy_(information_theory)
    # SE = - ∑ P_i * log_b( P_i )
    # where
    #  P_i is the probability of sample i
    #  log_b is log-base-b
    # It is not clear in the Pereyra paper what units were used, so we
    # try the two most common bases: 2 and e
    entropy = 0.0
    for probability in probabilities:
        # base e: entropy measured in nats
        entropy -= probability * bigfloat.log(probability)
    return entropy


def compute_boltzmann(energies, debug):
    # create Boltzmann distribution
    # https://en.wikipedia.org/wiki/Boltzmann_distribution
    # p_i = exp( -E_i / kT ) / ∑ exp( -E_j / kT )
    # where:
    #    p_i is the probability of state i occuring
    #    E_i is the energy of state i
    #    E_j is the energy of state j, where the ∑ in the denominator
    #    iterates over all states j
    # first we calculate exp( -E_i / kT ) for all states
    if len(energies) == 0:
        return []
    if debug:
        print("Boltzmann Distribution")
    divisor = BigFloat.exact(0.0)
    for energy in energies:
        ep = bigfloat.exp(-energy)
        if debug:
            print("energy, ep = %g, %g" % (energy, ep))
        # divisor = ∑ exp( -E_j / kT )
        divisor += ep
    probabilities = []
    for energy in energies:
        # p_i = exp( -E_i / kT ) / divisor
        numerator = bigfloat.exp(-energy)
        probability = numerator / divisor
        # save probability to dictionary
        probabilities.append(float(probability))
        if debug:
            print("%s / %s = %s" % (numerator, divisor, probability))

    return probabilities


if __name__ == '__main__':
    main()
