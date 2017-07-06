#!/usr/bin/python

# This script takes the output of `compute_shannon_entropy.py --scan xx --csv yyy` and creates a Swiss PDB Viewer
# attribute file containing the Structural Entropies for each site.

import argparse, csv

parser = argparse.ArgumentParser()
parser.add_argument('input',help='input site_entropies.csv file which was created by compute_shannon_entropy --scan --csv')
parser.add_argument('output',nargs='?',default='structural_entropy.txt',help='output filename')
args =  parser.parse_args()


with open(args.input) as input_file, open(args.output,'wb') as output_file:
	input = csv.reader(input_file)
	output_file.write('attribute: structuralEntropy\nrecipient: residues\n')
	for protein,site,entropy in input:
		output_file.write('\t:%(site)s.A,%(site)s.C,%(site)s.E\t%(entropy)s\n'%{'site':site,'entropy':entropy})
