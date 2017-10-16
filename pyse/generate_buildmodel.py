#!/usr/bin/python
# coding=utf-8
#
# Generate FoldX job directories for parallel runs.  Each directory contains one WT -> amino replacement for a protein, on all corresponding chains.

from proteins import pdbs, codes, chains, proteins
import argparse
import sys
import uuid
import os


def parse_seqfile(filename):
    """Parse a sequence file generated by FoldX's SequenceOnly command to get the amino sequence found in FoldX's internal representation.
    :param filename: the filename of the sequence file generated by FoldX.
    :returns: dictionary of chains with a list of aminos for that chain as the value. """
    seqfile = open("../" + filename)
    line_number = 0

    chain_aminos = {}

    for line in seqfile:
        line = line.strip().rstrip()

        if line_number >= 11:
            if line != "":
                line_array = line.split("\t")
                if line_array[1] in chain_aminos:
                    chain_aminos[line_array[1]].append(line_array[3])
                else:
                    chain_aminos[line_array[1]] = [line_array[3]]

        line_number = line_number + 1

    return chain_aminos


def makeJobDir(protein_name, line):
    """Generate a job directory for a work item.
    :param protein_name: the name of the protein that this work item is used with.
    :param line: the mutation(s) that this work item should perform. 
    :returns: None.  Generates a folder and FoldX work files on disk. """
    pdb_name = pdbs[protein_name]
    dir_name = "foldxbm-" + str(uuid.uuid4())
    full_loc = protein_name + "/" + dir_name
    os.makedirs(full_loc)
    f = open(full_loc + '/individual_list.txt', 'w+')
    f2 = open(full_loc + '/list.txt', 'w+')
    f2.write(os.path.basename(pdb_name) + "\n")
    f2.close()
    f.write(line + "\n")
    f.close()


def dumpChains(chain_aminos, protein_name):
    """ Helper function for dumping chain amino data in a python-friendly format. """

    print protein_name + ":"
    print "{",

    for group in chains[protein_name]:
        for chain in group:
            sys.stdout.write("'" + chain + "':" + '"')

            for amino in chain_aminos[chain]:
                sys.stdout.write(amino[0])

            print '",',

    print "}"


def genDirs(group, chain_aminos, protein_name):
    """Generates many job directories for every possible WT mutation for a set of chains on a protein.
    :param group: list of chains that are expected to be identical within a given protein.
    :param chain_aminos: dictionary where keys are the chain letters, and values are the amino sequence for that chain as a string. 
    :param protein_name: the name of the protein, such as p24 or REV.
    :returns: None.  Generates many folders and FoldX work files on disk. """

    # select representitive chain.
    rep_chain = group[0]
    print "Selecting representitive chain: " + rep_chain

    print "keys in chain aminos:"
    for key, value in chain_aminos.iteritems():
        print key

    for amino in chain_aminos[rep_chain]:
        # pull apart data.
        amino_val = amino[0]
        chain_val = amino[1]
        pos_val = amino[2:]

        # ignore missing locations.
        if amino_val not in ['x', '_', '-']:
            if amino_val == 'x':
                print "Found a missing or bad amino"
            for code in codes:
                individual_list_string = ""
                for chain in group:
                    individual_list_string = individual_list_string + \
                        amino_val + chain + pos_val + codes[code] + ","

                individual_list_string = individual_list_string[:-1] + ";"
                makeJobDir(protein_name, individual_list_string)


def main():
    parser = argparse.ArgumentParser(
        description='Generate distributable foldx buildmodel dirs')
    parser.add_argument('protein', help='protein name - for example rt or p24')
    args = parser.parse_args()
    protein = args.protein.upper()

    filename = pdbs[protein]
    filename = ("data/foldxseq/SequenceOnly_" +
                filename.split("/")[-1])[0:-4] + ".fxout"

    print "Chain sets for this protein:"
    print chains[protein]

    # dict of chain: list of amino and pos data in format EA32
    print "Search for seqres in: " + filename
    chain_aminos = parse_seqfile(filename)

    # each list in chains[protein] is a set of matching chains in the protein.
    # for p24, for example, this is a set with 6 entries, A through F.
    # for other proteins this is AB match vs CD match.
    # most of the time we will only have one group.
    # select the first element as the representitive chain for generation.
    # tweaking this here would result in less errors, but changes in overall energy.
    rep_chain_number = ord(proteins[protein].values()[0]['chains'][0]) - 65

    print "Selecting Rep chain number: " + str(rep_chain_number)

    genDirs(chains[protein][rep_chain_number], chain_aminos, protein)


if __name__ == "__main__":
    main()