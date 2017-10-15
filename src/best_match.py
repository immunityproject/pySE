#!/usr/bin/python
# coding=utf-8
#
# Searches proteins for specific epitope occurrences and outputs the chain relative location.
# Can fuzzy match partially mutated epitopes, or sections with missing data.

from proteins import all_aminos, chains, chain_bases
from argparse import ArgumentParser


DEBUG = False


class Result:
    def __init__(self, proteinName, chain, offset, sequence, score):
        self.proteinName = proteinName
        self.chain = chain
        self.offset = offset
        self.sequence = sequence
        self.score = score


def getScore(seq1, seq2):
    """Compare two sequences of equal length to determine if they are a match.
    :param seq1: the first amino sequence to compare.
    :param seq2: the second amino sequence to compare.
    :returns: a number representing the similarity between sequences."""

    total_match = 0.0

    if(len(seq1) != len(seq2)):
        print "Error, seq len mismatch:" + seq1 + ":" + seq2
    else:
        for i in range(len(seq1)):
            if (seq1[i - 1] == seq2[i - 1]):
                total_match = total_match + 1

    return total_match / len(seq1)


def scanSingleProteinForMatch(protein, proteinArray, epitope, cutoff):
    """Scan a protein for an epitope.
    :param protein: the name of the protein to scan.
    :param proteinArray: protein/chain aminos for each protein to search.
    :param epitope: the amino sequence to scan for.
    :param cutoff: minimal confidence score to include results.
    :returns: array of Result objects for any qualifying matches"""
    resultArray = []

    for chain in proteinArray[protein]:
        chain_base = chain_bases[protein][chain]
        for startLoc in range(len(proteinArray[protein][chain]) - len(epitope) + 1):
            testEpitope = proteinArray[protein][chain][startLoc:startLoc +
                                                       len(epitope)]
            scoreResult = getScore(epitope, testEpitope)
            if (scoreResult >= cutoff):
                resultArray.append(
                    Result(protein, chain, startLoc + chain_base, testEpitope, scoreResult))

            if(DEBUG):
                print "Added " + protein + "-" + chain + "-" + str(startLoc) + " (" + testEpitope + ") : " + str(scoreResult)

    return resultArray


def scanProteinsForMatch(proteinArray, epitope, cutoff):
    """Scan all proteins for an epitope.
    :param proteinArray: protein/chain aminos for each protein to search.
    :param epitope: the amino sequence to scan for.
    :param cutoff: minimal confidence score to include results.
    :returns: sorted array of Result objects for any qualifying matches"""

    resultsArray = []

    for protein in proteinArray:
        proteinResult = scanSingleProteinForMatch(
            protein, proteinArray, epitope, cutoff)
        resultsArray.extend(proteinResult)

    resultsArray.sort(key=lambda result: 1 - result.score)
    return resultsArray


def printTopResults(resultsArray, eptiope, numberOfResults):
    """Print details for the best matches for a given epitope, ranked by confidence.
    :param resultsArray: list of all result objects from a previous search step, usually from running scanProteinsForMatch
    :param epitope: the amino sequence that was searched for.
    :param numberOfResults: the number of results to show, max.
    :returns: None."""
    print "Best " + str(numberOfResults) + " matches for:" + eptiope

    for i in range(numberOfResults):
        if(i >= len(resultsArray)):
            print str(i) + ": ---"
        else:
            result = resultsArray[i]
            print str(i) + ": " + result.proteinName + "-" + result.chain + "-" + str(result.offset),
            print " (" + result.sequence + ") " + "SCORE:" + str(result.score)


def findChainsFromBestHit(bestResult):
    """ Helper function that returns the full set of chains for a high confidence hit on a single chain on a given protein."""
    for chain in chains[bestResult.proteinName]:
        if bestResult.chain in chain:
            return chain


def constructResultDict(first, last, chain_array, conf, ref):
    """Collect variables relevant to a epitope match to a dictionary.
    :param first: the starting position of the match in the matching chain aminos.
    :param last: the ending position of the match in the matching chain aminos.
    :param chain_array: an array containing the set of chains upon which this match was found.
    :param conf: how confident the match was.
    :param ref: a reference value that we expect, such as a published result.
    :returns: a dictionary containing these values with keys used in other scripts."""
    return_dict = {}
    return_dict['first'] = first
    return_dict['last'] = last
    return_dict['chains'] = chain_array
    return_dict['conf'] = conf
    return_dict['reference_value'] = float(ref)

    return return_dict


def findLocalMatchFromList(epitope_array):
    """Scans all proteins for a list of epitopes.
    :param epitope_array: array of amino sequences that should be seeked.
    :returns: a dictionary of proteins, each protein a dictionary of epitopes that were found there, each with relevant epitope match data.  Can be inserted into proteins.py file to be permanently referenced by other scripts."""
    # initialize protein dict.
    protein_result = {}
    for protein in all_aminos:
        protein_result[protein] = {}

    for epitope_data in epitope_array:
        epitope = epitope_data[0]
        shortname = epitope[0] + epitope[-1] + str(len(epitope))
        resultsArray = scanProteinsForMatch(all_aminos, epitope, 0.2)
        # use best result and look up chain data.
        bestResult = resultsArray[0]
        chain_array = findChainsFromBestHit(bestResult)

        # if the confidence is less than 2/3, assume site is missing from data or too mutated to recognize.
        if bestResult.score < .665:
            print "Bad conf for: " + shortname + " " + str(bestResult.score)

        else:
            protein_result[bestResult.proteinName][shortname] = constructResultDict(
                bestResult.offset, bestResult.offset + len(bestResult.sequence) - 1, chain_array, bestResult.score, epitope_data[1])

    print protein_result
    return protein_result


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument('-e', dest='epitope', default=None,
                        help="list an epitope chain to search for. ex ALFALFQ")
    parser.add_argument('-c', dest='cutoff', default=None,
                        help="cutoff for if an match is interesting")
    parser.add_argument('-f', dest='filename', default=None,
                        help="specify a file to open that contains newline separated list of epitopes to search for")
    args = parser.parse_args()

    if args.epitope is not None:
        args.epitope = args.epitope.upper()

        print "I will search for: " + args.epitope

        print "I have the following proteins loaded into memory:"
        for protein in all_aminos:
            print "\t*" + protein

        resultsArray = scanProteinsForMatch(all_aminos, args.epitope, 0.2)
        printTopResults(resultsArray, args.epitope, 10)

    elif args.filename is not None:
        f = open(args.filename, 'r')
        epitope_array = []
        for line in f:
            split_line = line.strip().split(',')
            epitope_array.append([split_line[0], split_line[1]])

        findLocalMatchFromList(epitope_array)

    else:
        print "Use -e or -f to specify epitope sequence to search for!"
