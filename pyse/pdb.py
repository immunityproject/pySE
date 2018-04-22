# -*- coding: utf-8 -*-

"""
Functions for processing pdb files
"""

from pyse.proteins import codes

def parse_pdb(pdb):
    """ Read the x,y,z positions for each atom in a PDB file """
    positions = []
    lctr = 0

    def stripstr(v):
        return str(v).strip()
    for line in pdb:
        lctr += 1

        if not line.startswith('ATOM'):
            continue

        # these are specified as inclusive ranges of 1-based
        # indexes to match the PDB spec
        fieldspec = [
            ('atom', 13, 16, stripstr),
            ('remnant', 18, 20, stripstr),
            ('chain', 22, 22, stripstr),
            ('position', 23, 27, stripstr),
            ('x', 31, 38, float),
            ('y', 39, 46, float),
            ('z', 47, 54, float),
        ]
        # The following fixes an inconsistency in the PDB spec when position
        # is over 999, because x shifts to the right
        if int(line[22:26]) > 999:
            fieldspec[4] = ('x', 32, 39, float)
            fieldspec[5] = ('y', 40, 47, float)
            fieldspec[6] = ('z', 48, 55, float)
        try:
            positions.append({f[0]: f[3](line[f[1] - 1:f[2]])
                              for f in fieldspec})
        except ValueError as e:
            print('Error on line {}: {}'.format(lctr, line))
            print(e)
    return positions

def get_peptide_chains(pdbdb):
    startsite = None
    chain = None
    peptide = ''
    endsite = None
    prevsite = None
    peptidechains = []
    for pdbentry in pdbdb:
        if pdbentry['atom'] != 'N':
            continue

        # Reset Case: The chain changes
        if chain != pdbentry['chain']:
            if chain:
                peptidechains.append({
                    'startsite': startsite,
                    'endsite': endsite,
                    'chain': chain,
                    'peptide': peptide
                })
            startsite = pdbentry['position']
            chain = pdbentry['chain']
            peptide = ''
            endsite = None
            prevsite = None

        peptide += codes[pdbentry['remnant']]
        try:
            endsite = int(pdbentry['position'])
            if prevsite != None and endsite - 1 != prevsite:
                print('Missing sites between {} and {}'.format(prevsite, endsite))
                for i in range(prevsite + 1, endsite + 1):
                    peptide += '-'
            prevsite = int(pdbentry['position'])
        except ValueError:
            print('Skipping {} because it is not an int'.format(
                pdbentry['position']))

    peptidechains.append({
        'startsite': startsite,
        'endsite': endsite,
        'chain': chain,
        'peptide': peptide
    })

    return peptidechains
