#!/usr/bin/env python
# # -*- coding: utf-8 -*-

"""
Functions for processing pdb files
"""

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
