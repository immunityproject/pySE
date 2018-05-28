from pyse.pdb import parse_pdb,get_peptide_chains

def test_get_peptide_chains():
    test_data = [{'chain': 'B', 'atom': 'N', 'remnant': 'LYS',
                  'x': -49.836, 'y': 78.687, 'position': '2', 'z': 6.238},
                 {'chain': 'B', 'atom': 'N', 'remnant': 'ILE',
                  'x': -50.865, 'y': 78.249, 'position': '3', 'z': 5.299},
                 {'chain': 'B', 'atom': 'N', 'remnant': 'PHE',
                  'x': -50.534, 'y': 78.611, 'position': '4', 'z': 3.838},
                 {'chain': 'B', 'atom': 'O', 'remnant': 'LYS',
                  'x': -50.243, 'y': 79.77, 'position': '5', 'z': 3.527},
                 # chains
                 {'chain': 'D', 'atom': 'N', 'remnant': 'PHE',
                  'x': -29.596, 'y': 40.891, 'position': '2', 'z': 20.714},
                 {'chain': 'D', 'atom': 'N', 'remnant': 'PHE',
                  'x': -30.421, 'y': 42.06, 'position': '2A', 'z': 20.397},
                 {'chain': 'D', 'atom': 'N', 'remnant': 'PHE',
                  'x': -30.371, 'y': 42.343, 'position': '12', 'z': 18.883},
                 {'chain': 'D', 'atom': 'N', 'remnant': 'PHE',
                  'x': -29.787, 'y': 43.336, 'position': '13', 'z': 18.412},
                 {'chain': 'D', 'atom': 'N', 'remnant': 'PHE',
                  'x': -30.001, 'y': 43.276, 'position': '14', 'z': 21.229},
                 {'chain': 'D', 'atom': 'N', 'remnant': 'PHE',
                  'x': -29.939, 'y': 43.099, 'position': '15', 'z': 22.73},
                 {'chain': 'D', 'atom': 'N', 'remnant': 'PHE',
                  'x': -31.074, 'y': 43.269, 'position': '16', 'z': 23.513},
                 {'chain': 'D', 'atom': 'N', 'remnant': 'PHE',
                  'x': -28.725, 'y': 42.878, 'position': '17', 'z': 23.371},
                 {'chain': 'D', 'atom': 'N', 'remnant': 'ILE',
                  'x': -30.999, 'y': 43.178, 'position': '18', 'z': 24.91}]
    expected_output = [{'peptide': '-KIF', 'startsite': 1, 'endsite': 4, 'chain': 'B'},
                       {'peptide': '-F---------FFFFFFI', 'startsite': 1, 'endsite': 18, 'chain': 'D'}]
    actual_output = get_peptide_chains(test_data)
    for i in range(0,len(expected_output)):
        for k in expected_output[i].keys():
            assert expected_output[i][k] == actual_output[i][k]
    
