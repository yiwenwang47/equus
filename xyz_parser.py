import numpy as np

def parse_xyz(xyzname):
    
    """
    Parse the coordinates from an xyz file.
    """

    f = open(xyzname, 'r')
    text = f.read().split('\n')
    natoms = int(text[0])
    text = text[2:2+natoms]
    atoms, coords_all = [], []
    for line in text:
        array = line.split()
        assert len(array) == 4
        atom = array[0]
        coords = np.array([float(array[i]) for i in range(1,4)])
        atoms.append(atom)
        coords_all.append(coords)
    coords_all = np.array(coords_all)
    return atoms, coords_all