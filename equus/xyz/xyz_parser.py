import numpy as np


def lines_to_atoms_coords(lines):
    atoms, coords_all = [], []
    for line in lines:
        array = line.strip().split()
        assert len(array) == 4
        atom = array[0]
        coords = np.array(array[1:]).astype(float)
        atoms.append(atom)
        coords_all.append(coords)
    coords_all = np.array(coords_all)
    return atoms, coords_all


def parse_xyz(xyzname):
    """Parse the coordinates from an xyz file."""

    f = open(xyzname, "r")
    lines = f.read().split("\n")
    natoms = int(lines[0])
    lines = lines[2 : 2 + natoms]
    return lines_to_atoms_coords(lines)
