from __future__ import annotations

# Copied from https://github.com/hjkgrp/molSimplify/blob/master/molSimplify/Classes/globalvars.py
# Dictionary containing atomic mass, atomic number, covalent radius, number of valence electrons
# Data from http://www.webelements.com/ (last accessed May 13th 2015)
elementdict = {
    "X": (1.0, 0, 0.77, 0),
    "H": (1.0079, 1, 0.37, 1),
    "He": (4.002602, 2, 0.46, 2),
    "Li": (6.94, 3, 1.33, 1),
    "Be": (9.0121831, 4, 1.02, 2),
    "B": (10.83, 5, 0.85, 3),
    "C": (12.0107, 6, 0.77, 4),
    "N": (14.0067, 7, 0.75, 5),
    "O": (15.9994, 8, 0.73, 6),
    "F": (18.9984, 9, 0.71, 7),
    "Ne": (20.1797, 10, 0.67, 8),
    "Na": (22.99, 11, 1.55, 1),
    "Mg": (24.30, 12, 1.39, 2),
    "Al": (26.98, 13, 1.26, 3),
    "Si": (28.08, 14, 1.16, 4),
    "P": (30.9738, 15, 1.06, 5),
    "S": (32.065, 16, 1.02, 6),
    "Cl": (35.453, 17, 0.99, 7),
    "Ar": (39.948, 18, 0.96, 8),
    "K": (39.10, 19, 1.96, 1),
    "Ca": (40.08, 20, 1.71, 2),
    "Sc": (44.96, 21, 1.7, 3),
    "Ti": (47.867, 22, 1.36, 4),
    "V": (50.94, 23, 1.22, 5),
    "Cr": (51.9961, 24, 1.27, 6),
    "Mn": (54.938, 25, 1.39, 7),
    "Fe": (55.84526, 26, 1.25, 8),
    "Co": (58.9332, 27, 1.26, 9),
    "Ni": (58.4934, 28, 1.21, 10),
    "Cu": (63.546, 29, 1.38, 11),
    "Zn": (65.39, 30, 1.31, 12),
    "Ga": (69.72, 31, 1.24, 3),
    "Ge": (72.63, 32, 1.21, 4),
    "As": (74.92, 33, 1.21, 5),
    "Se": (78.96, 34, 1.16, 6),
    "Br": (79.904, 35, 1.14, 7),
    "Kr": (83.798, 36, 1.17, 8),
    "Rb": (85.47, 37, 2.10, 1),
    "Sr": (87.62, 38, 1.85, 2),
    "Y": (88.91, 39, 1.63, 3),
    "Zr": (91.22, 40, 1.54, 4),
    "Nb": (92.91, 41, 1.47, 5),
    "Mo": (95.96, 42, 1.38, 6),
    "Tc": (98.9, 43, 1.56, 7),
    "Ru": (101.1, 44, 1.25, 8),
    "Rh": (102.9, 45, 1.25, 9),
    "Pd": (106.4, 46, 1.20, 10),
    "Ag": (107.9, 47, 1.28, 11),
    "Cd": (112.4, 48, 1.48, 12),
    "In": (111.818, 49, 1.42, 3),
    "Sn": (118.710, 50, 1.40, 4),
    "Sb": (121.760, 51, 1.40, 5),
    "Te": (127.60, 52, 1.99, 6),
    "I": (126.90447, 53, 1.40, 7),
    "Xe": (131.293, 54, 1.31, 8),
    "Cs": (132.9055, 55, 2.32, 1),
    "Ba": (137.327, 56, 1.96, 2),
    "La": (138.9, 57, 1.69, 3),
    "Ce": (140.116, 58, 1.63, 4),
    "Pr": (140.90766, 59, 1.76, 5),
    "Nd": (144.242, 60, 1.74, 6),
    "Pm": (145, 61, 1.73, 7),
    "Sm": (150.36, 62, 1.72, 8),
    "Eu": (151.964, 63, 1.68, 9),
    "Gd": (157.25, 64, 1.69, 10),
    "Tb": (158.92535, 65, 1.68, 11),
    "Dy": (162.500, 66, 1.67, 12),
    "Ho": (164.93033, 67, 1.66, 13),
    "Er": (167.259, 68, 1.65, 14),
    "Tm": (168.93422, 69, 1.64, 15),
    "Yb": (173.045, 70, 1.70, 16),
    "Lu": (174.9668, 71, 1.62, 3),
    "Hf": (178.5, 72, 1.50, 8),
    "Ta": (180.9, 73, 1.38, 5),
    "W": (183.8, 74, 1.46, 6),
    "Re": (186.2, 75, 1.59, 7),
    "Os": (190.2, 76, 1.28, 8),
    "Ir": (192.2, 77, 1.37, 9),
    "Pt": (195.1, 78, 1.23, 10),
    "Au": (197.0, 79, 1.24, 11),
    "Hg": (200.6, 80, 1.49, 2),
    "Tl": (204.38, 81, 1.44, 3),
    "Pb": (207.2, 82, 1.44, 4),
    "Bi": (208.9804, 83, 1.51, 5),
    "Po": (208.98, 84, 1.90, 6),
    "At": (209.99, 85, 2.00, 7),
    "Rn": (222.6, 86, 142, 4),
    "Fr": (223.02, 87, 3.48, 8),
    "Ra": (226.03, 88, 2.01, 2),
    "Ac": (277, 89, 1.86, 3),
    "Th": (232.0377, 90, 1.75, 4),
    "Pa": (231.04, 91, 2.00, 5),
    "U": (238.02891, 92, 1.70, 6),
    "Np": (237.05, 93, 1.90, 7),
    "Pu": (244.06, 94, 1.75, 8),
    "Am": (243.06, 95, 1.80, 9),
    "Cm": (247.07, 96, 1.69, 10),
    "Bk": (247.07, 97, 1.68, 11),
    "Cf": (251.08, 98, 1.68, 12),
}

metalslist = [
    "Li",
    "li",
    "LI",
    "lithium",
    "Be",
    "be",
    "BE",
    "beryllium",
    "Na",
    "na",
    "NA",
    "sodium",
    "Mg",
    "mg",
    "MG",
    "magnesium",
    "Al",
    "al",
    "AL",
    "aluminum",
    "aluminium",
    "K",
    "k",
    "potassium",
    "Ca",
    "ca",
    "CA",
    "calcium",
    "Rb",
    "rb",
    "RB",
    "rubidium",
    "Sr",
    "sr",
    "SR",
    "strontium",
    "Cs",
    "cs",
    "CS",
    "cesium",
    "Ba",
    "ba",
    "BA",
    "barium",
    "Fr",
    "fr",
    "FR",
    "francium",
    "Ra",
    "ra",
    "RA",
    "radium",
    "Sc",
    "sc",
    "SC",
    "scandium",
    "Ti",
    "ti",
    "TI",
    "titanium",
    "V",
    "v",
    "vanadium",
    "Cr",
    "cr",
    "CR",
    "chromium",
    "Mn",
    "mn",
    "MN",
    "manganese",
    "Fe",
    "fe",
    "FE",
    "iron",
    "Co",
    "co",
    "CO",
    "cobalt",
    "Ni",
    "ni",
    "NI",
    "nickel",
    "Cu",
    "cu",
    "CU",
    "copper",
    "Zn",
    "zn",
    "ZN",
    "zinc",
    "Ga",
    "ga",
    "GA",
    "gallium",
    "Y",
    "y",
    "yttrium",
    "Zr",
    "zr",
    "ZR",
    "zirconium",
    "Nb",
    "nb",
    "NB",
    "niobium",
    "Mo",
    "mo",
    "MO",
    "molybdenum",
    "Tc",
    "tc",
    "TC",
    "technetium",
    "Ru",
    "ru",
    "RU",
    "ruthenium",
    "Rh",
    "rh",
    "RH",
    "rhodium",
    "Pd",
    "pd",
    "PD",
    "palladium",
    "Ag",
    "ag",
    "AG",
    "silver",
    "Cd",
    "cd",
    "CD",
    "cadmium",
    "In",
    "in",
    "IN",
    "indium",
    "Sn",
    "sn",
    "SN",
    "tin",
    "Hf",
    "hf",
    "HF",
    "hafnium",
    "Ta",
    "ta",
    "TA",
    "tantalum",
    "W",
    "w",
    "tungsten",
    "Re",
    "re",
    "RE",
    "rhenium",
    "Os",
    "os",
    "OS",
    "osmium",
    "Ir",
    "ir",
    "IR",
    "iridium",
    "Pt",
    "pt",
    "PT",
    "platinum",
    "Au",
    "au",
    "AU",
    "gold",
    "Hg",
    "hg",
    "HG",
    "mercury",
    "X",
    "Tl",
    "tl",
    "TL",
    "thallium",
    "Pb",
    "pb",
    "PB",
    "lead",
    "Bi",
    "bi",
    "BI",
    "bismuth",
    "Po",
    "po",
    "PO",
    "polonium",
    "La",
    "la",
    "LA",
    "lanthanum",
    "Ce",
    "ce",
    "CE",
    "cerium",
    "Pr",
    "pr",
    "PR",
    "praseodymium",
    "Nd",
    "nd",
    "ND",
    "neodymium",
    "Pm",
    "pm",
    "PM",
    "promethium",
    "Sm",
    "sm",
    "SM",
    "samarium",
    "Eu",
    "eu",
    "EU",
    "europium",
    "Gd",
    "gd",
    "GD",
    "gadolinium",
    "Tb",
    "tb",
    "TB",
    "terbium",
    "Dy",
    "dy",
    "DY",
    "dysprosium",
    "Ho",
    "ho",
    "HO",
    "holmium",
    "Er",
    "er",
    "ER",
    "erbium",
    "Tm",
    "tm",
    "TM",
    "thulium",
    "Yb",
    "yb",
    "YB",
    "ytterbium",
    "Lu",
    "lu",
    "LU",
    "lutetium",
    "Ac",
    "ac",
    "AC",
    "actinium",
    "Th",
    "th",
    "TH",
    "thorium",
    "Pa",
    "pa",
    "PA",
    "proactinium",
    "U",
    "u",
    "uranium",
    "Np",
    "np",
    "NP",
    "neptunium",
    "Pu",
    "pu",
    "PU",
    "plutonium",
    "Am",
    "am",
    "AM",
    "americium",
    "Cu",
    "cu",
    "CU",
    "curium",
    "Bk",
    "bk",
    "BK",
    "berkelium",
    "Cf",
    "cf",
    "CF",
    "californium",
    "Es",
    "es",
    "ES",
    "einsteinium",
    "Fm",
    "fm",
    "FM",
    "fermium",
    "Md",
    "md",
    "MD",
    "mendelevium",
    "No",
    "no",
    "NO",
    "nobelium",
    "Lr",
    "lr",
    "LR",
    "lawrencium",
]


def ismetal(atom: str) -> bool:
    return atom in metalslist


covalent_radius = lambda atom: elementdict[atom][2]

import copy
from collections import defaultdict
from dataclasses import dataclass

import numpy as np
from scipy.spatial import distance_matrix

from equus.xyz.xyz_parser import parse_xyz


@dataclass
class simple_mol:
    """
    A very simple molecule object.
    Parses information from a xyz file. Now also deals with generic
    cases.
    """

    atoms: list[str]
    coords: np.ndarray

    @staticmethod
    def from_xyzfile(xyzname: str) -> simple_mol:
        atoms, coords = parse_xyz(xyzname)
        return simple_mol(atoms, coords)

    def __post_init__(self):
        self.natoms = len(self.atoms)
        self.graph = get_graph_full_scope(self)
        self.init_distances()
        self.get_all_distances()

    def copy(self) -> simple_mol:
        return copy.deepcopy(self)

    def init_distances(self):
        self.distances = self.graph.copy()

    def get_all_distances(self, depth: int = 100):
        """
        Only calculates shortest-path distances (topological) up to depth.
        """

        for ind in range(self.natoms):
            bfs_distances(self, ind, depth)

    def get_bonded_atoms(self, atom_index: int) -> np.ndarray:
        con = self.graph[atom_index]
        return np.where(con == 1)[0]

    def get_bonded_atoms_multiple(self, atom_ind: list[int]) -> list[int]:
        indices = set(atom_ind)
        bonded = set()
        for i in indices:
            bonded.update(self.get_bonded_atoms(i))
        return list(bonded)


# Connectivities and distances.
# Following the conventional approach of setting bond length cutoffs.

_get = lambda l, indices: [l[i] for i in indices]


def _connect(graph, i, j):
    graph[i][j] = 1
    graph[j][i] = 1


def get_bond_cutoff(a1: str, a2: str) -> float:
    a1, a2 = sorted([a1, a2])
    r1, r2 = covalent_radius(a1), covalent_radius(a2)
    cutoff = 1.15 * (r1 + r2)
    if a1 == "C" and a2 != "H":
        cutoff = min(2.75, cutoff)
    if a1 == "C" and a2 == "Cl":
        cutoff = 2.1
    if a1 == "H" and ismetal(a2):
        cutoff = 1.1 * (r1 + r2)
    if a2 == "H" and ismetal(a1):
        cutoff = 1.1 * (r1 + r2)
    # Strict cutoff for Iodine
    if a1 == "I" and a2 == "I":
        cutoff = 3
    if a1 == "Ir" and a2 == "N":
        cutoff = 2.5
    return cutoff


def get_cutoffs(atoms: list) -> defaultdict[dict]:
    """
    Calculate bond length cutoffs for all possible atom pairs.
    """

    unique = list(set(atoms))
    cutoffs = defaultdict(dict)
    n = len(unique)
    for i in range(n):
        for j in range(i, n):
            a1, a2 = unique[i], unique[j]
            cutoff = get_bond_cutoff(a1, a2)
            cutoffs[a1][a2] = cutoff
            cutoffs[a2][a1] = cutoff
    return cutoffs


def get_graph_full_scope(mol: simple_mol) -> np.ndarray:
    """
    Creates adjacency matrix.
    """

    n = mol.natoms
    graph = np.zeros((n, n))
    cutoffs = get_cutoffs(mol.atoms)

    mol.matrix = distance_matrix(mol.coords, mol.coords)

    for i in range(n):
        for j in range(i + 1, n):
            atom_i, atom_j = mol.atoms[i], mol.atoms[j]
            if mol.matrix[i][j] <= cutoffs[atom_i][atom_j]:
                _connect(graph, i, j)
    return graph


def bfs_distances(mol: simple_mol, origin: int, depth: int):
    """
    A breadth-first search algorithm to find the shortest-path distances
    between any two atoms.
    Only searches for distances up to the given depth.
    """

    all_active = set([origin])
    current_active = set([origin])

    for distance in range(1, depth + 1):
        new_active = set(mol.get_bonded_atoms_multiple(list(current_active)))
        new_active -= all_active
        if not new_active:
            break
        if distance > 1:
            for atom in new_active:
                mol.distances[origin][atom] = distance
        all_active.update(new_active)
        current_active = new_active
