import numpy as np
from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from equus.backbone.utils import find_carbons_by_degree


def replace_atom(mol: Mol, idx: int, atomic_num: int) -> Mol:
    """
    idx: index of atom to be replaced
    atomic_num: atomic number of the replacement atom

    Note: This function will not deal with failed cases!
    """

    editable_molecule = Chem.RWMol(mol)
    editable_molecule.GetAtomWithIdx(idx).SetAtomicNum(newNum=atomic_num)

    return editable_molecule.GetMol()


def replace_carbon_atom(mol: Mol, idx: int | None, atomic_num: int | None) -> Mol:

    """
    idx: index of carbon atom to be replaced
    atomic_num: atomic number of the replacement atom

    Note: This function will not deal with failed cases!
    """

    # if not specified, picks a random new atom from {N, O, S}
    if atomic_num is None:
        candidates = [7, 8, 16]
        weights = np.array([4, 4, 1])
        weights = weights / weights.sum()
        atomic_num = np.random.choice(candidates, 1, p=weights)[0]

    # if not specified, picks a random carbon atom
    # neutral molecules only -> secondary carbons for O/S, secondary/tertiary carbons for N
    if idx is None:
        carbons = find_carbons_by_degree(mol=mol, degree=2)
        if atomic_num == 7:
            carbons += find_carbons_by_degree(mol=mol, degree=3)
        idx = np.random.choice(carbons, 1)[0]

    return replace_atom(mol=mol, idx=idx, atomic_num=atomic_num)
