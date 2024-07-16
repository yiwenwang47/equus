import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.rdchem import Mol

from equus.backbone.utils import find_carbons_by_degree
from equus.diamine.utils import clean


def replace_atom(mol: Mol, idx: int, atomic_num: int) -> Mol:
    """
    idx: index of atom to be replaced
    atomic_num: atomic number of the replacement atom

    Note: This function will not deal with failed cases!
    """

    editable_molecule = Chem.RWMol(mol)
    atom = editable_molecule.GetAtomWithIdx(idx)
    atom.SetAtomicNum(atomic_num)

    # deals with the special case where Sanitize thinks S needs 2 Hs
    if atomic_num == 16:
        atom.SetAtomMapNum(1)

    # to maintain some sanity
    new_mol = rdmolops.RemoveHs(editable_molecule)
    Chem.SanitizeMol(new_mol)
    new_mol = rdmolops.AddHs(new_mol)

    if atomic_num == 16:
        editable_mol = Chem.RWMol(new_mol)
        for atom in editable_mol.GetAtoms():
            if atom.GetAtomMapNum() == 1:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 1:
                        editable_mol.RemoveAtom(neighbor.GetIdx())
                atom.SetAtomMapNum(0)
                break
        new_mol = editable_mol.GetMol()
        new_mol.UpdatePropertyCache()

    return new_mol


def replace_carbon_atom(
    mol: Mol,
    idx: int | None = None,
    atomic_num: int | None = None,
    verbose: bool = False,
) -> Mol:

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
        atomic_num = int(np.random.choice(candidates, 1, p=weights)[0])

    # if not specified, picks a random carbon atom
    # neutral molecules only -> secondary carbons for O/S, secondary/tertiary carbons for N
    if idx is None:
        carbons = find_carbons_by_degree(mol=mol, degree=2)
        if atomic_num == 7:
            carbons = find_carbons_by_degree(
                mol=mol, degree=2
            ) + find_carbons_by_degree(mol=mol, degree=3)

        # the carbon to be replaced should not have strange neighbors (halogens, for instance)
        real_carbons = []
        for i in carbons:
            to_keep = True
            atom = mol.GetAtomWithIdx(i)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() not in [1, 6, 7, 8, 16]:
                    to_keep = False
                    break
            if to_keep:
                real_carbons.append(i)

        if len(real_carbons) == 0:
            if verbose:
                print("Nothing to replace.")
            return mol
        idx = int(np.random.choice(real_carbons, 1)[0])

    return replace_atom(mol=mol, idx=idx, atomic_num=atomic_num)
