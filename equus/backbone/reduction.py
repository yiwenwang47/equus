import random

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops
from rdkit.Chem.rdchem import Mol
from scipy.stats import norm

from equus.edit import clean, read_smiles, to_smiles


def gaussian_random(b: int, mu: float | None = None, sigma: float | None = None) -> int:
    """
    samples an integer using a Gaussian dist within [1, b]
    """
    assert b >= 1
    if b > 1:
        if mu is None:
            mu = (b + 1) / 2  # default mean is the midpoint of the range
        if sigma is None:
            sigma = (b - 1) / 3  # default standard deviation is 1/3 of the range
        candidates = np.arange(b) + 1
        weights = norm.pdf(candidates, mu, sigma)
        weights /= weights.sum()
        return np.random.choice(candidates, 1, p=weights)[0]
    else:
        return 1


def reduce_one_bond(mol: Mol, bond_idx: int) -> Mol:
    """
    Reduces the bonds indexed by bond_idx to a single bond
    """
    editable_mol = Chem.RWMol(mol)
    bond = editable_mol.GetBondWithIdx(bond_idx)
    bond.SetBondType(Chem.rdchem.BondType.SINGLE)
    new_mol = editable_mol.GetMol()
    Chem.SanitizeMol(new_mol)
    return new_mol


def find_unsaturated_bonds(mol: Mol) -> list[int]:

    # Find double and triple bonds in aromatic rings
    bonds_to_reduce = []
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() > 1:
            bonds_to_reduce.append(bond.GetIdx())

    return bonds_to_reduce


def randomly_reduce_bonds(mol: Mol, num_bonds_to_reduce: int | None = None) -> Mol:
    """
    Randomly selects double/triple bonds to reduce.
    """

    smi = to_smiles(mol)
    mol = read_smiles(
        smi, no_aromatic_flags=False, hydrogens=False
    )  # to maintain some sanity

    num_of_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

    # to preserve aromaticity
    if num_of_aromatic_rings <= 1:
        print(
            f"Only {num_of_aromatic_rings} aromatic ring(s) present. Not going to reduce any bonds."
        )
        return mol

    # sample num_bonds_to_reduce within range of (1, num_of_aromatic_rings-1)
    if num_bonds_to_reduce is None:
        num_bonds_to_reduce = gaussian_random(num_of_aromatic_rings - 1)

    smi = to_smiles(mol)
    mol = read_smiles(
        smi, no_aromatic_flags=True, hydrogens=False
    )  # to maintain some sanity

    for _ in range(num_bonds_to_reduce):
        # Sample bonds to reduce
        candidates = find_unsaturated_bonds(mol=mol)
        idx = random.choice(candidates)
        mol = reduce_one_bond(mol, idx)

    return rdmolops.AddHs(mol)
