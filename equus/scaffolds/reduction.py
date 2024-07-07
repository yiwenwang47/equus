import random

import numpy as np
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from scipy.stats import norm

from equus.edit import read_smiles, to_smiles


def gaussian_random(b: int, mu: float | None = None, sigma: float | None = None) -> int:
    """
    samples an integer using a Gaussian dist within [1, b]
    """
    if mu is None:
        mu = (b + 1) / 2  # default mean is the midpoint of the range
    if sigma is None:
        sigma = (b - 1) / 3  # default standard deviation is 1/3 of the range
    candidates = np.arange(b) + 1
    weights = norm.pdf(candidates, mu, sigma)
    weights /= weights.sum()
    return np.random.choice(candidates, 1, p=weights)[0]


def reduce_bonds(mol: Mol, list_of_bond_idx: list[int]) -> Mol:
    """
    Reduces the bonds indexed by list_of_bond_idx to single bonds
    """
    editable_mol = Chem.RWMol(mol)
    for bond_idx in list_of_bond_idx:
        bond = editable_mol.GetBondWithIdx(bond_idx)
        bond.SetBondType(Chem.rdchem.BondType.SINGLE)
    Chem.SanitizeMol(editable_mol)
    return editable_mol.GetMol()


def randomly_reduce_bonds(smi: str, num_bonds_to_reduce: int | None = None) -> str:
    """
    Randomly selects double/triple bonds to reduce.
    """

    mol = read_smiles(smi, no_aromatic_flags=True, hydrogens=False)

    # Find double and triple bonds in aromatic rings
    bonds_to_reduce = []
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() > 1:
            bonds_to_reduce.append(bond.GetIdx())
    number = len(bonds_to_reduce)

    # sample num_bonds_to_reduce within range of (1, number-1)
    if num_bonds_to_reduce is None:
        num_bonds_to_reduce = gaussian_random(number - 1)

    # Sample bonds to reduce
    selected_bonds = random.sample(bonds_to_reduce, num_bonds_to_reduce)
    new_mol = reduce_bonds(mol, selected_bonds)

    return to_smiles(new_mol)
