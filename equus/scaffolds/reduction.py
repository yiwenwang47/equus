import random

import numpy as np
from rdkit import Chem
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


def reduce_bonds(smi: str) -> str:

    mol = read_smiles(smi, no_aromatic_flags=True, hydrogens=False)

    # Find double and triple bonds in aromatic rings
    bonds_to_reduce = []
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() > 1:
            bonds_to_reduce.append(bond.GetIdx())
    number = len(bonds_to_reduce)

    # sample num_bonds_to_reduce within range of (1, number-1)
    num_bonds_to_reduce = gaussian_random(number - 1)

    # Randomly select bonds to reduce
    selected_bonds = random.sample(bonds_to_reduce, num_bonds_to_reduce)

    # Reduce the selected bonds
    emol = Chem.RWMol(mol)
    for bond_idx in selected_bonds:
        bond = emol.GetBondWithIdx(bond_idx)
        bond.SetBondType(Chem.rdchem.BondType.SINGLE)

    # Sanitize the molecule to update the valences
    Chem.SanitizeMol(emol)

    return to_smiles(emol.GetMol())
