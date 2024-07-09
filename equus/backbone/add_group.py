from random import choice

from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from equus.edit.utils import connect_base_mol_and_deuterium, read_smiles

# Carbons


def add_methyl(mol: Mol, idx: int) -> Mol:
    methyl = read_smiles("[2H]C")
    return connect_base_mol_and_deuterium(mol, idx, methyl)


def add_ethyl(mol: Mol, idx: int) -> Mol:
    ethyl = read_smiles("[2H]CC")
    return connect_base_mol_and_deuterium(mol, idx, ethyl)


def add_isopropyl(mol: Mol, idx: int) -> Mol:
    isopropyl = read_smiles("[2H]C(C)C")
    return connect_base_mol_and_deuterium(mol, idx, isopropyl)


# Halogens


def add_halogen(mol: Mol, idx: int, halogen: str) -> Mol:
    assert halogen in ["F", "Cl", "Br", "I"], "Incorrect halogen symbol."
    halogen_D = read_smiles(f"[2H]{halogen}")
    return connect_base_mol_and_deuterium(mol, idx, halogen_D)


# Miscellaneous


def add_group(mol: Mol, idx: int, key: str | None) -> Mol:

    groups = {
        "CF3": "[2H]C(F)(F)F",
        "CCl3": "[2H]C(Cl)(Cl)Cl",
        "CN": "[2H]C#N",
        "NO2": "[2H][N+](=O)[O-]",
        "CHO": "[2H]C=O",
        "COOCH3": "[2H]C(=O)OC",
        "OCH3": "[2H]OC",
        "SCH3": "[2H]C=S",
    }

    if key is None:
        key = choice(groups.keys())

    group = read_smiles(groups[key])
    return connect_base_mol_and_deuterium(mol, idx, group)
