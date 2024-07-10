import random
import warnings

from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from equus.edit import clean, read_smiles, to_smiles


def break_single_bond(mol: Mol, bond_idx: int) -> Mol:

    bond_to_break = mol.GetBondWithIdx(bond_idx)

    # Get the atom indices of the atoms connected by the bond
    atom1_idx = bond_to_break.GetBeginAtomIdx()
    atom2_idx = bond_to_break.GetEndAtomIdx()

    # Break the bond
    editable_mol = Chem.EditableMol(mol)
    editable_mol.RemoveBond(atom1_idx, atom2_idx)
    mol = editable_mol.GetMol()

    # Add hydrogen atoms to the resulting fragments
    mol = clean(Chem.AddHs(mol))

    return mol


def break_random_single_bond(mol: Mol) -> Mol | None:
    """
    Randomly breaks one single bond.
    Does not allow the generated molecule to have a ring with more than 7 atoms.
    User is responsible for dealing with failed cases.
    """

    smi = to_smiles(mol)
    mol = read_smiles(
        smi, no_aromatic_flags=True, hydrogens=False
    )  # to maintain some sanity

    # Get all single bonds in the molecule
    single_bonds = [
        bond
        for bond in mol.GetBonds()
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and bond.IsInRing()
    ]

    if not single_bonds:
        raise ValueError("No single bonds found in the molecule.")

    # Randomly select one single bond to break
    bond_to_break = random.choice(single_bonds)
    bond_idx = bond_to_break.GetIdx()
    mol = break_single_bond(mol=mol, bond_idx=bond_idx)

    # big rings are not allowed
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) > 7:
            warnings.warn("There is a big ring!")
            return None

    return mol
