import random
import warnings

from rdkit import Chem

from equus.edit import clean, read_smiles, to_smiles


def break_random_single_bond(smi: str) -> str | None:
    """
    Randomly breaks one single bond.
    Does not allow the generated molecule to have a ring with more than 7 atoms.
    User is responsible for dealing with failed cases.
    """

    mol = read_smiles(smi, no_aromatic_flags=True, hydrogens=False)

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

    # Get the atom indices of the atoms connected by the bond
    atom1_idx = bond_to_break.GetBeginAtomIdx()
    atom2_idx = bond_to_break.GetEndAtomIdx()

    # Break the bond
    em = Chem.EditableMol(mol)
    em.RemoveBond(atom1_idx, atom2_idx)
    mol = em.GetMol()

    # Add hydrogen atoms to the resulting fragments
    mol = clean(Chem.AddHs(mol))
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) > 7:
            warnings.warn("There is a big ring!")
            return None
    return to_smiles(mol)
