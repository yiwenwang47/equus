from rdkit.Chem.rdchem import Atom, Bond, Mol

from equus.edit.utils import connect_base_mol_and_deuterium, num_of_Hs, read_smiles


def aliphatic_ring_double_bonds(smi: str) -> list[Bond]:

    """
    Finds all non-aromatic double bonds in rings.
    Returns a list of Bond objects.
    """

    mol = read_smiles(smi, hydrogens=True)
    candidates = []
    bonds = mol.GetBonds()
    for bond in bonds:
        if (
            bond.GetBondTypeAsDouble() == 2
            and bond.IsInRing()
            and not bond.GetIsAromatic()
        ):
            if num_of_Hs(bond.GetBeginAtom()) > 0 and num_of_Hs(bond.GetEndAtom()) > 0:
                candidates.append(bond)
    return candidates


def aromatic_ring_double_bonds(smi: str) -> list[Bond]:

    """
    Finds all aromatic double bonds in rings.
    Returns a list of Bond objects.
    """

    mol = read_smiles(smi, no_aromatic_flags=False, hydrogens=True)
    candidates = []
    bonds = mol.GetBonds()
    for bond in bonds:
        if bond.GetBondTypeAsDouble() > 1 and bond.IsInRing() and bond.GetIsAromatic():
            if num_of_Hs(bond.GetBeginAtom()) > 0 and num_of_Hs(bond.GetEndAtom()) > 0:
                candidates.append(bond)
    return candidates


__mol_NH2 = read_smiles("N", hydrogens=True)
__N_idx = next(
    atom for atom in __mol_NH2.GetAtoms() if atom.GetAtomicNum() == 7
).GetIdx()


def add_two_NH2_to_bond(bond: Bond) -> Mol:

    """
    Adds two -NH2 groups to the carbons of bond.
    """

    atom = bond.GetBeginAtom()
    H = next(H for H in atom.GetNeighbors() if H.GetAtomicNum() == 1)
    H.SetIsotope(2)

    atom = bond.GetEndAtom()
    H = next(H for H in atom.GetNeighbors() if H.GetAtomicNum() == 1)
    H.SetIsotope(2)

    mol = bond.GetOwningMol()
    mol = connect_base_mol_and_deuterium(__mol_NH2, __N_idx, mol)
    mol = connect_base_mol_and_deuterium(__mol_NH2, __N_idx, mol)

    return mol
