from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.rdchem import Mol, RWMol

from equus.edit.utils import num_of_Hs


def oxidize_bonds(mol: Mol, list_of_bond_idx: list[int]) -> Mol:
    """
    Oxidizes the bonds indexed by list_of_bond_idx to double bonds
    """
    editable_mol = Chem.RWMol(mol)
    for bond_idx in list_of_bond_idx:
        bond = editable_mol.GetBondWithIdx(bond_idx)
        bond.SetBondType(Chem.rdchem.BondType.DOUBLE)
    Chem.SanitizeMol(editable_mol)
    return editable_mol.GetMol()


def prep_mol(mol: Mol, idx: int, to_remove: int) -> tuple[RWMol, int]:
    """
    Helper function for oxidation reactions.

    idx: index of atom where oxidation happens
    to_remove: number of hydrogens to be removed
    """

    # marks the atom where oxidation happens
    mol = Chem.Mol(mol)
    mol.GetAtomWithIdx(idx).SetAtomMapNum(1)
    mol = rdmolops.AddHs(mol)

    editable_mol = Chem.RWMol(mol)
    atom = editable_mol.GetAtomWithIdx(idx)
    neighbors = atom.GetNeighbors()

    assert num_of_Hs(atom) >= to_remove, "Not enough hydrogen atoms!"

    for neighbor in neighbors:
        if neighbor.GetAtomicNum() == 1:
            editable_mol.RemoveAtom(neighbor.GetIdx())
            to_remove -= 1
            if to_remove == 0:
                break

    # finds the marked atom
    for atom in editable_mol.GetAtoms():
        if atom.GetAtomMapNum() == 1:
            idx = atom.GetIdx()
            atom.SetAtomMapNum(0)
            break

    return editable_mol, idx


def CH2_to_CO(mol: Mol, idx: int) -> Mol:
    """
    Oxidizion of CH2 into CO
    idx: index of C atom where the oxidation happens
    """

    oxygen = Chem.MolFromSmarts("[O]")

    # removes 2 hydrogens
    to_remove = 2
    editable_mol, idx = prep_mol(mol, idx, to_remove)

    # connects C with O
    combo = Chem.RWMol(Chem.CombineMols(editable_mol, oxygen))
    combo.AddBond(idx, editable_mol.GetNumAtoms(), order=Chem.rdchem.BondType.DOUBLE)
    final_mol = combo.GetMol()

    return final_mol


def CH2_to_CCH2(mol: Mol, idx: int) -> Mol:
    """
    Oxidizion of CH2 into C=CH2
    idx: index of C atom where the oxidation happens
    """

    CH2 = Chem.MolFromSmarts("[CH2]")
    for atom in CH2.GetAtoms():
        if atom.GetAtomicNum() == 6:
            idx2 = atom.GetIdx()

    # removes 2 hydrogens
    to_remove = 2
    editable_mol, idx = prep_mol(mol, idx, to_remove)

    # connects C with CH2
    combo = Chem.RWMol(Chem.CombineMols(editable_mol, CH2))
    combo.AddBond(
        idx, idx2 + editable_mol.GetNumAtoms(), order=Chem.rdchem.BondType.DOUBLE
    )
    final_mol = combo.GetMol()

    return final_mol


def CH3_to_CN(mol: Mol, idx: int) -> Mol:
    """
    Oxidizion of -CH3 into -CN
    idx: index of C atom where the oxidation happens
    """

    nitrogen = Chem.MolFromSmarts("[N]")

    # removes 3 hydrogens
    to_remove = 3
    editable_mol, idx = prep_mol(mol, idx, to_remove)

    # connects C with N
    combo = Chem.RWMol(Chem.CombineMols(editable_mol, nitrogen))
    combo.AddBond(idx, editable_mol.GetNumAtoms(), order=Chem.rdchem.BondType.TRIPLE)
    final_mol = combo.GetMol()

    return final_mol


def CH3_to_CCH(mol: Mol, idx: int) -> Mol:
    """
    Oxidizion of -CH3 into -CCH
    idx: index of C atom where the oxidation happens
    """

    CH = Chem.MolFromSmarts("[CH]")
    for atom in CH.GetAtoms():
        if atom.GetAtomicNum() == 6:
            idx2 = atom.GetIdx()

    # removes 3 hydrogens
    to_remove = 3
    editable_mol, idx = prep_mol(mol, idx, to_remove)

    # connects C with CH
    combo = Chem.RWMol(Chem.CombineMols(editable_mol, CH))
    combo.AddBond(
        idx, idx2 + editable_mol.GetNumAtoms(), order=Chem.rdchem.BondType.TRIPLE
    )
    final_mol = combo.GetMol()

    return final_mol
