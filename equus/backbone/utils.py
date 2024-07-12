from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.Scaffolds import MurckoScaffold

from equus.edit import read_smiles, to_smiles


def find_carbons_by_num_of_Hs(
    mol: Mol, min_num_of_Hs: int | None = None, max_num_of_Hs: int | None = None
) -> list[int]:

    # defines the range of number of Hs allowed
    if min_num_of_Hs is None and max_num_of_Hs is None:
        min_num_of_Hs, max_num_of_Hs = 0, 3
    else:
        if min_num_of_Hs is None:
            min_num_of_Hs = 0
        if max_num_of_Hs is None:
            max_num_of_Hs = 3
    assert min_num_of_Hs in [0, 1, 2, 3], "Incorrect min_num_of_Hs!"
    assert max_num_of_Hs in [0, 1, 2, 3], "Incorrect max_num_of_Hs!"
    assert max_num_of_Hs >= min_num_of_Hs, "Incorrect (min_num_of_Hs, max_num_of_Hs)!"

    # defines the smarts expression
    smarts = f"[#6H{min_num_of_Hs}"
    if max_num_of_Hs > min_num_of_Hs:
        for i in range(min_num_of_Hs + 1, max_num_of_Hs + 1):
            smarts += f",#6H{i}"
    smarts += "]"

    # finds the carbons
    pattern = Chem.MolFromSmarts(smarts)
    matches = mol.GetSubstructMatches(pattern)
    return [i[0] for i in matches]


def find_carbons_by_degree(mol: Mol, degree: int) -> list[int]:
    """
    Finds carbon indices.
    Example: primary carbons, degree=1, pattern=[#6H3]
    """
    assert degree in [1, 2, 3, 4]
    num_of_Hs = 4 - degree
    return find_carbons_by_num_of_Hs(
        mol=mol, min_num_of_Hs=num_of_Hs, max_num_of_Hs=num_of_Hs
    )


def find_carbon_carbon_single_bonds(
    mol: Mol, min_num_of_Hs: int = 1, filter_for_oxidation: bool = True
) -> list[int]:
    """
    min_num_of_Hs: minimum number of hydrogen atoms each carbon atom should have
    filter_for_oxidation: if True, carbon atoms connected by the single bonds should not have any pi bonds

    returns a list of bond indices
    """

    # defines the smarts pattern to match C-C single bonds
    sub_smarts = f"[#6H{min_num_of_Hs}"
    for i in range(min_num_of_Hs + 1, 4):
        sub_smarts += f",#6H{i}"
    sub_smarts += "]"
    smarts = f"{sub_smarts}-{sub_smarts}"
    pattern = Chem.MolFromSmarts(SMARTS=smarts)
    matches = mol.GetSubstructMatches(pattern)

    # finds the bond indices
    list_of_bond_idx = []
    for match in matches:
        idx_1, idx_2 = match
        to_keep = True

        # if the bond is to be oxidized, the carbons should not have pi bonds
        if filter_for_oxidation:
            C1, C2 = mol.GetAtomWithIdx(idx_1), mol.GetAtomWithIdx(idx_2)
            for bond in C1.GetBonds():
                if bond.GetBondTypeAsDouble() > 1:
                    to_keep = False
                    break
            if to_keep:
                for bond in C2.GetBonds():
                    if bond.GetBondTypeAsDouble() > 1:
                        to_keep = False
                        break
        if to_keep:
            bond = mol.GetBondBetweenAtoms(idx_1, idx_2)
            list_of_bond_idx.append(bond.GetIdx())

    return list_of_bond_idx


# copied from https://github.com/rdkit/rdkit/discussions/6844
PATT = Chem.MolFromSmarts("[$([D1]=[*])]")
REPL = Chem.MolFromSmarts("[*]")


def get_scaffold(
    mol: Mol, real_bm: bool = False, use_csk: bool = False, use_bajorath: bool = False
) -> Mol:
    Chem.RemoveStereochemistry(mol)  # important for canonization of CSK!
    scaff = MurckoScaffold.GetScaffoldForMol(mol)
    if use_bajorath:
        scaff = AllChem.DeleteSubstructs(scaff, PATT)
    if real_bm:
        scaff = AllChem.ReplaceSubstructs(scaff, PATT, REPL, replaceAll=True)[0]
    if use_csk:
        scaff = MurckoScaffold.MakeScaffoldGeneric(scaff)
        if real_bm:
            scaff = MurckoScaffold.GetScaffoldForMol(scaff)
    return scaff


def validate(smi: str) -> bool:
    """
    To make sure the generated molecule is reasonable.
    Written for sf_tokens.permute_sf_tokens
    """

    mol = read_smiles(smi, no_aromatic_flags=True, hydrogens=True)

    # no fragments
    if len(rdmolops.GetMolFrags(mol, asMols=True)) > 1:
        return False

    # no triple bonds in rings
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() == 3:
            if bond.IsInRing():
                return False

    # no X=X=X sub structure
    smarts = "*=*=*"
    pattern = Chem.MolFromSmarts(smarts)
    if any(mol.GetSubstructMatches(pattern)):
        return False
    smarts = "*=*#*"
    pattern = Chem.MolFromSmarts(smarts)
    if any(mol.GetSubstructMatches(pattern)):
        return False

    # no radicals
    for atom in mol.GetAtoms():
        if atom.GetNumRadicalElectrons() > 0:
            return False

    # no charges on atoms unless the atoms belong to a -NO2 group
    nitro_smarts = "[N+](=O)[O-]"
    nitro_pattern = Chem.MolFromSmarts(nitro_smarts)
    matches = mol.GetSubstructMatches(nitro_pattern)
    nitro_atom_idx = [i for match in matches for i in match]
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0:
            idx = atom.GetIdx()
            if idx not in nitro_atom_idx:
                return False

    return True
