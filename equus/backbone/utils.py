from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.Scaffolds import MurckoScaffold

from equus.diamine import read_smiles


def find_carbons_by_num_of_Hs(
    mol: Mol, min_num_of_Hs: int | None = None, max_num_of_Hs: int | None = None
) -> list[int]:

    # defines the range of number of Hs allowed
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


def validate(smi: str, verbose: bool = False) -> bool:
    """
    To make sure the generated molecule is reasonable.
    Written for sf_tokens.permute_sf_tokens
    """

    mol = read_smiles(smi, no_aromatic_flags=True, hydrogens=True)
    mol_aromatic = read_smiles(smi, no_aromatic_flags=False, hydrogens=True)

    # no fragments
    if len(rdmolops.GetMolFrags(mol, asMols=True)) > 1:
        if verbose:
            print("Multiple fragments!")
        return False

    # no triple bonds in rings
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() == 3:
            if bond.IsInRing():
                if verbose:
                    print("Found a triple bond in a ring!")
                return False

    # no X=X=X sub structure
    smarts = "*=*=*"
    pattern = Chem.MolFromSmarts(smarts)
    if any(mol.GetSubstructMatches(pattern)):
        if verbose:
            print("Consecutive double/triple bonds!")
        return False
    smarts = "*=*#*"
    pattern = Chem.MolFromSmarts(smarts)
    if any(mol.GetSubstructMatches(pattern)):
        if verbose:
            print("Consecutive double/triple bonds!")
        return False

    # no O-C-N or S-C-N sub structure
    pattern1 = Chem.MolFromSmarts("[#7,#16]-C-N")
    pattern2 = Chem.MolFromSmarts("[#7,#16]-C=N")
    pattern3 = Chem.MolFromSmarts("[#7,#16]-C#N")
    if (
        any(mol.GetSubstructMatches(pattern3))
        or any(mol_aromatic.GetSubstructMatches(pattern1))
        or any(mol_aromatic.GetSubstructMatches(pattern2))
    ):
        if verbose:
            print("O/S connected to CN!")
        return False

    # no C=N-H sub structure
    smarts = "[#6]=[#7H1]"
    pattern = Chem.MolFromSmarts(smarts)
    if any(mol.GetSubstructMatches(pattern)):
        if verbose:
            print("Found a C=N-H sub structure!")
        return False

    # no H-N-N-H
    smarts = "[#7H1,#7H2]~[#7H1,#7H2]"
    pattern = Chem.MolFromSmarts(smarts)
    if any(mol.GetSubstructMatches(pattern)):
        if verbose:
            print("Found a H-N-N-H sub structure!")
        return False

    # no B=C bonds
    smarts = "[#5]=[#6]"
    pattern = Chem.MolFromSmarts(smarts)
    if any(mol_aromatic.GetSubstructMatches(pattern)):
        if verbose:
            print("Found a B=C double bond!")
        return False

    # no 3-hetero-atoms-in-a-row sub structure
    smarts = "[#7,#8,#16,#9,#17,#35,#53]~[#7,#8,#16,#9,#17,#35,#53]~[#7,#8,#16,#9,#17,#35,#53]"  # N, O, S, F, Cl, Br, I
    pattern = Chem.MolFromSmarts(smarts)
    if any(mol.GetSubstructMatches(pattern)):
        if verbose:
            print("Found 3 hetero atoms in a row!")
        return False

    # no quaternary carbon (4 single bonds) connected to a hetero atom
    mol_implicit_H = read_smiles(smi, no_aromatic_flags=True, hydrogens=False)
    carbons = find_carbons_by_degree(mol=mol_implicit_H, degree=4)
    for carbon in carbons:
        atom = mol.GetAtomWithIdx(carbon)
        if atom.GetDegree() == 4:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() not in [1, 6]:
                    if verbose:
                        print("Found a quaternary carbon connected to a hetero atom!")
                    return False

    # no radicals
    for atom in mol.GetAtoms():
        if atom.GetNumRadicalElectrons() > 0:
            if verbose:
                print("Found a radical!")
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
                if verbose:
                    print("Found a formal charge!")
                return False

    return True
