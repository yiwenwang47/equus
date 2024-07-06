from rdkit import AllChem, Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.Scaffolds import MurckoScaffold

from equus.edit import read_smiles, to_smiles

# copied from https://github.com/rdkit/rdkit/discussions/6844
PATT = Chem.MolFromSmarts("[$([D1]=[*])]")
REPL = Chem.MolFromSmarts("[*]")


def get_scaffold(mol, real_bm=False, use_csk=False, use_bajorath=False):
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
