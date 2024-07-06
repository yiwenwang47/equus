from rdkit import AllChem, Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

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


# TO DO
# note: no ==, no triple bond in ring, no charges on atoms unless the atoms belong to a -NO2 group

# def validate(smi):
