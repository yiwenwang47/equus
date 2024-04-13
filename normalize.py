from rdkit import Chem
from rdkit.Chem import rdmolops
from .edit import clean, find_atom_indices, find_naked_atom_idx
from .edit import num_of_Hs, remove_atom, remove_unconnected_Hs

'''
Note: all molecules should have explicit Hydrogens!
'''

mol_H = Chem.MolFromSmiles('[H]', sanitize=False)

def add_one_H(mol_naked):
    '''
    mol_naked: molecule with only one naked atom
    '''
    i = find_naked_atom_idx(mol_naked)[0]
    combo = Chem.EditableMol(Chem.CombineMols(mol_naked, mol_H))
    combo.AddBond(i, mol_naked.GetNumAtoms(), order=Chem.rdchem.BondType.SINGLE)
    mol = clean(combo.GetMol())
    return mol

def remove_OH(mol):

    '''
    mol: molecule with explicit Hs
    replaces all -OH groups with Hydrogens
    '''

    oxygens = find_atom_indices(mol, atomic_number=8)
    OH = [i for i in oxygens if num_of_Hs(mol.GetAtomWithIdx(i))==1]
    for O in OH:
        mol = remove_atom(mol, O)
        mol = remove_unconnected_Hs(mol)
        mol = add_one_H(mol)

    return mol

def remove_SH(mol):

    '''
    mol: molecule with explicit Hs
    replaces all -SH groups with Hydrogens
    '''

    sulfurs = find_atom_indices(mol, atomic_number=16)
    SH = [i for i in sulfurs if num_of_Hs(mol.GetAtomWithIdx(i))==1]
    for S in SH:
        mol = remove_atom(mol, S)
        mol = remove_unconnected_Hs(mol)
        mol = add_one_H(mol)

    return mol

def remove_unwanted_NH2(mol, N_idx):

    '''
    mol: molecule with explicit Hs
    N_idx: idx of N atom in the -NH2 to be removed
    '''

    mol = remove_atom(mol, N_idx)
    mol = remove_unconnected_Hs(mol)
    mol = add_one_H(mol)

    return mol