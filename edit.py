from rdkit import Chem
from rdkit.Chem import rdmolops, Descriptors, rdDetermineBonds
from .iso import *
import copy

'''
All editing functions are based on rdkit utilities.
'''

# helper functions

def clean(mol):
    smiles = Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    return mol

def create_smiles(mol):
    return Chem.CanonSmiles(Chem.MolToSmiles(mol))

def embed(mol):
    '''
    generates one conformer
    very helpful because connectivities are well established this way
    '''
    if mol.GetNumConformers()==0:
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)

def read_smiles(smiles_string, no_aromatic_flags=True, hydrogens=True):
    mol = Chem.MolFromSmiles(smiles_string)
    if no_aromatic_flags:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    if hydrogens:
        mol = rdmolops.AddHs(mol=mol)
    return mol

def add_Hs(mol):
    if sum(atom.GetImplicitValence() for atom in mol.GetAtoms()) > 0:
        return clean(rdmolops.AddHs(mol))
    else:
        return mol

def find_atom_indices(mol, atomic_number):
    idx_list = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum()==atomic_number]
    return idx_list

def num_of_Hs(atom):
    explicit = len([H for H in atom.GetNeighbors() if H.GetAtomicNum()==1])
    if atom.GetImplicitValence()>0:
        return explicit + atom.GetNumImplicitHs()
    else:
        return explicit

def find_carbons(mol, num_of_amines=1):
    '''
    finds idx of all carbons as candidates for substitution
    '''
    skip = 1
    carbons = find_atom_indices(mol, 6)
    carbons = [carbon for carbon in carbons if num_of_Hs(mol.GetAtomWithIdx(carbon))>0]
    if len(carbons)==1:
        skip = 0
    nitrogens = find_atom_indices(mol, 7)
    primary_amine_nitrogens = []
    for nitrogen in nitrogens:
        atom = mol.GetAtomWithIdx(nitrogen)
        if num_of_Hs(atom)==2:
            primary_amine_nitrogens.append(nitrogen)
    assert len(primary_amine_nitrogens)==num_of_amines
    starting_points = primary_amine_nitrogens
    to_skip = []
    for i in range(skip):
        new_starting_points = []
        for j in starting_points:
            atom = mol.GetAtomWithIdx(j)
            for another_atom in atom.GetNeighbors():
                another_idx = another_atom.GetIdx()
                if another_idx in carbons:
                    new_starting_points.append(another_idx)
                    to_skip.append(another_idx)
        starting_points = new_starting_points.copy()
    return list(set(carbons)-set(to_skip))

def find_primary_amine_pos(mol):
    nitrogens = find_atom_indices(mol, 7)
    return [i for i in nitrogens if num_of_Hs(mol.GetAtomWithIdx(i))==2]

def find_secondary_amine_pos(mol):
    nitrogens = find_atom_indices(mol, 7)
    return [i for i in nitrogens if num_of_Hs(mol.GetAtomWithIdx(i))==1]

def find_imine_pos(mol):

    '''
    Please make sure that mol has explicit Hs!!!!!
    '''
    
    nitrogens = find_atom_indices(mol, 7)
    _list = []
    for i in nitrogens:
        atom = mol.GetAtomWithIdx(i)
        if atom.GetDegree()==2 and num_of_Hs(atom)==1:
            _list.append(i)
    return _list

def remove_atom(mol, idx):
    new_mol = Chem.EditableMol(mol)
    new_mol.RemoveAtom(idx)
    return new_mol.GetMol()

def remove_unconnected_Hs(mol):
    Hs = [H for H in find_atom_indices(mol, 1) if len(mol.GetAtomWithIdx(H).GetNeighbors())==0]
    while len(Hs)>0:
        mol = remove_atom(mol, Hs[0])
        Hs = [H for H in find_atom_indices(mol, 1) if len(mol.GetAtomWithIdx(H).GetNeighbors())==0]
    return mol

def remove_NH2(mol):
    nitrogen = find_primary_amine_pos(mol)[0]
    mol_naked = remove_atom(mol, nitrogen)
    mol_naked = remove_unconnected_Hs(mol_naked)
    return mol_naked

def remove_one_H_from_NH2(mol):
    nitrogen = find_primary_amine_pos(mol)[0]
    H = [H.GetIdx() for H in mol.GetAtomWithIdx(nitrogen).GetNeighbors() if H.GetAtomicNum()==1][0]
    mol_naked = remove_atom(mol, H)
    return mol_naked

def remove_one_H_from_NH(mol):

    '''
    only deals with imines with no subs
    '''

    nitrogen = find_imine_pos(mol)[0]
    H = [H.GetIdx() for H in mol.GetAtomWithIdx(nitrogen).GetNeighbors() if H.GetAtomicNum()==1][0]
    mol_naked = remove_atom(mol, H)
    return mol_naked

def get_num_of_bonds(atom):
    bonds = atom.GetBonds()
    return sum([bond.GetBondTypeAsDouble() for bond in bonds])

def find_naked_atom_idx(mol):

    '''
    Only deals with naive cases. 
    H: 1 bond
    C: 4 bonds
    N: 3 bonds
    O: 2 bonds
    S: 2 or 6 bonds
    Halides: 1 bond
    P: 5 bonds
    '''

    num_of_bonds_dict = {1:[1], 6:[4], 7:[3], 8:[2], 16:[2,6], 
                        9:[1], 17:[1], 35:[1], 53:[1], 15:[5]}
    list_of_idx = []
    for atom in mol.GetAtoms():
        naked = False
        num = atom.GetAtomicNum()
        num_of_bonds = get_num_of_bonds(atom)
        charge = atom.GetFormalCharge()
        val = num_of_bonds - charge
        if val not in num_of_bonds_dict[num]:
            naked = True
        if num in [8, 16] and val==3: # this deals with furans and thiophenes
            naked = False
        if naked:
            list_of_idx.append(atom.GetIdx())
    return list_of_idx

# diamine -> diimine

def find_carbon_nitrogen_bond(nitrogen_atom):
    bonds = nitrogen_atom.GetBonds()
    atomic_num_pairs = [(bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()) for bond in bonds]
    for i, pair in enumerate(atomic_num_pairs):
        if 6 in pair:
            break
    bond = bonds[i]
    if pair[0] == 6:
        carbon = bond.GetBeginAtom()
    else:
        carbon = bond.GetEndAtom()
    return carbon, bond

def naive_diimine_fix(mol):

    '''
    re-calculates the bond orders simply by turning N-C=C-N into N=C-C=N
    '''

    mol = copy.deepcopy(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    nitrogens = find_naked_atom_idx(mol)
    assert len(nitrogens)==2
    nitrogen_atoms = [mol.GetAtomWithIdx(i) for i in nitrogens]
    carbon1, bond1 = find_carbon_nitrogen_bond(nitrogen_atoms[0])
    carbon2, bond2 = find_carbon_nitrogen_bond(nitrogen_atoms[1])
    CC_bond = mol.GetBondBetweenAtoms(carbon1.GetIdx(), carbon2.GetIdx())
    assert CC_bond.GetBondTypeAsDouble()==2
    CC_bond.SetBondType(Chem.rdchem.BondType.SINGLE)
    bond1.SetBondType(Chem.rdchem.BondType.DOUBLE)
    bond2.SetBondType(Chem.rdchem.BondType.DOUBLE)
    return mol

def rm_radicals(smi):

    '''
    Attempts to remove double radicals by generating resonance Lewis structures via kekulization. 
    Building a resonance supplier seems to be the only graceful way.
    '''

    # initial Lewis struture
    mol = Chem.MolFromSmiles(smi)
    Chem.Kekulize(mol, clearAromaticFlags=True) # aromatic bonds not allowed

    # resonance structure supplier
    supplier = Chem.ResonanceMolSupplier(mol, Chem.KEKULE_ALL)
    for res_mol in supplier:
        res_smi = Chem.MolToSmiles(res_mol, kekuleSmiles=True)
        new_mol = Chem.MolFromSmiles(res_smi)
        if Descriptors.NumRadicalElectrons(new_mol)==0: # breaks the moment a feasible SMILES string is found
            return res_smi
        
    try:
        return Chem.MolToSmiles(fix_NO2(mol))
    except:
        pass
    print('Ouch! What a tough one.\n' + smi)
    return smi   

def fix_NO2(mol):
    '''
    deals with the special case where multiple radicals exist, half of them on the N in NO2
    the other half on C atoms
    '''
    mol = copy.deepcopy(mol)
    radicals = [atom for atom in mol.GetAtoms() if atom.GetNumRadicalElectrons() > 0]
    atomic_num = sorted(set([atom.GetAtomicNum() for atom in radicals]))
    k = len(radicals)
    assert atomic_num==[6,7] and k%2==0
    radicals = sorted(radicals, key=lambda atom: atom.GetAtomicNum())
    cut = int(k/2)
    C_atoms, N_atoms = radicals[:cut], radicals[cut:]
    for C_atom, N_atom in zip(C_atoms, N_atoms):
        O_atom = None
        for atom in N_atom.GetNeighbors():
            if atom.GetAtomicNum()==8:
                O_atom = atom
                break
        bond = mol.GetBondBetweenAtoms(N_atom.GetIdx(), O_atom.GetIdx())
        bond.SetBondType(Chem.rdchem.BondType.DOUBLE)
        N_atom.SetNumRadicalElectrons(0)
        N_atom.SetFormalCharge(1)
        O_atom.SetFormalCharge(0)
        C_atom.SetNumRadicalElectrons(0)
        C_atom.SetFormalCharge(0)
    return mol

def primary_diamine_to_diimine(smi: str) -> str:

    '''
    smi: SMILES of primary diamine
    di_smi: SMILES of converted diimine

    '''

    mol = Chem.MolFromSmiles(smi)
    mol = rdmolops.AddHs(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)

    # generates conformer, this is absolutely needed to correctly determine connectivity
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)

    # removes two hydrogens
    mol = remove_one_H_from_NH2(mol)
    mol = remove_one_H_from_NH2(mol)

    try:
        mol = naive_diimine_fix(mol)
    except:
        rdDetermineBonds.DetermineBonds(mol, charge=0)

    # generates SMILES for diimine
    di_smi = Chem.MolToSmiles(mol)
    di_smi = Chem.CanonSmiles(di_smi)
    if Descriptors.NumRadicalElectrons(mol)>0:
        return Chem.CanonSmiles(rm_radicals(di_smi))
    return di_smi

# code for assembling

def connect_diamine_with_amine(mol_diamine, mol_amine):
    '''
    mol_diamine: a primary diamine
    mol_amine: a primary amine
    '''
    embed(mol_diamine)
    mol_diamine_naked = remove_one_H_from_NH2(mol=mol_diamine)
    mol_amine_naked = remove_NH2(mol=mol_amine)
    idx1, idx2 = find_naked_atom_idx(mol_diamine_naked)[0], find_naked_atom_idx(mol_amine_naked)[0]
    combo = Chem.EditableMol(Chem.CombineMols(mol_diamine_naked, mol_amine_naked))
    combo.AddBond(idx1, idx2 + mol_diamine_naked.GetNumAtoms(), order=Chem.rdchem.BondType.SINGLE)
    final_mol = combo.GetMol()
    return final_mol

def connect_diimine_with_amine(mol_diimine, mol_amine):
    '''
    mol_diimine: diimine with no subs
    mol_amine: a primary amine
    '''
    embed(mol_diimine)
    mol_diimine_naked = remove_one_H_from_NH(mol=mol_diimine)
    mol_amine_naked = remove_NH2(mol=mol_amine)
    idx1, idx2 = find_naked_atom_idx(mol_diimine_naked)[0], find_naked_atom_idx(mol_amine_naked)[0]
    combo = Chem.EditableMol(Chem.CombineMols(mol_diimine_naked, mol_amine_naked))
    combo.AddBond(idx1, idx2 + mol_diimine_naked.GetNumAtoms(), order=Chem.rdchem.BondType.SINGLE)
    final_mol = combo.GetMol()
    return final_mol

test_mol = Chem.MolFromSmiles('[2H]F')

def connect_amine_and_deuterium(mol1, idx1, mol2):

    '''
    please make sure that mol1 and mol2 already have explicit hydrogens!!!!!!!!!
    mol1: a Hydrogen atom will be removed from the atom indexed by idx1
    idx1: index of the carbon to be connected to mol2
    mol2: a D atom will be removed
    '''

    # removes one Hydrogen from mol1
    mol1_explicit_H = mol1
    for atom in mol1_explicit_H.GetAtomWithIdx(idx1).GetNeighbors():
        if (atom.GetAtomicNum() == 1):
            mol1_explicit_H = remove_atom(mol1_explicit_H, atom.GetIdx())
            break

    # removes D from mol2
    mol2_explicit_H = mol2
    for atom in mol2_explicit_H.GetAtoms():
        if atom.GetIsotope()>0:
            mol2_explicit_H = remove_atom(mol2_explicit_H, atom.GetIdx())
            break
    
    idx1 = find_naked_atom_idx(mol1_explicit_H)[0]
    idx2 = find_naked_atom_idx(mol2_explicit_H)[0]
    combo = Chem.EditableMol(Chem.CombineMols(mol1_explicit_H, mol2_explicit_H))
    combo.AddBond(idx1, idx2 + mol1_explicit_H.GetNumAtoms(), order=Chem.rdchem.BondType.SINGLE)
    final_mol = combo.GetMol()
    
    return clean(final_mol)
            
def find_unique_carbons(mol, num_of_amines):
    '''
    finds unique carbons for substitution
    '''
    carbons = find_carbons(mol=mol, num_of_amines=num_of_amines)
    list_of_mols = [connect_amine_and_deuterium(mol, i, test_mol) for i in carbons]
    idx = find_unique_mols(list_of_mols)
    return [carbons[i] for i in idx]

def find_connected_carbons(atom):

    '''
    finds idx of all carbons connected to atom
    '''

    atoms = atom.GetNeighbors()
    carbons = [i.GetIdx() for i in atoms if i.GetAtomicNum()==6]
    return sorted(carbons)

def bfs_find_carbons(mol, starting_idx, depth):

    '''
    naive bfs algorithm that searches for carbons that are i bonds away from the input starting_idx
    '''

    explored = set()
    list_of_carbons_per_level = []
    for i in range(depth):
        this_level = []
        for i in starting_idx:
            atom = mol.GetAtomWithIdx(i)
            this_level += find_connected_carbons(atom)
        this_level = sorted(set(this_level)-explored)
        list_of_carbons_per_level.append(this_level)
        explored.update(this_level)
        starting_idx = this_level
    return list_of_carbons_per_level

def filter_carbons(mol, list_of_carbons):
    '''
    finds carbon candidates for substitution
    list_of_carbons: atom indices
    '''
    carbons = [carbon for carbon in list_of_carbons if num_of_Hs(mol.GetAtomWithIdx(carbon))>0]
    return carbons