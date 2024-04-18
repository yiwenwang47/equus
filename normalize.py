from rdkit import Chem
from rdkit.Chem import rdmolops
from .edit import clean, find_atom_indices, find_naked_atom_idx, find_primary_amine_pos
from .edit import num_of_Hs, remove_atom, remove_unconnected_Hs
from .edit import pure_mol_to_nx
import networkx as nx
from .iso import node_match
from networkx.algorithms import isomorphism

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

def find_OH(mol):
    oxygens = find_atom_indices(mol, atomic_number=8)
    oxygens = [i for i in oxygens if num_of_Hs(mol.GetAtomWithIdx(i))==1]
    return oxygens

def remove_OH(mol):

    '''
    mol: molecule with explicit Hs
    replaces all -OH groups with Hydrogens
    '''

    oxygens = find_OH(mol)
    while len(oxygens)>0:
        mol = remove_atom(mol, oxygens[0])
        mol = remove_unconnected_Hs(mol)
        mol = add_one_H(mol)
        oxygens = find_OH(mol)
    return mol

def find_SH(mol):
    sulfurs = find_atom_indices(mol, atomic_number=16)
    sulfurs = [i for i in sulfurs if num_of_Hs(mol.GetAtomWithIdx(i))==1]
    return sulfurs

def remove_SH(mol):

    '''
    mol: molecule with explicit Hs
    replaces all -SH groups with Hydrogens
    '''

    sulfurs = find_SH(mol)
    while len(sulfurs)>0:
        mol = remove_atom(mol, idx=sulfurs[0])
        mol = remove_unconnected_Hs(mol)
        mol = add_one_H(mol)
        sulfurs = find_SH(mol)
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

def remove_all_NH2(mol):

    '''
    mol: molecule with explicit Hs
    all -NH2 groups will be removed
    '''

    nitrogens = find_primary_amine_pos(mol)
    
    while len(nitrogens)>0:
        mol = remove_unwanted_NH2(mol, nitrogens[0])
        nitrogens = find_primary_amine_pos(mol)
    
    return mol

def remove_all(mol):
    mol = remove_all_NH2(mol)
    mol = remove_OH(mol)
    mol = remove_SH(mol)
    return mol

# graph of N-C-C-N
bridge_atoms = [7, 6, 6, 7]
bridge_graph = nx.Graph()
for i, atom in enumerate(bridge_atoms):
    bridge_graph.add_node(i, atom=atom)
    if i>0:
        bridge_graph.add_edge(i-1, i)

def find_bridges(mol):

    '''
    Finds all N-C-C-N bridges. Hydrogens are not considered.
    '''

    # mappings
    graph = pure_mol_to_nx(mol)
    matcher = isomorphism.GraphMatcher(graph, bridge_graph, node_match)
    assert matcher.subgraph_is_isomorphic()
    mappings = [i for i in matcher.subgraph_isomorphisms_iter()]

    # sorts each mapping in the order of N-C-C-N
    helper = lambda mapping: [i[0] for i in sorted(mapping.items(), key=lambda x: x[1])]
    bridges_found = [helper(mapping) for mapping in mappings]

    # finds unique mappings
    unique_bridges = []
    sorted_bridges = []
    for i in bridges_found:
        sorted_bridge = sorted(i)
        if sorted_bridge not in sorted_bridges:
            unique_bridges.append(i)
            sorted_bridges.append(sorted_bridge)

    return unique_bridges