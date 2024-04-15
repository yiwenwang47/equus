from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
import networkx as nx
from networkx.algorithms import isomorphism

def pure_mol_to_nx(mol):
    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                   atom=atom.GetAtomicNum())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx())
    return G

def mol_to_nx(mol):
    if sum(atom.GetImplicitValence() for atom in mol.GetAtoms()) > 0:
        mol_H = rdmolops.AddHs(mol)
    else:
        mol_H = mol
    G = pure_mol_to_nx(mol_H)
    return G

node_match = lambda node1, node2: node1['atom']==node2['atom']

def find_unique_mols(list_of_mols):
    idx = []
    list_of_graphs = [pure_mol_to_nx(mol) for mol in list_of_mols]
    unique_graphs = []
    for i, graph in enumerate(list_of_graphs):
        keep = True
        for other_graph in unique_graphs:
            if isomorphism.GraphMatcher(graph, other_graph, node_match=node_match).is_isomorphic():
                keep = False
                break
        if keep:
            unique_graphs.append(graph)
            idx.append(i)
    return idx

def find_unique_smiles(list_of_smiles):
    list_of_mols = [Chem.MolFromSmiles(smiles) for smiles in list_of_smiles]
    idx = find_unique_mols(list_of_mols)
    return idx