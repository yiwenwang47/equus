import networkx as nx
from networkx.algorithms import isomorphism
from rdkit import Chem
from rdkit.Chem import rdmolops


def pure_mol_to_nx(mol):
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(
            atom.GetIdx(),
            atom=atom.GetAtomicNum(),
            atom_stereo=str(atom.GetChiralTag()),
            atom_hybrid=str(atom.GetHybridization()),
            atom_charge=atom.GetFormalCharge(),
            num_Hs=atom.GetTotalNumHs(),
            num_radical_electrons=atom.GetNumRadicalElectrons(),
        )
    for bond in mol.GetBonds():
        G.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            order=bond.GetBondTypeAsDouble(),
            bond_stereo=str(bond.GetStereo()),
        )
    return G


def node_match(node1, node2):
    if node1["atom"] != node2["atom"]:
        return False
    if node1["atom_stereo"] != node2["atom_stereo"]:
        return False
    if node1["atom_hybrid"] != node2["atom_hybrid"]:
        return False
    if node1["atom_charge"] != node2["atom_charge"]:
        return False
    if node1["num_Hs"] != node2["num_Hs"]:
        return False
    if node1["num_radical_electrons"] != node2["num_radical_electrons"]:
        return False
    return True


def edge_match(edge1, edge2):
    if edge1["order"] != edge2["order"]:
        return False
    if edge1["bond_stereo"] != edge2["bond_stereo"]:
        return False
    return True


def nx_isomorphism(G1: nx.Graph, G2: nx.Graph) -> bool:
    return isomorphism.GraphMatcher(
        G1, G2, node_match=node_match, edge_match=edge_match
    ).is_isomorphic()


def find_unique_mols(list_of_mols):
    idx = []
    list_of_graphs = [pure_mol_to_nx(mol) for mol in list_of_mols]
    unique_graphs = []
    for i, graph in enumerate(list_of_graphs):
        keep = True
        for other_graph in unique_graphs:
            if isomorphism.GraphMatcher(
                graph, other_graph, node_match=node_match
            ).is_isomorphic():
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
