from collections.abc import Iterator

import networkx as nx
from networkx.algorithms import isomorphism
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.rdchem import Atom, Bond, Mol

from equus.edit.iso import node_match
from equus.edit.utils import (
    clean,
    find_atom_indices,
    find_naked_atom_idx,
    find_primary_amine_pos,
    get_num_of_bonds,
    num_of_Hs,
    pure_mol_to_nx,
    read_smiles,
    remove_atom,
    remove_unconnected_Hs,
)

"""
Note: all molecules should have explicit Hydrogens!
"""

mol_H = Chem.MolFromSmiles("[H]", sanitize=False)


def add_one_H(mol_naked: Mol) -> Mol:
    """
    mol_naked: molecule with at least one naked atom
    """
    i = next(find_naked_atom_idx(mol_naked))
    combo = Chem.EditableMol(Chem.CombineMols(mol_naked, mol_H))
    combo.AddBond(i, mol_naked.GetNumAtoms(), order=Chem.rdchem.BondType.SINGLE)
    mol = clean(combo.GetMol())
    return mol


def find_OH(mol: Mol) -> Iterator[int]:
    oxygens = find_atom_indices(mol, atomic_number=8)
    oxygens = (i for i in oxygens if num_of_Hs(mol.GetAtomWithIdx(i)) == 1)
    return oxygens


def remove_OH(mol: Mol) -> Mol:
    """
    mol: molecule with explicit Hs
    replaces all -OH groups with Hydrogens
    """
    while True:
        try:
            O = next(find_OH(mol))
            mol = remove_atom(mol, O)
            mol = remove_unconnected_Hs(mol)
            mol = add_one_H(mol)
        except:
            break
    return mol


def find_SH(mol) -> Iterator[int]:
    sulfurs = find_atom_indices(mol, atomic_number=16)
    sulfurs = (
        i
        for i in sulfurs
        if num_of_Hs(mol.GetAtomWithIdx(i)) == 1
        and get_num_of_bonds(mol.GetAtomWithIdx(i)) == 2
    )
    return sulfurs


def remove_SH(mol: Mol) -> Mol:
    """
    mol: molecule with explicit Hs
    replaces all -SH groups with Hydrogens
    """
    while True:
        try:
            S = next(find_SH(mol))
            mol = remove_atom(mol, idx=S)
            mol = remove_unconnected_Hs(mol)
            mol = add_one_H(mol)
        except:
            break
    return mol


def remove_unwanted_NH2(mol: Mol, N_idx: int) -> Mol:
    """
    mol: molecule with explicit Hs
    N_idx: idx of N atom in the -NH2 to be removed
    """
    mol = remove_atom(mol, N_idx)
    mol = remove_unconnected_Hs(mol)
    mol = add_one_H(mol)

    return mol


def remove_all_NH2(mol: Mol) -> Mol:
    """
    mol: molecule with explicit Hs
    all -NH2 groups will be removed
    """

    while True:
        try:
            N = next(find_primary_amine_pos(mol))
            mol = remove_unwanted_NH2(mol, N)
        except:
            break

    return mol


def remove_all(mol: Mol) -> Mol:
    mol = remove_all_NH2(mol)
    mol = remove_OH(mol)
    mol = remove_SH(mol)
    return mol


# graph of N-C-C-N
bridge_atoms = [7, 6, 6, 7]
bridge_graph = nx.Graph()
for i, atom in enumerate(bridge_atoms):
    bridge_graph.add_node(i, atom=atom)
    if i > 0:
        bridge_graph.add_edge(i - 1, i)

helper = lambda mapping: [i[0] for i in sorted(mapping.items(), key=lambda x: x[1])]


def verify_N_atom(atom: Atom) -> bool:

    """
    Criterion: at least one single bond. For this purpose, aromatic flags are needed.
    """

    return any(bond.GetBondTypeAsDouble() == 1 for bond in atom.GetBonds())


def find_bridges(mol: Mol) -> list:

    """
    Finds all N-C-C-N bridges.
    Hydrogens are not considered.
    """

    # mappings
    graph = pure_mol_to_nx(mol)
    matcher = isomorphism.GraphMatcher(graph, bridge_graph, node_match)
    assert matcher.subgraph_is_isomorphic()
    mappings = [i for i in matcher.subgraph_isomorphisms_iter()]

    # sorts each mapping in the order of N-C-C-N
    bridges_found = [helper(mapping) for mapping in mappings]

    # finds unique mappings
    unique_bridges = []
    sorted_bridges = []
    for i in bridges_found:
        N1, N2 = mol.GetAtomWithIdx(i[0]), mol.GetAtomWithIdx(i[-1])
        if verify_N_atom(N1) and verify_N_atom(N2):
            sorted_bridge = sorted(i)
            if sorted_bridge not in sorted_bridges:
                unique_bridges.append(i)
                sorted_bridges.append(sorted_bridge)
    return unique_bridges


def normalize_primary_diamine(mol: Mol) -> Mol:

    """
    Keeps the N-C-C-N bridge substructure.
    Removes all other -NH2 groups.
    """

    smi = Chem.MolToSmiles(mol)
    mol = read_smiles(smi, no_aromatic_flags=False, hydrogens=True)

    def helper(mol):
        bridges = find_bridges(mol)
        nitrogens = [N for N in find_primary_amine_pos(mol)]
        true_bridges = []
        for bridge in bridges:
            if bridge[0] in nitrogens and bridge[-1] in nitrogens:
                bond = mol.GetBondBetweenAtoms(bridge[1], bridge[2])
                N1, N2 = mol.GetAtomWithIdx(bridge[0]), mol.GetAtomWithIdx(bridge[-1])
                if (
                    bond.GetBondTypeAsDouble() > 1
                    and N1.GetFormalCharge() == 0
                    and N2.GetFormalCharge() == 0
                ):
                    true_bridges.append(bridge)
        assert len(true_bridges) == 1
        bridge = true_bridges[0]
        unwanted = list(set(nitrogens) - set(bridge))
        return unwanted

    unwanted = helper(mol)

    while len(unwanted) > 0:
        mol = remove_unwanted_NH2(mol, unwanted[0])
        unwanted = helper(mol)

    smi = Chem.CanonSmiles(Chem.MolToSmiles(mol=mol))
    return read_smiles(smi, no_aromatic_flags=False, hydrogens=True)


acyl_atoms = [1, 1, 7, 6, 8]
acyl_graph = nx.Graph()
for i, atom in enumerate(acyl_atoms):
    acyl_graph.add_node(i, atom=atom)
acyl_graph.add_edge(0, 2)
for i in range(1, 4):
    acyl_graph.add_edge(i, i + 1)


def is_acyl(mol: Mol) -> bool:

    """
    Checks whether mol has an acyl group right next to -NH2.
    mol must have explicit hydrogens.
    """

    graph = pure_mol_to_nx(mol)
    matcher = isomorphism.GraphMatcher(graph, acyl_graph, node_match)
    if not matcher.subgraph_is_isomorphic():
        return False
    mappings = [i for i in matcher.subgraph_isomorphisms_iter()]
    for mapping in mappings:
        C, O = helper(mapping)[-2:]
        bond = mol.GetBondBetweenAtoms(C, O)
        if bond.GetBondTypeAsDouble() == 2:
            return True
    return False


urea_atoms = [1, 7, 6, 7, 1, 8]
urea_graph = nx.Graph()
for i, atom in enumerate(urea_atoms):
    urea_graph.add_node(i, atom=atom)
    if i > 0 and i != 5:
        urea_graph.add_edge(i, i - 1)
    if i == 5:
        urea_graph.add_edge(2, i)


def is_urea(mol: Mol) -> bool:

    """
    Checks whether mol has a urea substructure.
    mol must have explicit hydrogens.
    """

    graph = pure_mol_to_nx(mol)
    matcher = isomorphism.GraphMatcher(graph, urea_graph, node_match)
    if not matcher.subgraph_is_isomorphic():
        return False
    mappings = [i for i in matcher.subgraph_isomorphisms_iter()]
    for mapping in mappings:
        plum = helper(mapping)
        C, O = plum[2], plum[-1]
        bond = mol.GetBondBetweenAtoms(C, O)
        if bond.GetBondTypeAsDouble() == 2:
            return True
    return False
