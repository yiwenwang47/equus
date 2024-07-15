import random

import numpy as np
from rdkit import Chem
from rdkit.Chem.rdchem import Mol

import equus.backbone as bb
import equus.diamine.normalize as norm
from equus.diamine import clean, read_smiles, to_smiles
from equus.diamine.data import Molecules


def lazy_editor(prob=0.5):
    def decorator(func):
        def helper(mol):
            try:
                if random.random() <= prob:
                    new_mol = func(mol)
                    assert new_mol is not None
                    return new_mol
                else:
                    return clean(mol)
            except:
                print("Didn't happen.")
                return clean(mol)

        return helper

    return decorator


@lazy_editor(prob=0.5)
def reduce(mol: Mol) -> Mol:
    return bb.randomly_reduce_bonds(mol)


@lazy_editor(prob=0.5)
def break_1_bond(mol: Mol) -> Mol:
    return bb.break_random_single_bond(mol)


@lazy_editor(prob=0.5)
def add_carbons(mol: Mol) -> Mol:

    carbon_idx = bb.find_carbons_by_num_of_Hs(mol, min_num_of_Hs=1, max_num_of_Hs=2)
    assert len(carbon_idx) > 0
    idx = random.choice(carbon_idx)

    funcs = [bb.add_methyl, bb.add_ethyl, bb.add_isopropyl]
    # e^(- num of carbons)
    weights = np.array([1, 2, 4])
    weights = np.exp(-weights)
    weights /= weights.sum()
    add_func = funcs[np.random.choice([0, 1, 2], 1, p=weights)[0]]

    return add_func(mol, idx)


@lazy_editor(prob=0.5)
def oxidize(mol: Mol, verbose: bool = False) -> Mol:

    methyls = bb.find_carbons_by_num_of_Hs(mol, min_num_of_Hs=3, max_num_of_Hs=3)
    CH2s = bb.find_carbons_by_num_of_Hs(mol, min_num_of_Hs=2, max_num_of_Hs=2)
    single_bonds = bb.find_carbon_carbon_single_bonds(
        mol, min_num_of_Hs=1, filter_for_oxidation=True
    )

    funcs = [
        bb.CH3_to_CCH,
        bb.CH3_to_CN,
        bb.CH2_to_CO,
        bb.CH2_to_CCH2,
        bb.oxidize_single_bond,
    ]

    masks = np.array([1, 1, 1, 1, 1])
    try:
        assert len(methyls) > 0 or len(CH2s) > 0 or len(single_bonds) > 0
        if len(methyls) == 0:
            masks[0] = 0
            masks[1] = 0
        if len(CH2s) == 0:
            masks[2] = 0
            masks[3] = 0
        if len(single_bonds) == 0:
            masks[4] = 0

        weights = masks / masks.sum()
        k = np.random.choice(np.arange(5), 1, p=weights)[0]
        func = funcs[k]
        if k in [0, 1]:
            idx = random.choice(methyls)
            return func(mol, idx)
        if k in [2, 3]:
            idx = random.choice(CH2s)
            return func(mol, idx)
        idx = random.choice(single_bonds)
        return func(mol, idx)
    except:
        if verbose:
            print("Nothing to oxidize.")
        return mol


@lazy_editor(prob=0.25)
def add_halogen(mol: Mol) -> Mol:
    carbon_idx = bb.find_carbons_by_num_of_Hs(mol, min_num_of_Hs=1, max_num_of_Hs=3)
    idx = random.choice(carbon_idx)
    halogen = random.choice(["F", "Cl", "Br", "I"])
    return bb.add_halogen(mol=mol, idx=idx, halogen=halogen)


@lazy_editor(prob=0.5)
def add_misc(mol: Mol) -> Mol:
    carbon_idx = bb.find_carbons_by_num_of_Hs(mol, min_num_of_Hs=1, max_num_of_Hs=3)
    idx = random.choice(carbon_idx)
    return bb.add_group(mol=mol, idx=idx)


@lazy_editor(prob=0.5)
def replace_carbon(mol: Mol) -> Mol:
    n_N = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    n_O = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_S = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    weights = np.array([n_N, n_O, n_S])
    weights = np.exp(-weights)
    weights /= weights.sum()
    atomic_num = int(np.random.choice([7, 8, 16], 1)[0])
    return bb.replace_carbon_atom(mol, atomic_num=atomic_num)


def pseudo_diffusion(mol: Mol, t: int = 5) -> Mol:
    """
    t: number of random operations
    Not all of them will happen!
    """
    func_candidates = [
        reduce,
        break_1_bond,
        add_carbons,
        oxidize,
        add_halogen,
        add_misc,
        replace_carbon,
    ]
    for _ in range(t):
        func = random.choice(func_candidates)
        mol = func(mol)
    mol = norm.remove_all(mol)
    return mol


def verify_backbone(smi: str) -> bool:
    """
    verifies if the backbone can be used to generate a diamine
    """
    mol_aromatic = read_smiles(smi, no_aromatic_flags=False)
    mol = read_smiles(smi, no_aromatic_flags=True)

    pattern_aromatic = Chem.MolFromSmarts("[#6H1]:[#6H1]")
    pattern = Chem.MolFromSmarts("[#6H1]=[#6H1]")

    if any(mol.GetSubstructMatches(pattern)) and any(
        mol_aromatic.GetSubstructMatches(pattern_aromatic)
    ):
        return True
    else:
        return False


def generate_random_permutations(smi: str, n: int = 10) -> list[str]:
    """
    attempts to generate n random new structures
    however, most of them will be invalid
    """

    mol = read_smiles(smi, no_aromatic_flags=True, hydrogens=True)
    pool = Molecules([], [])
    pool.add("0", smi)

    for i in range(40):
        new_mol = pseudo_diffusion(mol, t=10)
        new_smi = to_smiles(new_mol)
        if bb.validate(new_smi) and verify_backbone(new_smi):
            _, found = pool.search(new_smi)
            if not found:
                pool.add(str(pool.n), new_smi)

    return pool.smiles[1:]
