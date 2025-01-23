from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from multiprocessing import Process, Queue

import pandas as pd
from pandas.core.frame import DataFrame
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, rdmolops
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from equus.diamine.iso import nx_isomorphism, pure_mol_to_nx

# def _has_substruct_match(q: Queue, mol: Chem.Mol, pattern: Chem.Mol, useChirality=True):
#     q.put(mol.HasSubstructMatch(pattern, useChirality=useChirality))


# def _time_out_mol_isomorphism_stereo(
#     mol: Chem.Mol, pattern: Chem.Mol, patience: float = 0.5
# ):
#     r"""
#     A time-out mechanism for the HasSubstructMatch check.

#     Returns None if the process is terminated after reaching the time limit (patience).

#     """
#     q = Queue()
#     p = Process(target=_has_substruct_match, args=(q, mol, pattern, True))
#     p.start()
#     p.join(timeout=patience)
#     res = None
#     if p.exitcode is None:
#         p.terminate()
#     else:
#         res = q.get()
#     return res


# def canonicalize_smiles(smi: str, verbose: bool = False, max_attempt: int = 3) -> str:
#     r"""
#     A naive self-consistent method that does canonicalization of SMILES strings.

#     Args:
#         smi: str, the SMILES string to be canonicalized.
#         verbose: bool, whether to print out the steps.
#         max_attempt: int, the maximum number of attempts to reach a canonical SMILES string.

#     Returns:
#         str, the canonicalized SMILES string.
#     """
#     counter = 0
#     while True:
#         new_smi = Chem.CanonSmiles(smi)
#         if counter >= max_attempt:
#             if verbose:
#                 print(f"Reached the maximum number of attempts ({max_attempt}).")
#             smi = new_smi
#             break
#         if new_smi == smi:
#             break
#         else:
#             if verbose:
#                 print(f"Oops, {smi} is not canonical. Let's try {new_smi}.")
#             smi = new_smi
#             counter += 1
#     return smi


_form = lambda mol: CalcMolFormula(mol)
_n_bits = 2048
_radius = 5
_fpgen_stereo = AllChem.GetMorganGenerator(
    radius=_radius, fpSize=_n_bits, includeChirality=True
)
_fpgen = AllChem.GetMorganGenerator(
    radius=_radius, fpSize=_n_bits, includeChirality=False
)


@dataclass
class Molecules:

    """
    Molecule pool.
    Has a search method.
    """

    names: list[str]
    smiles: list[str]
    use_stereo: bool = True

    def __post_init__(self):
        self.smiles = [Chem.CanonSmiles(smi) for smi in self.smiles]
        self.mols = [Chem.MolFromSmiles(smi) for smi in self.smiles]

        self.mol_form_dict = defaultdict(list)
        for i, mol in enumerate(self.mols):
            formula = _form(mol)
            self.mol_form_dict[formula].append(i)

        if self.use_stereo:
            self.fp_generator = _fpgen_stereo
        else:
            self.fp_generator = _fpgen

        self.fps = [self.fp_generator.GetFingerprint(mol) for mol in self.mols]
        self.n = len(self.smiles)
        self.graphs = [pure_mol_to_nx(mol) for mol in self.mols]

    def __getitem__(self, i: int | str) -> tuple[str, str]:
        if type(i) == str:
            i = self.names.index(i)
        return (self.names[i], self.smiles[i])

    def __len__(self) -> int:
        return self.n

    def head(self, k: int = 10):
        return Draw.MolsToGridImage(
            mols=[Chem.MolFromSmiles(smi) for smi in self.smiles[:k]]
        )

    def add(self, name: str, smi: str):
        self.names.append(name)
        self.smiles.append(smi)
        mol = Chem.MolFromSmiles(smi)
        self.mols.append(mol)
        self.fps.append(self.fp_generator.GetFingerprint(mol))
        formula, i = _form(mol), self.n
        self.mol_form_dict[formula].append(i)
        self.graphs.append(pure_mol_to_nx(mol))
        self.n += 1

    def search(self, smi: str, verbose=False) -> tuple[bool, str, str]:
        r"""
        Returns:
        True, name, smi if molecule is found.
        False, "", "" if not found.
        """

        # naive smiles matching
        smi = Chem.CanonSmiles(smi)
        if smi in self.smiles:
            i = self.smiles.index(smi)
            return True, self.names[i], self.smiles[i]

        if verbose:
            print(
                "Did not find exact SMILES match. Moving on to molecular formula matching."
            )

        # find candidates by molecular formula matching
        mol = Chem.MolFromSmiles(smi)
        formula = _form(mol)
        if formula not in self.mol_form_dict:
            return False, "", ""

        if verbose:
            print(
                "Found candidates by molecular formula matching. Moving on to fingerprint matching."
            )

        # filter the list of candidates by Tanimoto Similarity of Morgan fingerprints
        candidates = self.mol_form_dict[formula]
        this_fp = self.fp_generator.GetFingerprint(mol)
        candidates = [
            i
            for i in candidates
            if DataStructs.TanimotoSimilarity(self.fps[i], this_fp) == 1.0
        ]
        if len(candidates) == 0:
            return False, "", ""

        if verbose:
            print(
                "Found candidates by fingerprint matching. Moving on to isomorphism matching."
            )
            print("Number of candidates:", len(candidates))

        this_graph = pure_mol_to_nx(mol)
        for i in candidates:
            ref_mol = Chem.MolFromSmiles(self.smiles[i])
            if self.use_stereo:
                ref_graph = self.graphs[i]
                result = nx_isomorphism(
                    this_graph, ref_graph
                )  # this seems to be much more accurate than mol1.HasSubstructMatch(mol2, useChirality=True)
            else:
                result = mol.HasSubstructMatch(ref_mol) and ref_mol.HasSubstructMatch(
                    mol
                )
            if result:
                return True, self.names[i], self.smiles[i]

        return False, "", ""

    @staticmethod
    def from_csv(
        filename: str,
        name_col_id: int = 0,
        smiles_col_id: int = 1,
        use_stereo: bool = True,
    ) -> Molecules:
        df = pd.read_csv(filename)
        return Molecules.from_df(
            df=df,
            name_col_id=name_col_id,
            smiles_col_id=smiles_col_id,
            use_stereo=use_stereo,
        )

    @staticmethod
    def from_df(
        df: DataFrame,
        name_col_id: int = 0,
        smiles_col_id: int = 1,
        use_stereo: bool = True,
    ) -> Molecules:
        names = list(df.values[:, name_col_id])
        smiles = list(df.values[:, smiles_col_id])
        return Molecules(names, smiles, use_stereo=use_stereo)

    def to_df(self, name_col: str = "name", smiles_col: str = "SMILES") -> DataFrame:
        df = pd.DataFrame(data={name_col: self.names, smiles_col: self.smiles})
        return df
