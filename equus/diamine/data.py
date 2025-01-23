from __future__ import annotations

from collections import defaultdict
from typing import Optional

import pandas as pd
from pandas.core.frame import DataFrame
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, rdmolops
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

Chem.SetUseLegacyStereoPerception(False)

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


_form = lambda mol: CalcMolFormula(mol)
_n_bits = 2048
_radius = 5
_fpgen_stereo = AllChem.GetMorganGenerator(
    radius=_radius, fpSize=_n_bits, includeChirality=True
)
_fpgen = AllChem.GetMorganGenerator(
    radius=_radius, fpSize=_n_bits, includeChirality=False
)


class Molecules(object):

    """
    Molecule pool.
    Has a search method.
    """

    def __init__(self, names: Optional[list] = None, smiles: Optional[list] = None):
        self.names = names if names else []
        self.smiles = [Chem.CanonSmiles(smi) for smi in smiles] if smiles else []
        self.mols = [Chem.MolFromSmiles(smi) for smi in self.smiles]

        self.mol_form_dict = defaultdict(list)
        for i, mol in enumerate(self.mols):
            formula = _form(mol)
            self.mol_form_dict[formula].append(i)

        self.fp_generator = _fpgen
        self.fp_generator_stereo = _fpgen_stereo

        self.fps = [self.fp_generator.GetFingerprint(mol) for mol in self.mols]
        self.fps_stereo = [
            self.fp_generator_stereo.GetFingerprint(mol) for mol in self.mols
        ]
        self.n = len(self.smiles)
        self.graphs = [pure_mol_to_nx(mol) for mol in self.mols]

    def __getitem__(self, i: int | str) -> tuple[str, str]:
        if type(i) == str:
            i = self.names.index(i)
        return (self.names[i], self.smiles[i])

    def __len__(self) -> int:
        return self.n

    def __repr__(self) -> str:
        return f"A pool containing {self.n} molecules."

    def head(self, k: int = 10):
        return Draw.MolsToGridImage(
            mols=[Chem.MolFromSmiles(smi) for smi in self.smiles[:k]]
        )

    def add(self, name: str, smi: str):
        smi = Chem.CanonSmiles(smi)
        self.names.append(name)
        self.smiles.append(smi)
        mol = Chem.MolFromSmiles(smi)
        self.mols.append(mol)
        self.fps.append(self.fp_generator.GetFingerprint(mol))
        self.fps_stereo.append(self.fp_generator_stereo.GetFingerprint(mol))
        formula, i = _form(mol), self.n
        self.mol_form_dict[formula].append(i)
        self.graphs.append(pure_mol_to_nx(mol))
        self.n += 1

    def search(
        self, smi: str, verbose=False, use_stereo: bool = True
    ) -> Optional[list]:
        r"""
        Returns:
        True, name, smi if molecule is found.
        False, "", "" if not found.
        """

        # naive smiles matching
        smi = Chem.CanonSmiles(smi)
        if use_stereo:
            if smi in self.smiles:
                i = self.smiles.index(smi)
                return [(self.names[i], self.smiles[i])]

        if verbose:
            print(
                "Did not find exact SMILES match. Moving on to molecular formula matching."
            )

        # find candidates by molecular formula matching
        mol = Chem.MolFromSmiles(smi)
        formula = _form(mol)
        if formula not in self.mol_form_dict:
            return None

        if verbose:
            print(
                "Found candidates by molecular formula matching. Moving on to fingerprint matching."
            )

        # filter the list of candidates by Tanimoto Similarity of Morgan fingerprints
        candidates = self.mol_form_dict[formula]

        if use_stereo:
            this_fp = self.fp_generator_stereo.GetFingerprint(mol)
            candidates = [
                i
                for i in candidates
                if DataStructs.TanimotoSimilarity(self.fps_stereo[i], this_fp) == 1.0
            ]
        else:
            this_fp = self.fp_generator.GetFingerprint(mol)
            candidates = [
                i
                for i in candidates
                if DataStructs.TanimotoSimilarity(self.fps[i], this_fp) == 1.0
            ]

        if len(candidates) == 0:
            return None

        if verbose:
            print(
                "Found candidates by fingerprint matching. Moving on to isomorphism matching."
            )
            print("Number of candidates:", len(candidates))

        # final results filtered by isomorphism
        results = []
        this_graph = pure_mol_to_nx(mol)
        for i in candidates:
            ref_mol = Chem.MolFromSmiles(self.smiles[i])
            if use_stereo:
                ref_graph = self.graphs[i]
                result = nx_isomorphism(
                    this_graph, ref_graph
                )  # this seems to be much more accurate than mol1.HasSubstructMatch(mol2, useChirality=True)
            else:
                result = mol.HasSubstructMatch(ref_mol) and ref_mol.HasSubstructMatch(
                    mol
                )
            if result:
                return results.append((self.names[i], self.smiles[i]))
        if len(results) > 0:
            return results
        else:
            return None

    @staticmethod
    def from_csv(
        filename: str,
        name_col_id: int = 0,
        smiles_col_id: int = 1,
    ) -> Molecules:
        df = pd.read_csv(filename)
        return Molecules.from_df(
            df=df,
            name_col_id=name_col_id,
            smiles_col_id=smiles_col_id,
        )

    @staticmethod
    def from_df(
        df: DataFrame,
        name_col_id: int = 0,
        smiles_col_id: int = 1,
    ) -> Molecules:
        names = list(df.values[:, name_col_id])
        smiles = list(df.values[:, smiles_col_id])
        return Molecules(names, smiles)

    def to_df(self, name_col: str = "name", smiles_col: str = "SMILES") -> DataFrame:
        df = pd.DataFrame(data={name_col: self.names, smiles_col: self.smiles})
        return df
