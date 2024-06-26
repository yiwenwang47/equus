from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from networkx.algorithms import isomorphism
from pandas.core.frame import DataFrame
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from equus.edit.iso import node_match, pure_mol_to_nx
from equus.edit.utils import read_smiles

__form = lambda mol: Chem.rdMolDescriptors.CalcMolFormula(mol)
__n_bits = 2048
__radius = 5
__fpgen = AllChem.GetMorganGenerator(radius=__radius, fpSize=__n_bits)


@dataclass
class Molecules:

    """
    Class for keeping track of molecules.
    Has a search method.
    """

    names: list[str]
    smiles: list[str]

    def __post_init__(self):
        self.mols = [
            read_smiles(smiles_string=smi, no_aromatic_flags=False, hydrogens=True)
            for smi in self.smiles
        ]
        self.mol_form = [__form(mol) for mol in self.mols]
        self.fps = [__fpgen.GetFingerprint(mol) for mol in self.mols]
        self.n = len(self.smiles)

    def __getitem__(self, i: int | str) -> tuple[str, str]:
        if type(i) == str:
            i = self.names.index(i)
        return (self.names[i], self.smiles[i])

    def __len__(self) -> int:
        return self.n

    def add(self, name: str, smi: str):
        self.names.append(name)
        self.smiles.append(smi)
        mol = read_smiles(smi, no_aromatic_flags=False, hydrogens=True)
        self.mols.append(mol)
        self.mol_form.append(__form(mol))
        self.fps.append(__fpgen.GetFingerprint(mol))
        self.n += 1

    def search(self, smi: str) -> tuple[bool, str]:
        """
        Returns:
        True, name if molecule is found.
        False, '' if not found.
        """

        # naive smiles matching
        smi = Chem.CanonSmiles(smi)
        try:
            i = self.smiles.index(smi)
            return True, self.name[i]
        except:
            pass

        # find candidates by molecular formula matching
        mol = read_smiles(smi, no_aromatic_flags=False, hydrogens=True)
        formula = __form(mol)
        if formula not in self.mol_form:
            return False, ""

        # filter the list of candidates by Tanimoto Similarity of Morgan fingerprints
        candidates = [i for i, x in enumerate(self.mol_form) if x == formula]
        this_fp = __fpgen.GetFingerprint(mol)
        candidates = [
            i
            for i in candidates
            if DataStructs.TanimotoSimilarity(self.fps[i], this_fp) > 0.99
        ]
        if len(candidates) == 0:
            return False, ""

        this_graph = pure_mol_to_nx(mol)
        for i in candidates:
            ref_graph = pure_mol_to_nx(self.mols[i])
            if isomorphism.GraphMatcher(
                this_graph, ref_graph, node_match=node_match
            ).is_isomorphic():
                return True, self.names[i]

        return False, ""

    @staticmethod
    def from_csv(
        filename: str, name_col_id: int = 0, smiles_col_id: int = 1
    ) -> Molecules:
        df = pd.read_csv(filename)
        return Molecules.from_df(
            df=df, name_col_id=name_col_id, smiles_col_id=smiles_col_id
        )

    @staticmethod
    def from_df(
        df: DataFrame, name_col_id: int = 0, smiles_col_id: int = 1
    ) -> Molecules:
        names = list(df.values[:, name_col_id])
        smiles = list(df.values[:, smiles_col_id])
        return Molecules(names, smiles)

    def to_df(self, name_col: str = "name", smiles_col: str = "SMILES") -> DataFrame:
        df = pd.DataFrame(data={name_col: self.names, smiles_col: self.smiles})
        return df
