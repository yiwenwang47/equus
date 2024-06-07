from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from networkx.algorithms import isomorphism
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from equus.edit.iso import node_match, pure_mol_to_nx
from equus.edit.utils import read_smiles

form = lambda mol: Chem.rdMolDescriptors.CalcMolFormula(mol)
n_bits = 2048
radius = 5
fpgen = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)


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
        self.mol_form = [form(mol) for mol in self.mols]
        self.fps = [fpgen.GetFingerprint(mol) for mol in self.mols]
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
        self.mol_form.append(form(mol))
        self.fps.append(fpgen.GetFingerprint(mol))
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
        formula = form(mol)
        if formula not in self.mol_form:
            return False, ""

        # filter the list of candidates by Tanimoto Similarity of Morgan fingerprints
        candidates = [i for i, x in enumerate(self.mol_form) if x == formula]
        this_fp = fpgen.GetFingerprint(mol)
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
        names = list(df.values[:, name_col_id])
        smiles = list(df.values[:, smiles_col_id])
        return Molecules[names, smiles]

    def to_csv(self, filename: str, name_col: str = "name", smiles_col: str = "SMILES"):
        df = pd.DataFrame(data={name_col: self.names, smiles_col: self.smiles})
        return df
