from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import pandas as pd
from pandas.core.frame import DataFrame
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, rdmolops
from rdkit.Chem.rdchem import Mol


def mol_isomorphism(mol1: Mol, mol2: Mol, use_stereo: bool = False) -> bool:
    r"""
    This should be used as a last resort, as it is slow.
    Assumes the two molecules have the same molecular formula, the same Morgan fingerprint.
    To prioritize speed over absolute accuracy, this function follows the following steps:
    1. Check if the canonical SMILES strings are the same.
    2. Check if the InChI strings are the same.
    3. Check if the two molecules are subgraphs of each other.

    Note: It is unclear if two stereoisomers can have the same canonical SMILES string or InChI string.
    """

    canon_match = Chem.CanonSmiles(Chem.MolToSmiles(mol1)) == Chem.CanonSmiles(
        Chem.MolToSmiles(mol2)
    )
    if not canon_match:
        return False
    inchi_match = Chem.MolToInchi(mol1) == Chem.MolToInchi(mol2)
    if not inchi_match:
        return False
    return any(mol1.GetSubstructMatches(mol2, useChirality=use_stereo)) and any(
        mol2.GetSubstructMatches(mol1, useChirality=use_stereo)
    )


_form = lambda mol: Chem.rdMolDescriptors.CalcMolFormula(mol)
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
    use_stereo: bool = False

    def __post_init__(self):
        self.smiles = [Chem.CanonSmiles(smi) for smi in self.smiles]
        self.mols = [rdmolops.AddHs(Chem.MolFromSmiles(smi)) for smi in self.smiles]

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
        smi = Chem.CanonSmiles(smi)
        self.names.append(name)
        self.smiles.append(smi)
        mol = rdmolops.AddHs(Chem.MolFromSmiles(smi))
        self.mols.append(mol)
        self.fps.append(self.fp_generator.GetFingerprint(mol))
        formula, i = _form(mol), self.n
        self.mol_form_dict[formula].append(i)
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
        mol = rdmolops.AddHs(Chem.MolFromSmiles(smi))
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

        # Note: this is slightly (~15%) faster than the networkx approach
        for i in candidates:
            ref_mol = rdmolops.AddHs(Chem.MolFromSmiles(self.smiles[i]))
            if mol_isomorphism(mol, ref_mol, use_stereo=self.use_stereo):
                return True, self.names[i], self.smiles[i]

        return False, "", ""

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
