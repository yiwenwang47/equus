from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from .edit import read_smiles
from .iso import pure_mol_to_nx, node_match
from networkx.algorithms import isomorphism


form = lambda mol: Chem.rdMolDescriptors.CalcMolFormula(mol)
n_bits = 2048
radius = 5
fpgen = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)

class Molecules(object):

    '''
    Class for keeping track of molecules.
    Has a search method.
    '''

    def __init__(self, names, smiles):

        self.names = names
        self.smiles = smiles
        self.mols = [read_smiles(smi, no_aromatic_flags=False, hydrogens=True) for smi in smiles]
        self.mol_form = [form(mol) for mol in self.mols]
        self.fps = [fpgen.GetFingerprint(mol) for mol in self.mols]

    def __getitem__(self, i) -> str:
        return self.smiles[i]

    def add(self, name, smi):

        self.names.append(name)
        self.smiles.append(smi)
        mol = read_smiles(smi, no_aromatic_flags=False, hydrogens=True)
        self.mols.append(mol)
        self.mol_form.append(form(mol))
        self.fps.append(fpgen.GetFingerprint(mol))

    def search(self, smi) -> list[bool, str]:

        '''
        Returns: 
        True, name if molecule is found.
        False, '' if not found.
        '''

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
            return False, ''
        
        # filter the list of candidates by Tanimoto Similarity of Morgan fingerprints
        candidates = [i for i, x in enumerate(self.mol_form) if x==formula]
        this_fp = fpgen.GetFingerprint(mol)
        candidates = [i for i in candidates if DataStructs.TanimotoSimilarity(self.fps[i], this_fp)>0.99]
        if len(candidates)==0:
            return False, ''
        
        this_graph = pure_mol_to_nx(mol)
        for i in candidates:
            ref_graph = pure_mol_to_nx(self.mols[i])
            if isomorphism.GraphMatcher(this_graph, ref_graph, node_match=node_match).is_isomorphic():
                return True, self.names[i]
            
        return False, ''
        
