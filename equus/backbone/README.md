# Carbon backbone generation methods

This is a submodule designed to generate backbones via permutation of input structures, heavily relying on SMARTS.

## List of individual permutation methods
Reduction.  
Oxidation.  
Breaking single bonds.  
Replacing atoms.  
Replacing SELFIES tokens.  
Adding small groups. (This is not necessarily only for backbone generation.)

To convert backbones to diamines/diimines, please use equus.edit.generate.

## Important

It is probably not a good idea to use these methods to create a dataset for model training. The distribution of the randomly generated molecules has little meaning, merely reflecting the bias of the person who defines the pseudo_diffusion method. 

Instead, this submodule will be used as a baseline method.
