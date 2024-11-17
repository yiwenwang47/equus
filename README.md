# Effective Quest into Bounded Chemical Space (EQUUS)

This is a project focusing on rapid virtual screening of aromatic 1,2-diamines as candidates for Proton-Coupled Electron Transfer processes. More specifically, the goal is to tune the average bond dissociation free energy (BDFE) of the N-H bonds. This repository includes the code we used to generate diamines/diimines so it is quite rdkit-heavy.

<p align="center">
    <img src="scheme_colored_small.jpg" alt="Scheme" title="Scheme" width="300"/>
</p>

Data, code for BDFE evaluation will be made public soon.

## Installation

Please create a virtual environment by running the following.
```sh
conda create -n equus python=3.11
```

Install in editable mode:
```sh
cd to/this/directory
pip install -v -e .
```

Or
```sh
cd to/this/directory
pip install -r requirements.txt
pip install .
```
