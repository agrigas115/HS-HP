# A stereochemically accurate all-atom protein model with hydrophobic interactions: Hard-sphere + hydrophobicity (HS+HP) protein model

This repository contains the necessary files to parameterize and run multiple versions the HS and HS+HP models under energy minimization and Langevin dynamics.

This project was initially developed to understand the jamming transition in folded proteins. Read more: https://arxiv.org/abs/2405.09646

## Adding hydrogens
Protein x-ray crystal structures are often deposited on the Protein Databank without all the appropriate hydrogens, as they can be difficult to resolve. Therefore, we need to appropriately add hydrogens to the structures before simulating them. Additionally, we need all the atom names to be consistent with what our parameter generation code expects. 

We use the Reduce software, as it is fairly standard in the protein structure field. A useful implementation of Reduce can be found in the larger protein structure software package Phenix. First, we trim all initial hydrogens from the PDB file PDBID.pdb, followed by adding all hydrogens, resulting in a structure file named PDBID_H.pdb
```
phenix.reduce -Trim -quiet PDBID.pdb > PDBID_noH.pdb

phenix.reduce -quiet PDBID_noH.pdb > PDBID_H.pdb
```
PDBID_H.pdb will be used to generate the simulation parameters and as the initial coordinates.

## Generating parameters

In order to run the simulation, we need to first identify all bonded and nonbonded interactions. This is handled by the python code included called get_params.py. The only non-standard library needed is Biopython: https://biopython.org. This library helps to parse the PDB file.

From the main directory, run the following command to generate the parameter files:
```
python PDBID_H.pdb
```
Example output for the PDBID 2f60 is included in the params directory.
