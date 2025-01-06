# The Hard-sphere (HS) and Hard-sphere + hydrophobicity (HS+HP) protein models

This repository contains the necessary files to parameterize and run multiple versions the HS and HS+HP models under energy minimization and Langevin dynamics.

This project was initially developed to understand the jamming transition in folded proteins. Read more: https://arxiv.org/abs/2405.09646

## Adding hydrogens
Protein x-ray crystal structures are often deposited on the Protein Databank without all the appropriate hydrogens, as they can be difficult to resolve. Therefore, we need to appropriately add hydrogens to the structures before simulating them. Additionally, we need all the atom names to be consistent with what our parameter generation code expects. 

We use the Reduce software, as it is fairly standard in the protein structure field. A useful implementation of Reduce can be found in the larger protein structure software package Phenix. First, we trim all initial hydrogens from the PDB file, followed by adding all hydrogens, resulting in a structure file named PDBID_H.pdb
```
'phenix.reduce -Trim -quiet '+file+' > '+PDBID+'_noH.pdb'

'phenix.reduce -quiet '+PDBID+'_noH.pdb > '+PDBID+'_H.pdb'
```
PDBID_H.pdb will be used to generate the simulation parameters and as the initial coordinates.
## Generating parameters

