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


## Running the HS+HP model
We include a version of the model in which the initial structure is breifly energy minimized using FIRE energy minimization. This is followed by a constant temperature integration of Langevin dynamics. The coordinates will be periodically output in the .xyz format, which can be visualized using Ovito: https://www.ovito.org

First, compile the code.
```
g++ run_run_xtal_HSHP_lang_T0-6.c -O3 -o MD
```
No C/C++ libraries are needed to compile the code. The -O3 compilation flag is for optimization and gives a modest speedup. The default temperature is T/\epsilon = 1e-6. The temperature can be changed manually, but make sure to adjust the step size appropriately.

argv[1] : PDBID
argv[2] : alpha | the attractive range of the hydrophobic interaction
argv[3] : beta | the attractive depth of the hydrophopbic interaction

The values of alpha and beta result in interaction distances r_alpha and r_beta. Alpha and beta must be chosen so that r_beta < r_alpha. In general, alpha^2 * beta / T ~ 1 samples the region of jamming onset, as described in the paper.

Here is an example run:
```
./MD 4f60 1.5 1e-6
```
