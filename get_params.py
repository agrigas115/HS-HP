#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate parameters for HS+HP simulation
argv[1] = PDBID.pdb
@author: agrigas115
"""
import glob
import numpy as np
import os
import sys
import json
from Bio.PDB import PDBParser

hydrophobicity_dict = {'ARG':0, 'ASP':0.09, 'GLU':0.16, 'LYS':0.16, 'ASN':0.25,
                       'GLN':0.29, 'PRO':0.39, 'HIS':0.4, 'SER':0.42, 'THR':0.48,
                       'GLY':0.52, 'TYR':0.64, 'ALA':0.67, 'CYS':0.74, 'MET':0.84,
                       'TRP':0.85, 'VAL':0.89, 'PHE':0.96, 'LEU':0.97, 'ILE':1.0}

p = PDBParser(PERMISSIVE=0)

N_list = []

pdb_file = sys.argv[1]

#try:
structure = p.get_structure('X', pdb_file)
atoms = structure.get_atoms()
atom_list = []
for atom in atoms:
    atom_list.append(atom)
    
resid_list = []
restype_list = []
resid_set = []
for model in structure:
    if len(model) == 1:
        for chain in model:
            for residue in chain:
                resid_set.append(residue._id[1])
                for atom in residue:
                    resid_list.append(residue._id[1])
                    restype_list.append(residue.resname)
        
        pdb_name = pdb_file.split('/')[-1].split('_H.pdb')[0]
        
        print(pdb_name)
        
        discont = 0
        for i in range(0,len(resid_set)-1):
            if resid_set[i]+1 != resid_set[i+1]:
                discont = 1
        
        if discont == 1:
            print('DISCONTINOUS!')
            
        
        if discont == 0:# and not os.path.exists('params/'+pdb_name+'_params/'):
            N = len(list(set(resid_list)))

            
            


            T = 1e-6
            
            coeff_bond = 1
            coeff_angle = 1
            coeff_omega = 1
            
            DTYPE = np.float64
            
            print('Loading Data...')
            
            box_size = 1000
            
            # Loading Mass Dictionary #
            file = './initialization/M_dict.txt'
            with open(file, 'r') as f:
                info = f.read()
            M_dict = json.loads(info)
            
            # Loading Atom Sizes #
            restype_sigma_dict = {}
            with open('./initialization/vdw_jennifer_dict.txt', 'r') as f:
                info = f.read()
            sourcelines = info.splitlines()
            sigma_dict = {}
            num_atoms_dict = {}
            num_atoms = 0
            restype = ''
            for line in sourcelines:
                if line.split()[0] == 'RESIDUE':
                    if len(sigma_dict) > 0:
                        restype_sigma_dict[restype] = sigma_dict
                        num_atoms_dict[restype] = num_atoms
                    restype = line.split()[2]
                    sigma_dict = {}
                    num_atoms = 0
                elif line.split()[0] == 'ATOM':
                    atom_type = line.split()[1]
                    sigma = float(line.split()[2])
                    sigma_dict[atom_type] = sigma
                    num_atoms += 1
            restype_sigma_dict[restype] = sigma_dict
            num_atoms_dict[restype] = num_atoms
            
            # Loading Bonded Dictionary #
            file = './initialization/bonded_dict.txt'
            with open(file, 'r') as f:
                info = f.read()
            bonded_dict = json.loads(info)
            
            # Loading Angle Dictionary #
            file = './initialization/angle_dict.txt'
            with open(file, 'r') as f:
                info = f.read()
            angle_dict = json.loads(info)
            
            def compute_ke(velocs, M):
                ke = ((M[:,0]/2.0) * np.square(np.linalg.norm(velocs,axis=1))).sum()#(velocs*velocs).sum(axis=1)).sum()#np.square(np.linalg.norm(velocs,axis=1)))
                return ke 
            
            # Loading Coordinates #
            print('Loading Coordinates...')
            with open(pdb_file, 'r') as f:
                info = f.read()
            sourcelines = info.splitlines()
            
            num_atoms = 0
            for i in range(0, len(sourcelines)):
                line = sourcelines[i]
                if line.split()[0] == 'ATOM':
                    num_atoms += 1
            
            atom_sequence = []
            resid_sequence = []
            restype_sequence = []
            epsilon_list = []
            coords = np.zeros((num_atoms,3), dtype=DTYPE)
            velocs =  np.zeros((num_atoms,3), dtype=DTYPE)
            M = np.zeros((num_atoms,3), dtype=DTYPE)
            num_residues = 0
            sequence_list = []
            prev_resid = -100
            count = 0
            for i in range(0, len(sourcelines)):
                line = sourcelines[i]
                if line.split()[0] == 'ATOM':
                    x = float(line.split()[6])
                    y = float(line.split()[7])
                    z = float(line.split()[8])
                    coords[count,:] = x, y, z
                    atom = line.split()[2]
                    if atom =='H1':
                        atom = 'H'
                    atom_sequence.append(atom)
                    resid = int(line.split()[5])
                    resid_sequence.append(resid)
                    restype = line.split()[3]
                    restype_sequence.append(restype)
                    epsilon = hydrophobicity_dict[restype]
                    epsilon_list.append(epsilon)
                    if resid != prev_resid:
                        num_residues += 1
                        sequence_list.append(restype)
                        prev_resid = resid
                        
                    M[count,0] = (M_dict[line.split()[2]])
                    M[count,1] = (M_dict[line.split()[2]])
                    M[count,2] = (M_dict[line.split()[2]])
                    count += 1
            
            # Centering in box #
            avg_x = np.average(coords[:,0])
            avg_y = np.average(coords[:,1])
            avg_z = np.average(coords[:,2])
            
            coords[:,0] -= (avg_x)# - (box_size/2))
            coords[:,1] -= (avg_y)# - (box_size/2))
            coords[:,2] -= (avg_z)# - (box_size/2))
            
            sigma_i_array = []
            for i in range(0, len(atom_sequence)):
                resid_i = resid_sequence[i]
                atom_i = atom_sequence[i]
                restype_i = restype_sequence[i]
                sigma_dict_i = restype_sigma_dict[restype_i]
                sigma_i_array.append(sigma_dict_i[atom_sequence[i]])
            
            # Draw initial velocities #
            # Removing Drift
            velocs = np.random.normal(-1,1,(num_atoms,3))
            velocs -= np.average(velocs,axis=0)
            # Removing rotation 
            len_origin = np.linalg.norm(coords, axis=1)[:, None]
            omega_velocs = np.cross(coords, velocs, axis=1) / np.square(len_origin)
            velocs -= np.cross(omega_velocs, coords)
            # Rescaling to T
            ke = compute_ke(velocs, M)
            temp = (2.0/3.0) * ke / num_atoms
            velocs *= np.sqrt(T/temp)
            
            # Loading bond length statistics #
            # Intra #
            print('Loading interaction stats...')
            bond_len_avg_dict = {}
            bond_len_k_dict = {}
            restypes = glob.glob('./initialization/bond_length_stats/intra/*')
            
            restype_bonded_dict = {}
            for res_dir in restypes:
                avg_dict = {}
                k_dict = {}
                bonded_list = []
                restype = res_dir.split('/')[-1]
                avg_file = res_dir+'/avg_bond_length_'+restype+'.txt'
                with open(avg_file, 'r') as f:
                    info = f.read()
                sourcelines = info.splitlines()
                for line in sourcelines:
                    bond = line.split()[0]
                    avg = float(line.split()[1])
                    avg_dict[bond] = avg
                    bonded_list.append(bond)
                
                std_file = res_dir+'/std_bond_length_'+restype+'.txt'
                with open(std_file, 'r') as f:
                    info = f.read()
                sourcelines = info.splitlines()
                for line in sourcelines:
                    bond = line.split()[0]
                    std = float(line.split()[1])
                    k = coeff_bond * T / (std**2)
                    k_dict[bond] = k
                
                bond_len_avg_dict[restype] = avg_dict
                bond_len_k_dict[restype] = k_dict
                restype_bonded_dict[restype] = bonded_list
            
            # Inter #
            with open('./initialization/bond_length_stats/inter/avg_bond_length_CN.txt') as f:
                info = f.read()
            sourcelines = info.splitlines()
            CN_bond_avg = float(sourcelines[0].split()[1])
            
            with open('./initialization/bond_length_stats/inter/std_bond_length_CN.txt') as f:
                info = f.read()
            sourcelines = info.splitlines()
            CN_bond_k = coeff_bond * T/(float(sourcelines[0].split()[1])**2)
            
            # Loading bond angle statistics #
            # Intra #
            bond_angle_avg_dict = {}
            bond_angle_k_dict = {}
            restypes = glob.glob('./initialization/bond_angle_stats/intra/*')
            #restype_bonded_dict = {}
            for res_dir in restypes:
                avg_dict = {}
                k_dict = {}
                angle_list = []
                restype = res_dir.split('/')[-1]
                avg_file = res_dir+'/avg_angle_'+restype+'.txt'
                with open(avg_file, 'r') as f:
                    info = f.read()
                sourcelines = info.splitlines()
                for line in sourcelines:
                    angle = line.split()[0]
                    avg = (float(line.split()[1])/180)*np.pi
                    avg_dict[angle] = avg
                    bonded_list.append(angle)
                
                std_file = res_dir+'/std_angle_'+restype+'.txt'
                with open(std_file, 'r') as f:
                    info = f.read()
                sourcelines = info.splitlines()
                for line in sourcelines:
                    angle = line.split()[0]
                    std = (float(line.split()[1])/180)*np.pi
                    k = coeff_angle * T / (std**2)
                    k_dict[angle] = k
                
                bond_angle_avg_dict[restype] = avg_dict
                bond_angle_k_dict[restype] = k_dict
                #restype_bonded_dict[restype] = bonded_list
            
            # Inter #
            inter_angle_avg_dict = {}
            inter_angle_k_dict = {}
            
            with open('./initialization/bond_angle_stats/inter/avg_angle_CNCA.txt') as f:
                info = f.read()
            inter_angle_avg_dict['CNCA'] = (float(info.splitlines()[0].split()[1])/180)*np.pi
            with open('./initialization/bond_angle_stats/inter/std_angle_CNCA.txt') as f:
                info = f.read()
            inter_angle_k_dict['CNCA'] = coeff_angle * T / ((float(info.splitlines()[0].split()[1])/180)*np.pi)**2
            
            with open('./initialization/bond_angle_stats/inter/avg_angle_CNH.txt') as f:
                info = f.read()
            inter_angle_avg_dict['CNH'] = (float(info.splitlines()[0].split()[1])/180)*np.pi
            with open('./initialization/bond_angle_stats/inter/std_angle_CNH.txt') as f:
                info = f.read()
            inter_angle_k_dict['CNH'] = coeff_angle * T / ((float(info.splitlines()[0].split()[1])/180)*np.pi)**2
            
            with open('./initialization/bond_angle_stats/inter/avg_angle_CACN.txt') as f:
                info = f.read()
            inter_angle_avg_dict['CACN'] = (float(info.splitlines()[0].split()[1])/180)*np.pi
            with open('./initialization/bond_angle_stats/inter/std_angle_CACN.txt') as f:
                info = f.read()
            inter_angle_k_dict['CACN'] = coeff_angle * T / ((float(info.splitlines()[0].split()[1])/180)*np.pi)**2
            
            with open('./initialization/bond_angle_stats/inter/avg_angle_OCN.txt') as f:
                info = f.read()
            inter_angle_avg_dict['OCN'] = (float(info.splitlines()[0].split()[1])/180)*np.pi
            with open('./initialization/bond_angle_stats/inter/std_angle_OCN.txt') as f:
                info = f.read()
            inter_angle_k_dict['OCN'] = coeff_angle * T / ((float(info.splitlines()[0].split()[1])/180)*np.pi)**2
            
            
            # Initializing bonds
            print('Setting up bond information...')
            nonbonded_list = []
            
            bond_length_indexes = []
            bond_length_indexes_avg = []
            bond_length_indexes_k = []
            va_bonded_index_array = []
            vb_bonded_index_array = []
            bonded_avg_array = []
            bonded_k_array = []
            
            va_nonbonded_index_array = []
            vb_nonbonded_index_array = []
            sigma_array = []
            epsilon_ij_array = []
            
            standard_bond_indexes_dict = {}
            bond_index_files = glob.glob('./initialization/standard_bond_indexes/*')
            for file in bond_index_files:
                name = file.split('/')[-1].split('_')[0]
                with open(file, 'r') as f:
                    info = f.read()
                sourcelines = info.splitlines()
                pair_res_list = []
                for line in sourcelines:
                    pair = [int(line.split()[0]), int(line.split()[1])]
                    pair_res_list.append(pair)
                standard_bond_indexes_dict[name] = pair_res_list
            
            print('\t bonds...')
            count = 0
            for k in range(0, num_residues):
                restype = restype_sequence[count+1]
                standard_bond_indexes = standard_bond_indexes_dict[restype]
                for pair in standard_bond_indexes:
                    i = pair[0]+count
                    j = pair[1]+count
                    atom_i = atom_sequence[i]
                    atom_j = atom_sequence[j]
                    bond_length_indexes.append([i,j])
                    va_bonded_index_array.append(i)
                    vb_bonded_index_array.append(j)
                    bonded_avg_array.append(bond_len_avg_dict[restype][atom_i+atom_j])
                    bonded_k_array.append(bond_len_k_dict[restype][atom_i+atom_j])
                count += num_atoms_dict[restype]
            
            count = 0 
            for i in range(0, len(atom_sequence)):
                atom_i = atom_sequence[i]
                if atom_i == 'C':
                    for j in range(i, len(atom_sequence)):
                        atom_j = atom_sequence[j]
                        if atom_j == 'N':
                            count += 1
                            va_bonded_index_array.append(i)
                            vb_bonded_index_array.append(j)
                            bonded_avg_array.append(CN_bond_avg)
                            bonded_k_array.append(CN_bond_k)
                            bond_length_indexes.append([i,j])
                            break
            
            bond_angle_indexes = []
            bond_angle_indexes_avg = []
            bond_angle_indexes_k = []
            bond_angle_lists = []
            va_angle_index_array = []
            vb_angle_index_array = []
            vc_angle_index_array = []
            angle_avg_array = []
            angle_k_array = []
            
            standard_angle_indexes_dict = {}
            angle_index_files = glob.glob('./initialization/standard_angle_indexes/*')
            for file in angle_index_files:
                name = file.split('/')[-1].split('_')[0]
                with open(file, 'r') as f:
                    info = f.read()
                sourcelines = info.splitlines()
                pair_res_list = []
                for line in sourcelines:
                    pair = [int(line.split()[0]), int(line.split()[1]), int(line.split()[2])]
                    pair_res_list.append(pair)
                standard_angle_indexes_dict[name] = pair_res_list
                
            
            print('\t angles...')
            count = 0
            for loop in range(0, num_residues):
                restype = restype_sequence[count+1]
                standard_angle_indexes = standard_angle_indexes_dict[restype]
                angle_count = 0
                for triple in standard_angle_indexes:
                    angle_count += 1
                    i = triple[0]+count
                    j = triple[1]+count
                    k = triple[2]+count
                    if k < 0:
                        k=0
                    atom_i = atom_sequence[i]
                    atom_j = atom_sequence[j]
                    atom_k = atom_sequence[k]
                    restype_i = restype_sequence[i]
                    restype_j = restype_sequence[j]
                    restype_k = restype_sequence[k]
                    if restype == 'PRO' and angle_count == 27:
                        while restype_k != sequence_list[loop-1] or atom_k != 'C':
                            k += 1
                            atom_k = atom_sequence[k]
                            restype_k = restype_sequence[k]
                    bond_angle_indexes.append([i, j, k])
                    va_angle_index_array.append(i)
                    vb_angle_index_array.append(j)
                    vc_angle_index_array.append(k)
                    angle_avg_array.append(bond_angle_avg_dict[restype_i][atom_i+atom_j+atom_k])
                    angle_k_array.append(bond_angle_k_dict[restype_i][atom_i+atom_j+atom_k])
                    bond_angle_indexes_avg.append(bond_angle_avg_dict[restype_i][atom_i+atom_j+atom_k])
                    bond_angle_indexes_k.append(bond_angle_k_dict[restype_i][atom_i+atom_j+atom_k])
                    bond_angle_lists.append([atom_i+atom_j+atom_k])
                    
                count += num_atoms_dict[restype]
             
            
            inter_residue_angles = [[[0,'CA'], [0,'C'], [1,'N']],
                                    [[0,'O'], [0,'C'], [1,'N']],
                                    [[0,'C'], [1,'N'], [1,'H']],
                                    [[0,'C'], [1,'N'], [1,'CA']]]
            
            got = 0
            for angle in inter_residue_angles:
                for i in range(0, len(atom_sequence)):
                    atom_i = atom_sequence[i]
                    resid_i = resid_sequence[i]
                    got = 0
                    if atom_i == angle[0][1]:
                        for j in range(0,len(atom_sequence)):
                            atom_j = atom_sequence[j]
                            resid_j = resid_sequence[j]
                            if atom_j == angle[1][1] and resid_j-resid_i == angle[1][0]:
                                for k in range(0,len(atom_sequence)):
                                    atom_k = atom_sequence[k]
                                    resid_k = resid_sequence[k]
                                    if atom_k == angle[2][1] and resid_k-resid_i == angle[2][0] and got == 0:
                                        got = 1
                                        bond_angle_indexes.append([i, j, k])
                                        bond_angle_indexes_avg.append(inter_angle_avg_dict[atom_i+atom_j+atom_k])
                                        bond_angle_indexes_k.append(inter_angle_avg_dict[atom_i+atom_j+atom_k])
                                        bond_angle_lists.append([atom_i+atom_j+atom_k])
                                        va_angle_index_array.append(i)
                                        vb_angle_index_array.append(j)
                                        vc_angle_index_array.append(k)
                                        angle_avg_array.append(inter_angle_avg_dict[atom_i+atom_j+atom_k])
                                        angle_k_array.append(inter_angle_k_dict[atom_i+atom_j+atom_k])
                                    
            
            va_angle_index_array = np.asarray(va_angle_index_array)
            vb_angle_index_array = np.asarray(vb_angle_index_array)
            vc_angle_index_array = np.asarray(vc_angle_index_array)
            angle_avg_array = np.reshape(np.asarray(angle_avg_array), (np.shape(angle_avg_array)[0],1))
            angle_k_array = np.reshape(np.asarray(angle_k_array), (np.shape(angle_avg_array)[0],1))
            
            # Sidechain Dihedrals
            print('\t planar Dihedrals...')
            planar_dihedral_indexes = []
            va_dihedral_index_array = []
            vb_dihedral_index_array = []
            vc_dihedral_index_array = []
            vd_dihedral_index_array = []
            standard_dihedral_indexes_dict = {}
            dihedral_index_files = glob.glob('./initialization/standard_sidechain_dihedral_indexes/*')
            for file in dihedral_index_files:
                name = file.split('/')[-1].split('_')[0]
                with open(file, 'r') as f:
                    info = f.read()
                sourcelines = info.splitlines()
                pair_res_list = []
                for line in sourcelines:
                    pair = [int(line.split()[0]), int(line.split()[1]), int(line.split()[2]), int(line.split()[3])]
                    pair_res_list.append(pair)
                standard_dihedral_indexes_dict[name] = pair_res_list
            
            count = 0
            for loop in range(0, num_residues):
                restype = restype_sequence[count+1]
                try:
                    standard_dihedral_indexes = standard_dihedral_indexes_dict[restype]
                    for quad in standard_dihedral_indexes:
                        i = quad[0]+count
                        j = quad[1]+count
                        k = quad[2]+count
                        l = quad[3]+count
                        atom_i = atom_sequence[i]
                        atom_j = atom_sequence[j]
                        atom_k = atom_sequence[k]
                        atom_l = atom_sequence[k]
                        restype_i = restype_sequence[i]
                        restype_j = restype_sequence[j]
                        restype_k = restype_sequence[k]
                        restype_l = restype_sequence[k]
                        planar_dihedral_indexes.append([i, j, k, l])
                        va_dihedral_index_array.append(i)
                        vb_dihedral_index_array.append(j)
                        vc_dihedral_index_array.append(k)
                        vd_dihedral_index_array.append(l)
                    count += num_atoms_dict[restype]
                except KeyError:
                    count += num_atoms_dict[restype]
                    pass
            
            va_dihedral_index_array = np.asarray(va_dihedral_index_array)
            vb_dihedral_index_array = np.asarray(vb_dihedral_index_array)
            vc_dihedral_index_array = np.asarray(vc_dihedral_index_array)
            vd_dihedral_index_array = np.asarray(vd_dihedral_index_array)
            
            print('\t improper Dihedrals...')
            improper_dihedral_indexes = []
            va_improper_index_array = []
            vb_improper_index_array = []
            vc_improper_index_array = []
            vd_improper_index_array = []
            standard_improper_indexes_dict = {}
            improper_index_files = glob.glob('./initialization/standard_sidechain_improper_indexes/*')
            for file in improper_index_files:
                name = file.split('/')[-1].split('_')[0]
                with open(file, 'r') as f:
                    info = f.read()
                sourcelines = info.splitlines()
                pair_res_list = []
                for line in sourcelines:
                    pair = [int(line.split()[0]), int(line.split()[1]), int(line.split()[2]), int(line.split()[3])]
                    pair_res_list.append(pair)
                standard_improper_indexes_dict[name] = pair_res_list
            
            count = 0
            for loop in range(0, num_residues):
                restype = restype_sequence[count+1]
                try:
                    standard_improper_indexes = standard_improper_indexes_dict[restype]
                    for quad in standard_improper_indexes:
                        i = quad[0]+count
                        j = quad[1]+count
                        k = quad[2]+count
                        l = quad[3]+count
                        atom_i = atom_sequence[i]
                        atom_j = atom_sequence[j]
                        atom_k = atom_sequence[k]
                        atom_l = atom_sequence[k]
                        restype_i = restype_sequence[i]
                        restype_j = restype_sequence[j]
                        restype_k = restype_sequence[k]
                        restype_l = restype_sequence[k]
                        improper_dihedral_indexes.append([i, j, k, l])
                        va_improper_index_array.append(i)
                        vb_improper_index_array.append(j)
                        vc_improper_index_array.append(k)
                        vd_improper_index_array.append(l)
                    count += num_atoms_dict[restype]
                except KeyError:
                    count += num_atoms_dict[restype]
                    pass
            
            va_improper_index_array = np.asarray(va_improper_index_array)
            vb_improper_index_array = np.asarray(vb_improper_index_array)
            vc_improper_index_array = np.asarray(vc_improper_index_array)
            vd_improper_index_array = np.asarray(vd_improper_index_array)
            
            # Get nonbonded atoms
            print('\t nonbonded...')
            print('\t inter-residue...')
            nonbonded_atoms = []
            num_inter = 0
            for i in range(0, len(atom_sequence)):
                resid_i = resid_sequence[i]
                atom_i = atom_sequence[i]
                epsilon_i = epsilon_list[i]
                for j in range(i+1, len(atom_sequence)):
                    resid_j = resid_sequence[j]
                    atom_j = atom_sequence[j]
                    epsilon_j = epsilon_list[j]
                    if i != j and resid_j > resid_i+1:
                        restype_i = restype_sequence[i]
                        restype_j = restype_sequence[j]
                        sigma_dict_i = restype_sigma_dict[restype_i]
                        sigma_dict_j = restype_sigma_dict[restype_j]
                        nonbonded_list.append([i,j])
                        nonbonded_atoms.append(atom_sequence[i]+atom_sequence[j])
                        va_nonbonded_index_array.append(i)
                        vb_nonbonded_index_array.append(j)
                        sigma_array.append(sigma_dict_i[atom_sequence[i]] + sigma_dict_j[atom_sequence[j]])
                        num_inter += 1
                        epsilon_ij = 0.5*(epsilon_i+epsilon_j)
                        epsilon_ij_array.append(epsilon_ij)
            print('\t intra-residue...')
            count = 0
            num_intra = 0
            for res in range(0, num_residues):
                restype = restype_sequence[count]
                num = num_atoms_dict[restype]
                for i in range(0, num):
                    index_i = i + count
                    atom_i = atom_sequence[index_i]
                    epsilon_i = epsilon_list[i]
                    for j in range(i+1, num):
                        index_j = j + count
                        atom_j = atom_sequence[index_j]
                        epsilon_j = epsilon_list[j]
                        standard_angle_indexes = standard_angle_indexes_dict[restype]
                        in_angle = 0
                        in_improper = 0
                        in_proper = 0
                        for angle in standard_angle_indexes:
                            if i in angle and j in angle:
                                in_angle = 1
                        if in_angle == 0:
                            try:
                                standard_dihedral_indexes = standard_dihedral_indexes_dict[restype]
                                for quad in standard_dihedral_indexes:
                                    if i in quad and j in  quad:
                                        in_proper = 1
                                if in_proper == 0:
                                    standard_improper_indexes = standard_improper_indexes_dict[restype]
                                    for quad in standard_improper_indexes:
                                        if i in quad and j in  quad:
                                            in_improper = 1
                            except KeyError:
                                pass
                        if in_angle == 0 and in_proper == 0 and in_improper == 0:
                            sigma_dict = restype_sigma_dict[restype]
                            nonbonded_list.append([index_i,index_j])
                            nonbonded_atoms.append(atom_j+atom_j)
                            va_nonbonded_index_array.append(index_i)
                            vb_nonbonded_index_array.append(index_j)
                            sigma_array.append(sigma_dict[atom_i] + sigma_dict[atom_j])
                            epsilon_ij = 0.5*(epsilon_i+epsilon_j)
                            epsilon_ij_array.append(epsilon_ij)
                            num_intra += 1
                count += num
                
            
            # nonbonded inter-residue by 1
            count_i = 0
            count_j = 0
            num_inter_1 = 0
            for res_i in range(0, num_residues-1):
                restype_i = restype_sequence[count_i]
                num_i = num_atoms_dict[restype_i]
                count_j += num_i
                res_j = res_i+1
                restype_j = restype_sequence[count_j]
                num_j = num_atoms_dict[restype_j]
                for i in range(0, num_i):
                    index_i = i + count_i
                    atom_i = atom_sequence[index_i]
                    epsilon_i = epsilon_list[i]
                    pair_i = [0, atom_i]
                    for j in range(0, num_j):
                        index_j = j + count_j
                        atom_j = atom_sequence[index_j]
                        epsilon_j = epsilon_list[j]
                        pair_j = [1, atom_j]
                        in_angle = 0
                        for angle in inter_residue_angles:
                            if pair_i in angle and pair_j in angle:
                                in_angle = 1
                        if in_angle == 0 and [atom_i, atom_j] != ['CA','CA'] and [atom_i, atom_j] != ['CA','H'] and [atom_i, atom_j] != ['O','H'] and [atom_i, atom_j] != ['O','CA']:
                            sigma_dict_i = restype_sigma_dict[restype_i]
                            sigma_dict_j = restype_sigma_dict[restype_j]
                            nonbonded_list.append([index_i,index_j])
                            nonbonded_atoms.append(atom_sequence[index_i]+atom_sequence[index_j])
                            va_nonbonded_index_array.append(index_i)
                            vb_nonbonded_index_array.append(index_j)
                            sigma_array.append(sigma_dict_i[atom_i] + sigma_dict_j[atom_j])
                            epsilon_ij = 0.5*(epsilon_i+epsilon_j)
                            epsilon_ij_array.append(epsilon_ij)
                            num_inter_1 += 1
                count_i += num_i
                
            
            va_bonded_index_array = np.asarray(va_bonded_index_array)
            vb_bonded_index_array = np.asarray(vb_bonded_index_array)
            bonded_avg_array = np.asarray(bonded_avg_array)
            bonded_k_array = np.asarray(bonded_k_array)        
            sigma_array = np.asarray(sigma_array)
            sigma_array = np.reshape(sigma_array,(np.shape(va_nonbonded_index_array)[0],1))
            
            va_nonbonded_index_array = np.asarray(va_nonbonded_index_array)
            vb_nonbonded_index_array = np.asarray(vb_nonbonded_index_array)
            
            #Initialize Dihedral omega
            print('\t dihedral omega...')
            CA1_omega_index_array = []
            C_omega_index_array = []
            N_omega_index_array = []
            CA2_omega_index_array = []
            for i in range(0, len(atom_sequence)):
                atom_i = atom_sequence[i]
                complete = 0
                if atom_i == 'CA':
                    CA1_index = i
                    CA1_resid = resid_sequence[i]
                    C_got = 0
                    N_got = 0
                    CA2_got = 0
                    for j in range (i+1, len(atom_sequence)):
                        atom_j = atom_sequence[j]
                        if atom_j == 'C' and C_got == 0:
                            C_index = j
                            C_resid = resid_sequence[j]
                            complete += 1
                            C_got = 1
                        if atom_j == 'N' and N_got == 0:
                            N_index = j
                            N_resid = resid_sequence[j] 
                            complete += 1
                            N_got = 1
                        if atom_j == 'CA' and CA2_got == 0:
                            CA2_index = j
                            CA2_resid = resid_sequence[j]
                            complete += 1
                            CA2_got = 1
                        if complete == 3:
                            continue
                                
                    try: 
                        if CA1_resid == C_resid == N_resid-1 == CA2_resid-1 and complete == 3:
                            CA1_omega_index_array.append(CA1_index)
                            C_omega_index_array.append(C_index)
                            N_omega_index_array.append(N_index)
                            CA2_omega_index_array.append(CA2_index)
                    except NameError:
                        pass
                        
            CA1_omega_index_array = np.asarray(CA1_omega_index_array)
            C_omega_index_array = np.asarray(C_omega_index_array)
            N_omega_index_array = np.asarray(N_omega_index_array)
            CA2_omega_index_array = np.asarray(CA2_omega_index_array)    
            
            with open('./initialization/dihedral_omega_stats/avg_omega.txt', 'r') as f:
                info = f.read()
            sourcelines = info.splitlines()
            avg_omega = float(sourcelines[0].split()[1])
            print(np.degrees(avg_omega))
            with open('./initialization/dihedral_omega_stats/std_omega.txt', 'r') as f:
                info = f.read()
            sourcelines = info.splitlines()
            std_omega = float(sourcelines[0].split()[1])
            print(np.degrees(std_omega))
            k_omega = coeff_omega * T / float(sourcelines[0].split()[1])**2
            
            
            sigma = np.zeros((len(atom_sequence),len(atom_sequence)), dtype=DTYPE)
            for i in range(0, len(atom_sequence)):
                for j in range(0,len(atom_sequence)):
                    restype_i = restype_sequence[i]
                    restype_j = restype_sequence[j]
                    sigma_dict_i = restype_sigma_dict[restype_i]
                    sigma_dict_j = restype_sigma_dict[restype_j]
                    sigma[i,j] = sigma_dict_i[atom_sequence[i]] + sigma_dict_j[atom_sequence[j]]

            
            #%%
            ####   Writing Data   ####
            
            param_dir = 'params/'+pdb_name+'_params/'
            
            command = 'mkdir '+param_dir
            os.system(command)
            
            with open(param_dir+'coords.txt', 'w') as f:
                for i in range(0, np.shape(coords)[0]):
                    for j in range(0, np.shape(coords)[1]):
                        f.write(str(coords[i,j])+'\t')
                    f.write('\n')
            
            with open(param_dir+'velocs.txt', 'w') as f:
                for i in range(0, np.shape(velocs)[0]):
                    for j in range(0, np.shape(velocs)[1]):
                        f.write(str(velocs[i,j])+'\t')
                    f.write('\n')
            
            with open(param_dir+'sigma_i_array.txt', 'w') as f:
                for i in range(0, len(sigma_i_array)):
                    f.write(str(sigma_i_array[i]))
                    f.write('\n')
            
            with open(param_dir+'sigma_array.txt', 'w') as f:
                for i in range(0, np.shape(sigma_array)[0]):
                    f.write(str(sigma_array[i,0]))
                    f.write('\n')
            
            with open(param_dir+'nonbonded_array1.txt', 'w') as f:
                for i in range(0, np.shape(va_nonbonded_index_array)[0]):
                    f.write(str(va_nonbonded_index_array[i]))
                    f.write('\n')
            
            with open(param_dir+'nonbonded_array2.txt', 'w') as f:
                for i in range(0, np.shape(vb_nonbonded_index_array)[0]):
                    f.write(str(vb_nonbonded_index_array[i]))
                    f.write('\n')
                    
            with open(param_dir+'atom_array.txt', 'w') as f:
                for i in range(0, len(atom_sequence)):
                    f.write(atom_sequence[i])
                    f.write('\n')
            
            with open(param_dir+'bonded_array1.txt', 'w') as f:
                for i in range(0, np.shape(va_bonded_index_array)[0]):
                    f.write(str(va_bonded_index_array[i]))
                    f.write('\n')
              
            with open(param_dir+'bonded_array2.txt', 'w') as f:
                for i in range(0, np.shape(vb_bonded_index_array)[0]):
                    f.write(str(vb_bonded_index_array[i]))
                    f.write('\n')  
              
            with open(param_dir+'bonded_avg_array.txt', 'w') as f:
                for i in range(0, np.shape(bonded_avg_array)[0]):
                    f.write(str(bonded_avg_array[i]))
                    f.write('\n')  
            
            with open(param_dir+'bonded_k_array.txt', 'w') as f:
                for i in range(0, np.shape(bonded_k_array)[0]):
                    f.write(str(bonded_k_array[i]))
                    f.write('\n')  
            
            with open(param_dir+'angle_array1.txt', 'w') as f:
                for i in range(0, np.shape(va_angle_index_array)[0]):
                    f.write(str(va_angle_index_array[i]))
                    f.write('\n')
                    
            with open(param_dir+'angle_array2.txt', 'w') as f:
                for i in range(0, np.shape(vb_angle_index_array)[0]):
                    f.write(str(vb_angle_index_array[i]))
                    f.write('\n')
             
            with open(param_dir+'angle_array3.txt', 'w') as f:
                for i in range(0, np.shape(vc_angle_index_array)[0]):
                    f.write(str(vc_angle_index_array[i]))
                    f.write('\n')
                    
            with open(param_dir+'angle_avg_array.txt', 'w') as f:
                for i in range(0, np.shape(angle_avg_array)[0]):
                    f.write(str(angle_avg_array[i][0]))
                    f.write('\n')
            
            with open(param_dir+'angle_k_array.txt', 'w') as f:
                for i in range(0, np.shape(angle_k_array)[0]):
                    f.write(str(angle_k_array[i][0]))
                    f.write('\n')
            
            with open(param_dir+'CA1_omega_array.txt', 'w') as f:
                for i in range(0, np.shape(CA1_omega_index_array)[0]):
                    f.write(str(CA1_omega_index_array[i]))
                    f.write('\n')
            
            with open(param_dir+'C_omega_array.txt', 'w') as f:
                for i in range(0, np.shape(C_omega_index_array)[0]):
                    f.write(str(C_omega_index_array[i]))
                    f.write('\n')
            
            with open(param_dir+'N_omega_array.txt', 'w') as f:
                for i in range(0, np.shape(N_omega_index_array)[0]):
                    f.write(str(N_omega_index_array[i]))
                    f.write('\n')
            
            with open(param_dir+'CA2_omega_array.txt', 'w') as f:
                for i in range(0, np.shape(CA2_omega_index_array)[0]):
                    f.write(str(CA2_omega_index_array[i]))
                    f.write('\n')
                    
            with open(param_dir+'avg_omega.txt', 'w') as f:
                f.write(str(avg_omega))
            
            with open(param_dir+'k_omega.txt', 'w') as f:
                f.write(str(k_omega))
            
            with open(param_dir+'box_size.txt', 'w') as f:
                f.write('500')
            
            with open(param_dir+'A_planar_array.txt', 'w') as f:
                for i in range(0, np.shape(va_dihedral_index_array)[0]):
                    f.write(str(va_dihedral_index_array[i]))
                    f.write('\n')
            
            with open(param_dir+'B_planar_array.txt', 'w') as f:
                for i in range(0, np.shape(vb_dihedral_index_array)[0]):
                    f.write(str(vb_dihedral_index_array[i]))
                    f.write('\n')
                
            with open(param_dir+'C_planar_array.txt', 'w') as f:
                for i in range(0, np.shape(vc_dihedral_index_array)[0]):
                    f.write(str(vc_dihedral_index_array[i]))
                    f.write('\n')
            
            with open(param_dir+'D_planar_array.txt', 'w') as f:
                for i in range(0, np.shape(vd_dihedral_index_array)[0]):
                    f.write(str(vd_dihedral_index_array[i]))
                    f.write('\n')
            
            with open(param_dir+'A_improper_array.txt', 'w') as f:
                for i in range(0, np.shape(va_improper_index_array)[0]):
                    f.write(str(va_improper_index_array[i]))
                    f.write('\n')
            
            with open(param_dir+'B_improper_array.txt', 'w') as f:
                for i in range(0, np.shape(vb_improper_index_array)[0]):
                    f.write(str(vb_improper_index_array[i]))
                    f.write('\n')
                
            with open(param_dir+'C_improper_array.txt', 'w') as f:
                for i in range(0, np.shape(vc_improper_index_array)[0]):
                    f.write(str(vc_improper_index_array[i]))
                    f.write('\n')
            
            with open(param_dir+'D_improper_array.txt', 'w') as f:
                for i in range(0, np.shape(vd_improper_index_array)[0]):
                    f.write(str(vd_improper_index_array[i]))
                    f.write('\n')
                    
            with open(param_dir+'inter_intra_idx.txt','w') as f:
                f.write(str(num_inter))
                f.write('\n')
                f.write(str(num_intra))
                f.write('\n')
                f.write(str(num_inter_1))
                f.write('\n')
            
            with open(param_dir+'resid_array.txt', 'w') as f:
                for i in range(0, len(resid_list)):
                    f.write(str(resid_list[i]))
                    f.write('\n')
            
            with open(param_dir+'epsilon_ij_array.txt', 'w') as f:
                for i in range(0, len(epsilon_ij_array)):
                    f.write(str(epsilon_ij_array[i]))
                    f.write('\n')
                    
    else:
        print('Multimer')
