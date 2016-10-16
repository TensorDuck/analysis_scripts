""" functions for analyzing a the Phi value 

From Best and Hummer PNAS 2016

Written for BPTI. Transition state is assumed to be sampled.
Calculates probabilities based off native contact probabilities

"""

import model_builder as mdb
import numpy as np
import mdtraj as md

def compute_phi_gradual(traj, relevant_native_potentials, periodic=False):
    #assumes traj is the correct traj
    #periodic is almost always False, unless for all-atom simulation
    native_pairs = []
    num_pairs = len(relevant_native_potentials)
    for pot in relevant_native_potentials:
        native_pairs.append([pot.atmi.index, pot.atmj.index])
    distances = md.compute_distances(traj, native_pairs, periodic=periodic)
    assert np.shape(distances)[1] == num_pairs
    
    total_prob = 0.
    for idx in range(num_pairs):
        prob = 0.
        r0 = relevant_native_potentials[idx].r0
        for jdx in range(traj.n_frames):
            if distances[jdx,idx] < r0:
                prob += 1.
            else:
                val = relevant_native_potentials[idx].dVdeps(distances[jdx,idx]) 
                prob -= val
        
        prob /= traj.n_frames
        total_prob += prob
    
    phi_value = total_prob / num_pairs
    return phi_value

def find_specific_native_potentials(phi_residue_idx, potentials_list, native_list):
    #assumes phi_residue_idx is already indexed the python way (-1)
    #potentials_list is a list of pairwise interactison from mdb
    #native_list has a list of native atom indices 
    
    affected_potentials = []
    for pot in potentials_list:
        if pot.atmi.residue.index == phi_residue_idx or pot.atmj.residue.index == phi_residue_idx:
            if check_idx_list(pot.atmi.index, pot.atmj.index, native_list):
                affected_potentials.append(pot)
    return affected_potentials

def check_idx_list(idx, jdx, native_list):
    #assumes idx, jdx and native_list are all python formatted (-1)
    found = False
    stop = False
    num_native = np.shape(native_list)[0]
    i = 0
    
    while not (found or stop):
        if native_list[i][0] == idx and native_list[i][1] == jdx:
            found = True
        elif native_list[i][1] == idx and native_list[i][0] == jdx: 
            found = True
        else:
            if i + 1 == num_native:
                stop = True
        i += 1  
           
    return found
