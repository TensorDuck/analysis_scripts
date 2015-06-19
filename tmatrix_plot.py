"""
Method for plotting transition matrices and comparing them with the data set
"""

import numpy as np
import os
import scipy.stats as stats
import matplotlib.pyplot as plt
import mdtraj as md

import fret_analysis.reader as reader
import analysis_scripts.pair_distance_calculator
import analysis_scripts.plot_package as plotter


def get_state_trace(y_data_set, ran_values, bins, weights=None):
    #ydata set is a 1-D arrray
    # set all the weights equal to one when None is specified
    if weights == None:
        weights = np.ones(np.shape(i)[0])
    
    statistic, bin_edges, slices = stats.binned_statistic(y_data_set, weights, statistic="sum", bins=bins, range=[ran_values])
    
    return statistic, bin_edges, state_traces


def get_T_matrix(num_states, state_trace, step=1, detailed_balance=False, slide=True):
    T_matrix = np.zeros((num_states, num_states))
    
    #get count matrix
    state_trace -= 1
    if slide:
        step_size = 1
    else:
        step_size = step
    
    for i in np.arange(0, np.shape(state_trace)[0]-step, step_size):
        #print states[i], states[i+1]
        T_matrix[states[i], states[i+1]] += 1
    
    if detailed_balance:
        import pyemma
        T_matrix_db = pyemma.msm.estimation.transition_matrix(T_matrix, reversible=True)
        T_masked = np.ma.masked_where(T_matrix_db == 0, T_matrix_db)
    else:
        #just normalize every row
        for i in range(np.shape(T_matrix)[0]):
            total_row = np.sum(T_matrix[i,:])
            if not total_row == 0:
                T_matrix[i,:] /= total_row
        T_masked = np.ma.masked_where(T_matrix == 0, T_matrix)
    return T_masked
    
    
    
def get_args():
    parser = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    parser.add_argument("--file", type=str, help="File, either .gro or .xtc")
    parser.add_argument("--top", type=str, help="topology file like a pdb")
    parser.add_argument("--savedir", type=str, default="%s/hist_plots"%os.getcwd(), help="Directory for saving the outputted data")
    parser.add_argument("--weights", type=str, help="File, .w")
    parser.add_argument("--save_name", type=str, default="plot", help="Name of file to be saved in")
    parser.add_argument("--pairs", nargs="+",type=int, default=[114,192], help="pairs for computing the y-distance")
    parser.add_argument("--spacing", type=float, default=0.1, help="spacing in nm for histogramming data")
    parser.add_argument("--y_shift", type=float, default=0, help="Specify the y-shift to the FRET distance data")
    parser.add_argument("--tmatrix_data", type=str, help="location of T-matrix data in 0-20 nm format")
    args = parser.parse_args()
    
    return args
    

if __name__ == "__main__":
    args = get_args()
    traj = md.load(args.file, top=args.top)
    distances = md.compute_distances(traj, args.pairs, periodic=False)
    
    
    
    
    
    
    
    
    
    
