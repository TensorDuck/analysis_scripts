"""Plot difference in transition matrices between two iterations"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import math
import argparse
import os

def plot_difference_matrix(args):

    subdir = args.subdir
    iter = args.iter

    cwd = os.getcwd()

    location_1 = "%s/iteration_%d"%(subdir, iter[0])
    os.chdir(location_1)
    Tmatrix_1 = np.loadtxt(args.files[0])
    if check_input(Tmatrix_1):
        Tmatrix_1 = unflatten_matrix(Tmatrix_1)

    os.chdir(cwd)

    location_2 = "%s/iteration_%d"%(subdir, iter[1])
    os.chdir(location_2)
    Tmatrix_2 = np.loadtxt(args.files[1])
    if check_input(Tmatrix_2):
        Tmatrix_2 = unflatten_matrix(Tmatrix_2)

    os.chdir(cwd)

    c_direc = "%s/%s/tmatrix_compare"%(cwd, subdir)
    
    if not os.path.isdir(c_direc):
        os.mkdir(c_direc)
    
    os.chdir(c_direc)
    
    diff_matrix = Tmatrix_1 - Tmatrix_2
    title = "Transition Difference Iteration %d and %d"%(iter[0], iter[1])
    plot_matrix(diff_matrix, title, -0.5, 0.5)
    
def unflatten_matrix(matrix):
    #Unflatten vector into square matrix

    dim = math.sqrt(np.shape(matrix)[0])

    T_matrix = np.zeros((dim, dim))

    matrix = np.transpose(matrix)
    for i in range(int(dim)):
        T_matrix[i,:] = matrix[dim*i:dim*(i+1)]

    return T_matrix

def check_input(matrix):
    # Check if matrix needs to be unflattened
    
    if np.shape(matrix)[0] == np.size(matrix):
        unflatten = True
    
    return unflatten

def plot_matrix(matrix, title, zmin, zmax):
    # Plot color map

    plt.figure()
    cmap = plt.pcolormesh(matrix, vmin=zmin, vmax=zmax)
    cbar = plt.colorbar(cmap)
    cbar.set_label("Transition Probability Difference")

    plt.axis([0, np.shape(matrix)[1], 0, np.shape(matrix)[0]])
    plt.xlabel("j")
    plt.ylabel("i")
    plt.title("%s"%title)

    plt.savefig("%s"%title)
    plt.close()

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("--subdir", type=str, help="Directory containing files")
    parser.add_argument("--iter", nargs=2, type=int, help="Iterations to compare")
    parser.add_argument("--files", default=("sim_feature.dat", "sim_feature.dat"), nargs=2, type=str, help="Tmatrix file names")

    args = parser.parse_args()

    return args

if __name__=="__main__":

    args = get_args()
    plot_difference_matrix(args)