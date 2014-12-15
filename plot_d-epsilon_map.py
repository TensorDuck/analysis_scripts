""" Plot the change in epislons of contacts map

Description: Takes two inputs, one for lower and one for upper triangle
Outputs a contact map based on the epsilons for a particular model

Justin Chen, November 17, 2014
"""


import numpy as np
import model_builder as mdb
import os as os
import matplotlib.pyplot as plt
import matplotlib.mpl as mpl

def plot_epsilons(top,bottom, namet, nameb):
    cwd = os.getcwd()
    fwd = "/home/jchen/projects/2014/10_31_1PB7/1PB7"

    os.chdir(fwd)
    modelA = mdb.check_inputs.load_model(fwd, False)
    os.chdir(cwd)
    
    pairs = modelA.contacts
    num_pairs = len(pairs)
    max_contact = np.max(pairs)
    
    colmax = np.max((np.max(np.max(np.abs(tophalf))), np.max(np.abs(bottomhalf))))    
    colmin = -1 * colmax
    
    x = np.arange(1, 294, 1)
    y = np.arange(1, 294, 1)
    square = np.zeros((max_contact,max_contact))
    
    for i in range(num_pairs):
        square[pairs[i][0]-1, pairs[i][1]-1] = bottom[i]
        square[pairs[i][1]-1, pairs[i][0]-1] = top[i]
    
    print np.shape(x)
    print np.shape(square)
    
    plt.pcolor(x, y, square ,cmap="RdBu", vmin=colmin, vmax=colmax)
    #plt.pcolor(x, y, square ,cmap="RdBu", vmin=-1, vmax=1)
    plt.xlabel("Residue i",fontsize=20)
    plt.ylabel("Residue j",fontsize=20)
    plt.title("Epsilon Comparison Between Top:"+namet+" and Bottom:"+nameb)
    plt.colorbar()
    plt.show()
    
    

if __name__ == "__main__":
    tophalf = np.loadtxt("epsilons_truncated_1PB7_135_0_trunc-1.dat")
    bottomhalf = np.loadtxt("epsilons_lambda_1PB7_135_0-l0.001.dat")
    
    
    plot_epsilons(tophalf, bottomhalf, "SVDT at 1", "SVDL at 0.001")

