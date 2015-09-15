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

def plot_epsilons(bottom, nameb, modelA):
    
    pairs = modelA.pairs
    pairs = pairs[modelA.fitting_params]
    num_pairs = modelA.n_fitting_params
    max_contact = np.max(pairs)
    colmax = np.max(np.abs(bottom))
    colmin = -1 * colmax
    
    x = np.arange(1, np.max(pairs)+1, 1)
    y = np.arange(1, np.max(pairs)+1, 1)
    square = np.zeros((max_contact,max_contact))
    
    for i in range(num_pairs):
        square[pairs[i][0]-1, pairs[i][1]-1] = bottom[i]
        square[pairs[i][1]-1, pairs[i][0]-1] = colmax
    
    print "Starting ", nameb
    plt.figure()
    plt.pcolor(x, y, square ,cmap="RdBu", vmin=colmin, vmax=colmax)
    #plt.pcolor(x, y, square ,cmap="RdBu", vmin=-1, vmax=1)
    plt.xlabel("Residue i",fontsize=20)
    plt.ylabel("Residue j",fontsize=20)
    plt.title("1PB7 - Native Contacts and "+nameb)
    plt.colorbar()
    
    plt.savefig(nameb+".png")
    #plt.savefig(nameb+".pdf")
    plt.close()
    
def plot_epsilons_bin(bottom,nameb,model):
    bottomhalf = np.zeros(len(bottom))
    for i in range(len(bottom)):
        if bottom[i] > 0:
            bottomhalf[i] = 1
        elif bottom[i] < 0:
            bottomhalf[i] = -1

    plot_epsilons(bottomhalf, "bins-"+nameb,model)
    
    
if __name__ == "__main__":  
    cwd = os.getcwd()
    model = mdb.check_inputs.load_model(cwd, False)
    fid = cwd + "/Mut_0/newton"
    temp = np.loadtxt("%s/temp-used-here.txt" % fid)
    nameb = "deps-1PB7-%d" % temp    
    indx = np.loadtxt("%s/Lambda_index.txt" % fid)
    bottom = np.loadtxt("%s/xp_%d.dat" % (fid,indx))
    
    plot_epsilons_bin(bottom,nameb,model)
    plot_epsilons(bottom, nameb, model)

