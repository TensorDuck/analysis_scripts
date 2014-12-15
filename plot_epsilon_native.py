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
    
    pairs = modelA.contacts
    num_pairs = len(pairs)
    max_contact = np.max(pairs)
    
    colmax = np.max(np.abs(bottom))
    colmin = -1 * colmax
    
    x = np.arange(1, 294, 1)
    y = np.arange(1, 294, 1)
    square = np.zeros((max_contact,max_contact))
    
    for i in range(num_pairs):
        square[pairs[i][0]-1, pairs[i][1]-1] = bottom[i]
        square[pairs[i][1]-1, pairs[i][0]-1] = colmax
    
    print np.shape(x)
    print np.shape(square)
    print "Starting ", nameb
    plt.figure()
    plt.pcolor(x, y, square ,cmap="RdBu", vmin=colmin, vmax=colmax)
    #plt.pcolor(x, y, square ,cmap="RdBu", vmin=-1, vmax=1)
    plt.xlabel("Residue i",fontsize=20)
    plt.ylabel("Residue j",fontsize=20)
    plt.title("1PB7 - Native Contacts and "+nameb)
    plt.colorbar()
    
    plt.savefig(nameb+".png")
    plt.savefig(nameb+".pdf")
    plt.close()
    
def plot_epsilons_bin(bottom,nameb):
    for i in range(len(bottomhalf)):
        if bottomhalf[i] > 0:
            bottomhalf[i] = 1
        elif bottomhalf[i] < 0:
            bottomhalf[i] = -1

    plot_epsilons(bottomhalf, "_bins-"+nameb)
    
    
if __name__ == "__main__":
'''
    #truncation = ["100", "1", "0.1", "0.01", "0.001"]
    truncation = ["0.01"]
    for i in truncation:
        bottomhalf = np.loadtxt("epsilons_truncated_1PB7_135_0_trunc-"+i+".dat")
        plot_epsilons(bottomhalf, "SVDT:"+i)

    lavenberg = ["20", "50", "0.01", "0.001"]
    for i in lavenberg:
        bottomhalf = np.loadtxt("epsilons_lambda_1PB7_135_0-l"+i+".dat")
        plot_epsilons_bin(bottomhalf, "SVDL:"+i)
 '''   
    cwd = os.getcwd()
    model = mdb.check_inputs.load_model(cwd, False)
    fid = cwd + "/Mut_0/newton"
    temp = np.loadtxt("%s/temp-used-here.txt" % fid)
    nameb = "deps-1PB7-%d" % temp    
    indx = np.loadtxt("%s/lambda_index.txt" % fid)
    bottom = np.loadtxt("%s/xp_%d.dat" % (fid,indx))
    
    plot_epsilons_bin(bottom,nameb,model)
    plot_epsilons(bottom, nameb, model)
    

