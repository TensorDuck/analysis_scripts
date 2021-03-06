""" traj_trace.py 
description: read in traj.xtc, output a trace of the FRET-probe distance

"""

import mdtraj as md
import numpy as np
import os as os
import matplotlib.pyplot as plt

def calc_trace(cwd, pairs):
    print "Starting calc_trace"
    traj = md.load(cwd+"/traj.xtc", top=cwd+"/Native.pdb")
    pair_distance = md.compute_distances(traj, pairs, periodic=False)
    numpairs = np.shape(pair_distance)[1]
    colors = ["k", "b", "r", "m", "c", "y"]
    
    print "Start plotting"
    plt.figure()
    for i in range(numpairs):
        xspace = np.arange(0.0,np.shape(pair_distance[:,i])[0],1.0)
        xspace *= 0.0005
        plt.plot(xspace, pair_distance[:,i], color=colors[i])
    
    plt.xlabel("time (ns)", fontsize=20)
    plt.ylabel("pair distance (nm)", fontsize=20)
    #plt.axis([0, xspace[len(xspace)-1], 0, 10])
    plt.savefig("%dns-run.png"%np.max(xspace))
    print "Finished Plotting"

if __name__ == "__main__":
    cwd = os.getcwd()
    pairs = [[114,192]]
    calc_trace(cwd, pairs)




