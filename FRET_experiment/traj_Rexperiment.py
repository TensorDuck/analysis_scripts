"""A subclass of traj_R that will contain arrays of experimental data
along with the functions for getting and working with the experimental data
"""


import numpy as np
import mdtraj as md
import os as os
from FRET_experiment.traj_R import traj_R

class traj_Rexperiment(traj_R):
    def __init__(self, p, bs, rs):
        print "There are four lights"
        self.bin_size = bs    #number of bins
        self.ran_size = rs    #range of values for the final distance
        self.pairs = p
        step = (float(self.ran_size[1]-self.ran_size[0])) / (2*self.bin_size)
        start = self.ran_size[0] + step
        end = self.ran_size[1] - step 
        self.xspace = np.linspace(start, end, self.bin_size)
        print "Shut Up Wesley."
        self.hist_data = [] 
        self.labels = []
        self.print_once = 0  
          
    def read_hist(self, directory, l):
        allpairs = []
        allpairs.append(np.loadtxt(directory))
        self.hist_data.append(allpairs)
        self.labels.append(l)
        