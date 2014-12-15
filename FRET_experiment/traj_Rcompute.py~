""" Class reading the traj.xtc and Native.pdb data from 
a set of temperature runs

Contains a numpy array that contains all the binned histograms


If ran as main:

Output a binned data set of the distances between pairs into files
named as their temperature.txt
Output is in the form of x column: center of bins, y column: counts
All files are written into the directory "Histograms_by_temperature"

Distances are stored as nanometers

Assumptions:
Assumes you are in the folder Tf_0, if using __main__. 
self.pairs are a 2-d numpy array, 1st dimenison is pair-index, second dimension is a pair
        pair indices are also decreased by 1 from actual residue index (python starts at 0)

Arguments:
output_data = indices are (temperature, pair_index, array-histogram values)
trace_data = indices are (temperature, pair_index, array-trace)
"""

import numpy as np
import mdtraj as md
import os as os
from FRET_experiment.traj_R import traj_R

class traj_Rcompute(traj_R):
    def __init__(self, p, cd, bs, rs):
        self.pairs = p
        self.file_dir = cd
        self.bin_size = bs    #number of bins
        self.ran_size = rs    #range of values for the final distance
        step = (float(self.ran_size[1]-self.ran_size[0])) / (2*self.bin_size)
        start = self.ran_size[0] + step
        end = self.ran_size[1] - step 
        self.xspace = np.linspace(start, end, self.bin_size)
        self.hist_data = []
        self.trace_data = []
        self.write_dir = self.file_dir + "/" + "Histograms_by_temperature"
        self.temps = np.loadtxt(self.file_dir+"/"+"T_array.txt", str)
        self.ensure_dir(self.write_dir)
        self.print_once = 0
    
    ##Methods for storing the histogram data into a matrix 
    def write_all(self):
        for q in self.temps:
            for w in self.pairs:
                output_data = np.array([self.xspace,self.make_hist(self.file_dir+"/"+q, w)])
                output_data = np.transpose(output_data)
                self.save_files(self.write_dir+"/"+q+"_rp"+str(w[0]+1)+"-"+str(w[1]+1)+".txt", output_data)
                print "Made file for ", q
    
    def store_all(self):
        if self.hist_data == []:
            if np.shape(self.temps) == ():
                allpairs_data = []
                q = self.temps
                for w in self.pairs:
                    output_data = np.array([self.xspace,self.make_hist(self.file_dir+"/"+str(q), w)])
                    output_data = np.transpose(output_data)
                    allpairs_data.append(output_data)
                self.hist_data.append(allpairs_data)
                print "computed histogram data for file " + str(q)
            else:
                for q in self.temps:
                    allpairs_data = []
                    for w in self.pairs:
                        output_data = np.array([self.xspace,self.make_hist(self.file_dir+"/"+q, w)])
                        output_data = np.transpose(output_data)
                        allpairs_data.append(output_data)
                    self.hist_data.append(allpairs_data)
                    print "computed histogram data for file " + q
        else:
            print "This traj_Rcompute object already contains a histogram set"
    
    def read_all(self):
        if self.hist_data == []:
            for q in self.temps:
                allpairs_data = []
                for w in self.pairs:
                    output_data = np.loadtxt(self.write_dir+"/"+q+"_rp"+str(w[0]+1)+"-"+str(w[1]+1)+".txt")
                    allpairs_data.append(output_data)
                    print "read data from file " + q +" pair" + str(w)
                self.hist_data.append(allpairs_data)    
        else:
            print "This traj_Rcompute object already contains a histogram set"
                
    def make_hist(self,direc, pairs):
        traj = md.load(direc+"/"+"traj.xtc", top=direc+"/"+"Native.pdb")
        pair_distance = self.compute_distances(traj, pairs)
        hist, edges = np.histogram(pair_distance, bins=self.bin_size, range=self.ran_size)
        return hist
    
    ##Methods for handling reading and writing and storing traces    
    def write_trace_all(self):
        for q in self.temps:
            traj = md.load(self.file_dir+"/"+q+"/"+"traj.xtc", top=self.file_dir+"/"+q+"/"+"Native.pdb")
            for w in self.pairs:
                output_data = np.array([traj.time, self.compute_distances(traj, w)])
                output_data = np.transpose(output_data)
                self.save_files(self.write_dir+"/"+q+"_trace"+"_rp"+str(w[0]+1)+"-"+str(w[1]+1)+".txt", output_data)
            print "Made trace file for ", q
            
    def store_trace_all(self):
        if self.trace_data == []:
            for q in self.temps:
                traj = md.load(self.file_dir+"/"+q+"/"+"traj.xtc", top=self.file_dir+"/"+q+"/"+"Native.pdb")
                allpairs_data = []
                for w in self.pairs:
                    output_data = np.array([traj.time, self.compute_distances(traj, w)])
                    output_data = np.transpose(output_data)
                    allpairs_data.append(output_data)
                self.trace_data.append(allpairs_data)
                print "Stored trace data for ", q 
        else:
            print "This traj_Rcompute object already contains a trace set"
    
    def read_trace_all(self):
        if self.trace_data == []:
            count = 0
            for q in self.temps:
                allpairs_data = []
                for w in self.pairs:
                    output_data = np.loadtxt(self.write_dir+"/"+q+"_trace"+"_rp"+str(w[0]+1)+"-"+str(w[1]+1)+".txt")
                    allpairs_data.append(output_data)
                self.trace_data.append(output_data)
                print "read data from trace file " + q
        else:
            print "This traj_Rcompute object already contains a histogram set"
    
    ##Methods for dealing with files: making directories and saving files
    def ensure_dir(self, f):
        if not os.path.exists(f):
            os.makedirs(f)

    def save_files(self, filename, arrays):
        np.savetxt(filename, arrays, fmt=["%.10f","%.10f"])
    
    def compute_distances(self, traj, prs):
        return md.compute_distances(traj, np.array([prs]), periodic=False )
        #acord = traj.xyz[:,prs[0],:]
        #bcord = traj.xyz[:,prs[1],:] 
        #return (np.sum(((acord-bcord)**2), axis=1)) ** 0.5
        
        
        
        
        
        
if __name__ == "__main__":
    print "Hello World!"
    pairs = np.array([[115,193]])    #Atom pairs for computing distances between
    bin_size = 100.0    #number of bins
    ran_size = (0,10)    #range of values for the final distance
    directory = "Histograms_by_temperature"
    current_dir = os.getcwd()
    rmodel = traj_Rcompute(pairs, current_dir, bin_size, ran_size)
    rmodel.write_all()
    print current_dir
    print "Bye World!"           
            
            

