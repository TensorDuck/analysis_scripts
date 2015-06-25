""" Script for calculating Tica from a weighted set of points
Does so by taking a random sampling of the points
"""

import numpy as np
try:
    import pyemma
    import pyemma.coordinates as coor
except:
    print "pyemma not loaded, tica methods will fail!"
import mdtraj as md
import time
import analysis_scripts.plot_package as pltpkg
import argparse
import time
import os

def run_sampling(args):
    topology = args.topfile
    ticadim = 10
    num_sample_frames = 10
    tica_lag_time = 50
    fn = args.filedir #file name
    wn = args.weights #weights name
    
    weights = np.loadtxt(wn)
    weights = weights/np.sum(weights)
    #first time
    time1 = time.clock()
    feat = coor.featurizer(topology)
    feat.add_distances_ca()
    selected_frames = np.random.choice(args.number_traj, size=num_sample_frames, replace=True, p=weights)
    
    selected_files = []
    selected_frames.sort()
    for i in selected_frames:
        selected_files.append("%s/traj%d.xtc"%(fn,i))
    time2 = time.clock()
    print "Took %f minutes to select new frames" % ((time2-time1)/60.0)
    sampled_frames = coor.load(selected_files, feat, stride=1)
    
    time3 = time.clock()
    print "Took %f minutes to load the new frames" % ((time3-time2)/60.0)
    
    tica_obj = coor.tica(sampled_frames, stride=1, lag=tica_lag_time, dim=ticadim)
    time4 = time.clock()
    print "Took %f minutes to calculate the tica_object" % ((time4-time3)/60.0)
    all_outputs = tica_obj.get_output()[0]
    for i in xrange(num_sample_frames-1):
        outputs = tica_obj.get_output()[i+1]
        all_outputs = np.append(all_outputs, outputs, axis=0)
    eigen = tica_obj.eigenvalues
    print "saving files"
    np.savetxt("output.dat", outputs)
    np.savetxt("eigenvalues.dat", eigen)
    print "files saved"
    time5 = time.clock()
    print "Took %f minutes to write the output files" % ((time5-time4)/60.0)
        
    

def get_args():
    parser = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    parser.add_argument("--filedir", default=os.getcwd(), type=str, help="Location of all xtc files")
    parser.add_argument("--number_traj", type=int, help="number of trajectories to look at")
    parser.add_argument("--topfile", type=str, help="top file for analysis")
    parser.add_argument("--weights", type=str, help="File, .w")
    args = parser.parse_args()
    
    return args
    

if __name__ == "__main__":
    args = get_args()
    run_sampling(args)



