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


def run_sampling(args):
    topology = "Native.pdb"
    ticadim = 10
    num_sample_frames = 10000
    
    fn = args.file #file name
    wn = args.weights #weights name
    
    weights = np.loadtxt(wn)
    weights = weights/np.sum(weights)
    
    feat = coor.featurizer(topology)
    feat.add_distances_ca()
    X1 = coor.load(fn, feat, stride=1)
    
    sampled_frames=np.zeros((num_sample_frames,np.shape(X1)[1]))
    
    selected_frames = np.random.choice(np.shape(X1)[0], size=num_sample_frames, replace=True, p=weights)
    
    for i in range(num_sample_frames):
        ##debug
        #print np.shape(sampled_frames)
        #print np.shape(X1)
        ##debugg
        sampled_frames[i,:] = X1[selected_frames[i],:]
    
    ##debug
    for j in sampled_frames:
        for i in j:
            if i == 0:
                print "ERROR, distance too short, something not written"
    ##debugg
    
    tica_obj = coor.tica(sampled_frames, stride=1, lag=1, dim=ticadim)
    outputs = tica_obj.get_output()[0]
    eigen = tica_obj.eigenvalues
    print "saving files"
    np.savetxt("output.dat", outputs)
    np.savetxt("eigenvalues.dat", eigen)
    print "files saved"
        
    

def get_args():
    parser = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    parser.add_argument("--file", type=str, help="File, either .gro or .xtc")
    parser.add_argument("--weights", type=str, help="File, .w")
    args = parser.parse_args()
    
    return args
    

if __name__ == "__main__":
    args = get_args()
    run_sampling(args)



