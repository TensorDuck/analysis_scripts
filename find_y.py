"""
Calculate the y c-alpha distance for a variety of data types
 
"""

import numpy as np
import matplotlib.pyplot as plt
import os as os
import argparse
import mdtraj as md

import analysis_scripts.gro_reader as gread
import analysis_scripts.pair_distance_calculator as pdistance

class traj(object):
    def __init__(self, xyz):
        self.xyz = xyz
        self.n_atoms = np.shape(xyz)[1]

def make_gro(args):
    distances = pdistance.compute_distances(traj(gread.read(args.file)), args.pairs)
    for i in range(np.shape(args.pairs)[0]):
        np.savetxt("%s/%s-y.out"%(args.savedir, args.save_name), distances[:,i])

def make_xtc(args):
    distances = pdistance.compute_distances(md.load(args.file, top=args.native), args.pairs)
    for i in range(np.shape(args.pairs)[0]):
        np.savetxt("%s/%s-y.out"%(args.savedir, args.save_name), distances[:,i])

def sanitize_args(args):
    ##set the pairs for fitting in the array format for the calculation
    pairs = np.array([[args.pairs[0], args.pairs[1]]])
    print "Number of pairs is: %d" % len(args.pairs)
    if len(args.pairs)>2:
        for i in np.arange(3, len(args.pairs), 2):
            pairs = np.append(pairs, np.array([[args.pairs[i-1], args.pairs[i]]]), axis=0)
    args.pairs = pairs
    
    ##sets the save_name to the file's name if it is not specified
    if args.save_name == None:
        args.save_name = args.file.split(".")[0]
        print "No save name specified, defaulting to %s" % args.save_name
        
    ##Throw an error if args.file is left as none
    if args.file == None or (args.method == "xtc" and args.native == None):
        print "ERROR: NECESSARY FILE NOT SPECIFIED, FAILING..."
    
    return args
        
def get_args():
    ##parent parser for shared parameters
    parser = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    parser.add_argument("--file", type=str, help="File, either .gro or .xtc")
    parser.add_argument("--savedir", type=str, default=os.getcwd(), help="Directory for saving the outputted data")
    parser.add_argument("--pairs", nargs="+",type=int, default=[114,192], help="pairs for computing the y-distance")
    parser.add_argument("--save_name", type=str, help="Name of file to be saved in")
    
    ##The Real Parser
    par = argparse.ArgumentParser(description="Options for finding the pair distance values of a trajectory")
    sub = par.add_subparsers(dest="method")
    
    ##Sub parsers:
    ##gro, for running on a .gro file
    gro_sub = sub.add_parser("gro", parents=[parser], help="Use option 'gro' for .gro files")
    
    ##xtc, for running on a compressed .xtc file
    xtc_sub = sub.add_parser("xtc", parents=[parser], help="Use option 'xtc' for .gro files")
    xtc_sub.add_argument("--native", help="Required native .pdb or conf.gro file")
    
    args = par.parse_args()
    
    args = sanitize_args(args)
    
    return args
    
if __name__ == "__main__":
    args = get_args()
    if args.method == "gro":
        make_gro(args)
    elif args.method == "xtc":
        calc_xtc(args)
    
    
