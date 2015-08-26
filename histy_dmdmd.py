"""
This is a script for histogramming the FRET y-distance for a dmdmd-set of files
"""

import numpy as np
import matplotlib.pyplot as plt
import os as os
import argparse
import mdtraj as md

import analysis_scripts.gro_reader as gread
import analysis_scripts.pair_distance_calculator as pdistance
import analysis_scripts.free_energy_plot_1d as fep1

from analysis_scripts.gro_reader import traj

kb = 6.022141*1.380650/(4.184*1000.0)

def handle_dmdmd(args):
    cwd = os.getcwd()
    
    #Load every array
    os.chdir(args.filedir)
    
    all_y, all_Q, all_w = get_qy(args)
    
    all_y = all_y + args.y_shift
    os.chdir(cwd)
    
    hist_y(all_y, all_Q, all_w, args)
    

def handle_dmaps(args):
    cwd = os.getcwd()
    
    #Load every array
    os.chdir(args.filedir)
    
    all_y = np.array([]) # FRET Probe C-alpha distance
    all_Q = np.array([]) # corresponding Q values
    all_w = np.array([]) # corresponding weights (.w) to the files
    
    for i in np.arange(args.iters[0], args.iters[1]+1, args.iters[2]):
        y = np.loadtxt("iter%d-y114-192.out" % i)
        Q = np.loadtxt("iter%d-Qclosed.out" % i)
        w = np.loadtxt("iter%d.w" % i)
        all_y = np.append(all_y, y, axis=0)
        all_Q = np.append(all_Q, Q, axis=0)
        all_w = np.append(all_w, w, axis=0)
    
    os.chdir(cwd)
    all_y = all_y + args.y_shift
    
    hist_y(all_y, all_Q, all_w, args)

def hist_y(all_y, all_Q, all_w, args):    
    cwd = os.getcwd()
    # Plot Q
    for Q_bound in args.QQ:
        q_larger = all_Q >= Q_bound

        plot_y = all_y[q_larger]    
        plot_w = all_w[q_larger]
        plot_Q = all_Q[q_larger] 
        
        ran_size = (np.floor(np.min(plot_y)/args.spacing)*args.spacing,(np.floor(np.max(plot_y)/args.spacing)*args.spacing) + args.spacing )
    
        nbins = (ran_size[1]-ran_size[0])/args.spacing
    
        hist, edges = np.histogram(plot_y, bins=nbins, range=ran_size, normed=True, density=True)

        edges = edges + (0.5*args.spacing)
        bincenters = edges[:-1]
        
        normalized_value = [[hist]]
        centers_of_bins = [[bincenters]]
        label = ["Q>%d" % Q_bound]
        title = "%s-Q-%d" % (args.save_name, Q_bound)
        
        print "Shape of the data is:"
        print np.shape(normalized_value)
        os.chdir(args.savedir)
        pdistance.plot_it(centers_of_bins, normalized_value, args.pairs, label, args.spacing, title, axis=args.axis)
        np.savetxt("%s-data.dat"%title, np.array([bincenters, hist]).transpose())
        os.chdir(cwd)
            
    # plot Q free Energy
    
    ran_sizeQ = (np.floor(np.min(all_Q)/args.Qspacing)*args.Qspacing, (np.floor(np.max(all_Q)/args.Qspacing)+1)*args.Qspacing)

    nbinsQ = (ran_sizeQ[1]-ran_sizeQ[0])/args.Qspacing
    histQ, edgesQ = np.histogram(all_Q, bins=nbinsQ, range=ran_sizeQ, density=True)
    Qcenters = edgesQ+(0.5*args.Qspacing)
    Qcenters = Qcenters[:-1]
    energy = np.log(histQ) * (-1*kb*args.temperature)
    
    energies= [energy]
    coord_centers = [Qcenters]
    titleQ = "%s_fep" % args.save_name
    labelQ = ["T=%d"%args.temperature]
    os.chdir(args.savedir)
    fep1.plot_Q(energies, coord_centers, labelQ, titleQ, axis=[0, 1000, 0, 6])
    os.chdir(cwd)
               
def get_qy(args):
    all_y = np.array([]) # FRET Probe C-alpha distance
    all_Q = np.array([]) # corresponding Q values
    all_w = np.array([]) # corresponding weights (.w) to the files

    for i in np.arange(args.iters[0], args.iters[1]+1, args.iters[2]):
        y = np.loadtxt("iter%d-y114-192.out" % i)
        Q = np.loadtxt("iter%d-Qclosed.out" % i)
        w = np.loadtxt("iter%d.w" % i)
        all_y = np.append(all_y, y, axis=0)
        all_Q = np.append(all_Q, Q, axis=0)
        all_w = np.append(all_w, w, axis=0)
    
    all_y = all_y + args.y_shift
    return all_y, all_Q, all_w


        
def sanitize_args(args):
    ##set the pairs for fitting in the array format for the calculation
    pairs = np.array([[args.pairs[0], args.pairs[1]]])
    #print "Number of pairs is: %d" % (len(args.pairs)/2)
    if len(args.pairs)>2:
        for i in np.arange(3, len(args.pairs), 2):
            pairs = np.append(pairs, np.array([[args.pairs[i-1], args.pairs[i]]]), axis=0)
    args.pairs = pairs
    
    ##sets the save_name to the file's name if it is not specified
    if args.save_name == None:
        args.save_name = "hist%d-%d"%(args.iters[0],args.iters[1])
        print "No save name specified, defaulting to %s" % args.save_name
    
    if not os.path.isdir(args.savedir):
        os.mkdir(args.savedir)
    '''
    if type(args.ran_size) == list:
        args.ran_size=(args.ran_size[0], args.ran_size[1])
       ''' 
    return args
        
def get_args():
    ##parent parser for shared parameters
    parser = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    parser.add_argument("--filedir", type=str, default=os.getcwd(), help="File, either .gro or .xtc")
    parser.add_argument("--savedir", type=str, default="%s/hist_plots"%os.getcwd(), help="Directory for saving the outputted data")
    parser.add_argument("--pairs", nargs="+",type=int, default=[114,192], help="pairs for computing the y-distance")
    parser.add_argument("--save_name", type=str, default="plot", help="Name of file to be saved in")
    parser.add_argument("--spacing", type=float, default=0.1, help="spacing in nm for histogramming data")
    #parser.add_argument("--ran_size", type=float, nargs=2, default=(0,10), help="range of y values to plot")
    parser.add_argument("--iters", type=int, nargs=3, default=[6, 8, 2], help="range of iter files to consider")
    parser.add_argument("--QQ", type=int, nargs="+", default=[800, 900], help="range of Q values to cutoff and use")
    parser.add_argument("--Qspacing", type=int, default=10, help="bin size for Q")
    parser.add_argument("--temperature", type=float, default=170, help="temperature of simulation")
    parser.add_argument("--y_shift", type=float, default=0, help="Specify the y-shift to the FRET distance data")
    parser.add_argument("--handle", type=str, default="dmdmd", help="Specify the directory structure")
    parser.add_argument("--axis", type=float, default=None, nargs=4, help="Specify the axis in 'x0 x1 y0 y1' format")
    ##The Real Parser
    par = argparse.ArgumentParser(description="Options for finding the pair distance values of a trajectory", parents=[parser])
    
    
    args = par.parse_args()

    args = sanitize_args(args)
    
    return args

if __name__ == "__main__":
    args = get_args()
    print args.save_name
    print args.filedir
    print args.savedir
    print args.pairs
    print args.spacing
    #print args.ran_size
    print args.iters
    print args.QQ
    
    handlers = {"dmdmd":handle_dmdmd, "dmaps":handle_dmaps}
    
    handle = handlers[args.handle]
    handle(args)
    
        
    
    
