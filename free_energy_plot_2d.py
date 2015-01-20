""" Script for computing the free energy plots (2-D) using densplot


"""

import numpy as np
from densplot import hist2d
from math import ceil
import os
import matplotlib.pyplot as plt
import argparse

def combine_iterations_rmsd(start, stop, file_location, output_location, w=False, nbins=50, axisr=None, axisq=None, rmsd_plot=True, Q_plot=True, scatter_only=False, temp=185):
    print "Starting iteration range from %d to %d" % (start, stop)
    #assumes rc1 = rmsd apo, rc2 = rmsd closed
    cwd = os.getcwd()
    rc1n = "rmsd-apo (nm)"
    rc2n = "rmsd-closed (nm)"
    rcqn = "Q-1PB7"
    
    if file_location == None:
        cfd = cwd
    else:
        cfd = file_location
    if output_location == None:
        if not os.path.isdir("plots"):
            os.mkdir("plots")
        csd = "%s/plots" % cwd
    else:
        csd = output_location
    
    #set final matrices, to append all the data onto
    rc1 = np.array([])
    rc2 = np.array([])
    rcq = np.array([])
    wsum = np.sum(np.loadtxt("%s/iter%d.w"%(cfd,start)))
    #add weights matrix if necessary
    if w:
        weights = np.array([])
        for i in np.arange(start, stop+1, 2):
            fw = np.loadtxt("%s/iter%d.w"%(cfd,i))
            weights = np.append(weights, fw, axis=0)
            if np.abs(np.sum(fw)-wsum)>1:
                 print "ERROR, sum of weights changes by more than 1! on iteration %d" % i
    else:
        weights = None
    
    for i in np.arange(start, stop+1, 2):
        f1 = np.loadtxt("%s/iter%d-rmsd-apo.xvg"%(cfd,i), skiprows=13)
        f2 = np.loadtxt("%s/iter%d-rmsd-closed.xvg"%(cfd,i), skiprows=13)
        fq = np.loadtxt("%s/iter%d-Q.out"%(cfd,i))
        rc1 = np.append(rc1, f1[:,1], axis=0)
        rc2 = np.append(rc2, f2[:,1], axis=0)
        rcq = np.append(rcq, fq, axis=0)
    
    os.chdir(csd)
    if Q_plot:
        plot_2D_Free_Energy(rcq, rc2, rcqn, rc2n, "iter%d-%d"%(start,stop), "scatter", weights, nbins, axisq, temp)
        if not scatter_only:
            plot_2D_Free_Energy(rcq, rc2, rcqn, rc2n, "iter%d-%d"%(start,stop), "contour", weights, nbins, axisq, temp)
    if rmsd_plot:
        plot_2D_Free_Energy(rc1, rc2, rc1n, rc2n, "iter%d-%d"%(start,stop), "scatter", weights, nbins, axisr, temp)
        if not scatter_only:
            plot_2D_Free_Energy(rc1, rc2, rc1n, rc2n, "iter%d-%d"%(start,stop), "contour", weights, nbins, axisr, temp)
    
    
    os.chdir(cwd)
    
def plot_2D_Free_Energy(rc1, rc2, rc1n, rc2n,  name, plot_style="scatter", weights=None,  nbins=50, axis=None, temp=300):
    plt.figure()
    label = [rc1n, rc2n]
    
    if weights==None:
        x, y, z = hist2d.make(rc1, rc2, nbins, nbins, temperature=temp, weight=weights, plot_style=plot_style, free_energy_plot=True, idx_smoothing=3)
    else:
        x, y, z = hist2d.make(rc1, rc2, nbins, nbins, temperature=temp, weight=weights, plot_style=plot_style, free_energy_plot=True, idx_smoothing=3)

    if plot_style == 'scatter':
        cp = plt.scatter(x, y, s=10, c=z, marker='o', linewidth=0., vmax=7)
        
    elif plot_style == 'contour':
        cp = plt.contourf(x, y, z, 13)
        plt.contour(x, y, z, cp.levels, colors='k', hold='on')

    if not axis == None:
        plt.axis(axis)
    plt.xlabel(label[0], size=16)
    plt.ylabel(label[1], size=16)

    cb = plt.colorbar(cp)
    pad = 10
    cb.set_label(r'Free Energy (kT units)',labelpad=pad)

    plt.savefig("%s_%s-%s-%s.png"%(name, label[0], label[1], plot_style))
    plt.close()

def check_bool(check):
    if check == "True":
        return True
    else:
        return False
        
if __name__=="__main__":
    ##Parent parser for universal parameters
    par = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    par.add_argument("--file_dir", default="/home/jchen/projects/2015/01-10-1PB7-dmdmd/analysis", type=str)
    par.add_argument("--range", default=[6,100], nargs=2, type=int)
    par.add_argument("--step", type=int)
    par.add_argument("--weight", type=bool, default=True)
    par.add_argument("--bins", type=int, default=50)
    par.add_argument("--scatter_only", action="store_true", default=False)
    par.add_argument("--Q_plot", action="store_false", default=True)
    par.add_argument("--rmsd_plot", action="store_false", default=True)
    
    ##the real parser
    parser = argparse.ArgumentParser(description="For Deciding how to plot the results")
    sub = parser.add_subparsers(dest="method")
    
    ##For defining an axis all the Q and Rs should be plotted on
    same_sub = sub.add_parser("same", parents=[par])
    same_sub.add_argument("--raxis", default=[0,8.5,0,8.5], nargs=4, type=float)
    same_sub.add_argument("--qaxis", default=[0, 1000, 0, 8.5], nargs=4, type=float)
    same_sub.add_argument("--save_dir", default="/home/jchen/analysis/2015/01-10-1PB7-dmdmd/same_axis", type=str)
        
    
    arb_sub = sub.add_parser("arb", parents=[par])
    arb_sub.add_argument("--save_dir", default="/home/jchen/analysis/2015/01-10-1PB7-dmdmd", type=str)

    args = parser.parse_args()
    
    if args.method == "same":
        raxis = args.raxis
        qaxis = args.qaxis
    elif args.method == "arb":
        raxis=None
        qaxis=None
    else:
        print "INVALID METHOD"

    file_location = args.file_dir
    output_location = args.save_dir

    if args.range[0] < 6:
        args.range[0] = 6
        print "Fixing range to lower bound of iteration 6"

    if args.step == None:
        combine_iterations_rmsd(args.range[0], args.range[1], file_location, output_location, w=args.weight, nbins=args.bins, axisr=raxis, axisq=qaxis, rmsd_plot=args.rmsd_plot, Q_plot=args.Q_plot, scatter_only=args.scatter_only )
    else:
        if args.range[0] < args.step:
            for stop in np.arange(args.step, args.range[1], args.step):
                combine_iterations_rmsd(args.range[0], stop, file_location, output_location, w=args.weight, nbins=args.bins, axisr=raxis, axisq=qaxis, rmsd_plot=args.rmsd_plot, Q_plot=args.Q_plot, scatter_only=args.scatter_only)
        else:
            for stop in np.arange(args.range[0]+args.step, args.range[1], args.step):
                combine_iterations_rmsd(args.range[0], stop, file_location, output_location, w=args.weight, nbins=args.bins, axisr=raxis, axisq=qaxis, rmsd_plot=args.rmsd_plot, Q_plot=args.Q_plot, scatter_only=args.scatter_only)
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
