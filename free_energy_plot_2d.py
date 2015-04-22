""" Script for computing the free energy plots (2-D) using densplot


"""

import numpy as np
from densplot import hist2d
from math import ceil
import os
import matplotlib.pyplot as plt
import argparse

def ensure(fname):
    if not os.path.isdir(fname):
        os.mkdir(fname)
    

def get_empty_arrays(weight):
    #return the empty arrays
    #make weights matrix, None if not weighted
    if weight:
        weights = np.array([])
    else:
        weights = None
    return np.array([]), np.array([]), weights

def handle_dmdmd(ext1, ext2, args):
    print "Plotting assuming a dmdmd structure"
    #set the necessary variables from args
    start = args.range[0]
    stop = args.range[1]
    cfd = args.file_dir
    step = args.step
    #set final matrices, to append all the data onto
    rc1, rc2, weights = get_empty_arrays(args.weight)
    rc1n, rc2n = get_labels(ext1, ext2)
    wsum = np.sum(np.loadtxt("%s/iter%d.w"%(cfd,start)))
    
    if args.step == None:
        rc1, rc2, weights = dmdmd_iteration(start, stop, weights, wsum, rc1, rc2, ext1, ext2, cfd)
        plot_2D_Free_Energy(rc1, rc2, rc1n, rc2n, "iter%d-%d"%(start,stop), args, weights=weights, temp=args.temps[0])
    else:
        ##multiple increments, set an array of things to go through and plot
        if start < step:
            fit_range = np.arange(step, stop, step)
        else:
            fit_range = np.arange(start+step, stop, step)
        if not stop==fit_range[-1]:
            fit_range = np.append(fit_range, stop)
        ##Do first step, it's unique
        rc1, rc2, weights = dmdmd_iteration(start, fit_range[0], weights, wsum, rc1, rc2, ext1, ext2, cfd)
        plot_2D_Free_Energy(rc1, rc2, rc1n, rc2n, "iter%d-%d"%(start,fit_range[0]), args, weights=weights, temp=args.temps[0])
        
        ##go through all the possibilities then    
        for i in range(np.shape(fit_range)[0]-1):    
            initial = start
            ##if flow flag is set, this will re=initialize the matrices being merged.
            if args.flow:
                initial = fit_range[i]+2
                rc1, rc2, weights = get_empty_arrays(args.weight)
            rc1, rc2, weights = dmdmd_iteration(fit_range[i]+2, fit_range[i+1], weights, wsum, rc1, rc2, ext1, ext2, cfd)
            plot_2D_Free_Energy(rc1, rc2, rc1n, rc2n, "iter%d-%d"%(initial,fit_range[i+1]), args, weights=weights, temp=args.temps[0])

def dmdmd_iteration(start, stop, weights, wsum, rc1, rc2, ext1, ext2, cfd):    
    ##subroutine handle_dmdmd 
    if not weights == None:
        for i in np.arange(start, stop+1, 2):
            fw = np.loadtxt("%s/iter%d.w"%(cfd,i))
            weights = np.append(weights, fw, axis=0)
            if np.abs(np.sum(fw)-wsum)>1:
                 print "ERROR, sum of weights changes by more than 1! on iteration %d" % i
    else:
        weights = None
    #loop through and combine all the data from every iteration.
    for i in np.arange(start, stop+1, 2):
        f1 = get_value("iter%d"%i, ext1, cfd)
        f2 = get_value("iter%d"%i, ext2, cfd)
        rc1 = np.append(rc1, f1, axis=0)
        rc2 = np.append(rc2, f2, axis=0)
    
    return rc1, rc2, weights


def handle_dmaps(ext1, ext2, args):
    #set the necessary variables from args
    start = args.range[0]
    stop = args.range[1]
    cfd = args.file_dir
    step = args.step
    #set final matrices, to append all the data onto
    rc1, rc2, weights = get_empty_arrays(args.weight)
    rc1n, rc2n = get_labels(ext1, ext2)
    
    
    for i in np.arange(start, stop+step, step):
        rc1 = get_value("iter%d"%i, ext1, cfd)
        rc2 = get_value("iter%d"%i, ext2, cfd)
        weights = np.loadtxt("%s/iter%d.w"%(cfd,i))
        plot_2D_Free_Energy(rc1, rc2, rc1n, rc2n, "iter%d"%i, args, weights=weights, temp=args.temps[0])
   
def handle_fret(ext1, ext2, args):
    temperatures = args.temps
    iteration_range = args.range
    cfd = args.file_dir
    rc1n, rc2n = get_labels(ext1, ext2)
    for T in temperatures:
        for j in np.arange(iteration_range[0], iteration_range[1]+1, 1):
            rc1 = get_value("T%dI%d"%(T,j), ext1, cfd)
            rc2 = get_value("T%dI%d"%(T,j), ext2, cfd)
            plot_2D_Free_Energy(rc1, rc2, rc1n, rc2n, "Temp-%d-Iter-%d"%(T, j+1), args, temp=T)
    
    
def handle_vanilla(ext1, ext2, args):
    temperature = args.temps[0]
    iteration_range = args.range
    cfd = args.file_dir
    
    rc1 = get_value(args.file_names[0], ext1, cfd)
    rc2 = get_value(args.file_names[1], ext2, cfd)
    
    rc1n, rc2n = get_labels(ext1, ext2)
    
    plot_2D_Free_Energy(rc1, rc2, rc1n, rc2n,  args.save_name, args, temp=temperature)
    

def get_labels(ext1, ext2):
    labels = {"-Qclosed.out":"Q closed", "-y114-192.out":"y between 115-193 (nm)", "-rmsd-closed.xvg":"rmsd-closed (nm)", "-rmsd-apo.xvg":"rmsd-apo (nm)", "-comA.xvg":"Lobe Center Distances (nm)", "ev1":"DC1", "ev2":"DC2"}
    return labels[ext1], labels[ext2]
   
def plot_2D_Free_Energy(rc1, rc2, rc1n, rc2n, name, args, weights=None, temp=300):
    print "Starting the Plot"
    nbins = args.bins
    if args.method == "arb":
        axis = None
    elif args.method == "same":
        axis = args.axis
    ##Plot a scatter plot
    if not args.contour_only:
        plt.figure()
        plot_style = "scatter"
        x, y, z = hist2d.make(rc1, rc2, nbins, nbins, temperature=temp, weight=weights, plot_style=plot_style, free_energy_plot=True, idx_smoothing=3)
        cp = plt.scatter(x, y, s=10, c=z, marker='o', linewidth=0.)
        #vmax=7
        if not axis == None:
            plt.axis(axis)
        plt.xlabel(rc1n, size=16)
        plt.ylabel(rc2n, size=16)

        cb = plt.colorbar(cp)
        pad = 10
        cb.set_label(r'Free Energy (kT units)',labelpad=pad)

        plt.savefig("%s/%s_%s-%s-%s.png"%(args.save_dir, name, rc1n, rc2n, plot_style))
        plt.close()
    
    
    ##Plot a contour plot   
    if not args.scatter_only:
        ensure("%s/contour"%args.save_dir)
        plt.figure()
        plot_style = "contour"
        xc, yc, zc = hist2d.make(rc1, rc2, nbins, nbins, temperature=temp, weight=weights, plot_style=plot_style, free_energy_plot=True, idx_smoothing=3)
        cp = plt.contourf(xc, yc, zc, 20)
        plt.contour(xc, yc, zc, cp.levels, colors='k', hold='on')
        if not axis == None:
            plt.axis(axis)
        plt.xlabel(rc1n, size=16)
        plt.ylabel(rc2n, size=16)

        cb = plt.colorbar(cp)
        pad = 10
        cb.set_label(r'Free Energy (kT units)',labelpad=pad)

        plt.savefig("%s/contour/%s_%s-%s-%s.png"%(args.save_dir, name, rc1n, rc2n, plot_style))
        plt.close()

def get_value(name, ext, cfd):
    if ext == "-comA.xvg":
        return np.loadtxt("%s/%s%s"%(cfd, name, ext), skiprows=22)[:,1]   
    elif ext[-4:] == ".xvg":
        return np.loadtxt("%s/%s%s"%(cfd,name, ext), skiprows=13)[:,1]   
    elif ext[-4:] == ".out":
        return np.loadtxt("%s/%s%s"%(cfd,name, ext))
    elif ext == "ev1":
        values = np.loadtxt("%s/%s.ev"%(cfd,name))
        return values[:,0]
    elif ext == "ev2":
        values = np.loadtxt("%s/%s.ev"%(cfd,name))
        return values[:,1] 
    else:
        print "ERROR: COULD NOT DETERMINE FILE TYPE, RETURNING NONE"
        return None
    
def sanitize_args(args):
    original_directory = os.getcwd()
    os.chdir(args.file_dir)
    ##If temps was not specified, will attempt to open Temparray.txt, Otherwise prints ERROR
    if not args.temps == None:
        pass
    elif os.path.isfile("Temparray.txt"):
        args.temps = np.loadtxt("Temparray.txt",dtype=int)
    else:
        print "ERROR: No Temperature Directories Specified"
        print "Defaulting to 185K"
        args.temps=[185.0]
    
    ##Check handles, make sure it's specified
    if args.handle == None:
        print "ERROR: No DIRECTORY STUCTURE SPECIFIED, NO DEFAULT AVAILABLE"
    elif args.handle == "dmdmd":
        if len(args.temps) > 1:
            print "Too many temperatures specified, defaulting to first temperature of %d" % args.temps[0]
    elif args.handle =="fret":
        pass
    elif args.handle == "vanilla":
        pass
    elif args.handle == "dmaps":
        pass
    else:
        print "ERROR: NO DIRECTORY STUCTURE SPECIFIED, NO DEFAULT AVAILABLE"
    
    ##check if axis specified. If not, specify the defaults per axeses pieces
    axeses = {"Q":[0,1000], "A":[0, 8.5], "C":[0, 8.5], "Y":[0, 50]}
    if args.method == "same" and args.axis == None:
        if (args.plot_type[0] not in axeses) and (args.plot_type[1] not in axeses):
            print "ERROR: INVALID PLOT TYPE SPECIFIED"
        else:
            a1 = axeses[args.plot_type[0]]
            a2 = axeses[args.plot_type[1]]
            args.axis = [a1[0], a1[1], a2[0], a2[1]]
    
    os.chdir(original_directory)
    
    ##make directories if they do not exist
    if not os.path.isdir(args.save_dir):
        os.mkdir(args.save_dir)
    
    ##make names at least a length 2 list
    if len(args.file_names)==0:
        if args.handle=="vanilla":
            raise IOError("Must specify file name")
    elif len(args.file_names)==1:
        args.file_names.append(args.names[0])
        
        
    return args
    
    
def get_args():
    ##Parent parser for universal parameters
    par = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    par.add_argument("--file_dir", default="%s/analysis"%os.getcwd(), type=str, help="location of files for analysis")
    par.add_argument("--range", default=[6,100], nargs=2, type=int, help="Range of iterations to plot for")
    par.add_argument("--step", type=int, help="Use for specifying what step-size to take in iterations for plotting")
    par.add_argument("--no_weight", dest="weight", action="store_false", default=True, help="Use if frames don't have a weight in a .w file")
    par.add_argument("--bins", type=int, default=50, help="Number of bins in each axis for binning data")
    par.add_argument("--scatter_only", action="store_true", default=False, help="use for plotting only a scatter plot")
    par.add_argument("--contour_only", action="store_true", default=False, help="use for plotting only a contour plot")
    par.add_argument("--plot_type", type=str, default="QC", help="specify the type of plot in xy format; default QC. C=RMSD-closed, A=RMSD-apo, Q=Q, Y=FRET probe distance, L=Lobes center distance, 1=DC1, 2=DC2")
    par.add_argument("--handle", type=str, help="specify either dmdmd, vanilla or fret")
    par.add_argument("--temps", type=float, nargs="+", help="specify the temperature for the data, can be an array")
    par.add_argument("--flow", action="store_true", default=False, help="Use if you want to plot iterations in intervals, i.e. 2-50, 52-60")
    ##for just handle: vanilla
    par.add_argument("--file_names", nargs="+", type=str, help="Specify the specific names of the DC containing files") ##for handle vanilla
    par.add_argument("--save_name", type=str, default="fep", help="Specify the specific names of the DC containing files") ##for handle vanilla
    ##the real parser
    parser = argparse.ArgumentParser(description="For Deciding how to plot the results")
    sub = parser.add_subparsers(dest="method")
    
    ##For defining an axis all the Q and Rs should be plotted on
    same_sub = sub.add_parser("same", parents=[par], help="for plotting all the plots on the same axis specified by 'raxis', 'qaxis'")
    same_sub.add_argument("--axis", nargs=4, type=float, help="specify axis. Defaults set for supported plot_type")
    same_sub.add_argument("--save_dir", default="%s/same_axis"%os.getcwd(), type=str, help="directory for saving the plots")
        
    
    arb_sub = sub.add_parser("arb", parents=[par], help="Plot all plots on arbitrary axis (python picks)")
    arb_sub.add_argument("--save_dir", default=os.getcwd(), type=str, help="directory for saving the plots")

    args = parser.parse_args()
    args = sanitize_args(args)
    return args
   
if __name__=="__main__":
    
    args = get_args()
    
    #if args.dir_structure ==  
    #keys of different methods, asign it to the handle which is then called to do the rest
    handlers = {"dmdmd":handle_dmdmd, "fret":handle_fret, "vanilla":handle_vanilla, "dmaps":handle_dmaps}  
    names = {"Q":"-Qclosed.out", "A":"-rmsd-apo.xvg", "C":"-rmsd-closed.xvg", "Y":"-y114-192.out", "L":"-comA.xvg", "1":"ev1", "2":"ev2"}
    
    handle = handlers[args.handle]
    rcn1 = names[args.plot_type[0]]
    rcn2 = names[args.plot_type[1]]
    
    handle(rcn1, rcn2, args)

 
    '''
    
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
        
    
    
    
    
    
    '''
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
