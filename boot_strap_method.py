"""
This is a set of methods for performing boot-strapping analysis
"""

import numpy as np
import argparse
import os

import analysis_scripts.plot_package as plt_pkg

def run_analysis(args):
    cwd = os.getcwd()
    savedir = args.savedir
    bootstrap_file = "%s/boot-strap_data.dat"%savedir #
    boostrap_results_file = "%s/boot-strap_results.dat"%savedir #contains the means and standard-deviations for each row
    
    if args.method == "single":
        data_set = np.loadtxt(args.file)
    elif args.method == "multi":
        data_set = get_multi_data_set(args)
    
    data_set = data_set + args.yshift
    #array for plotting
    x = []
    y = []
    yerr = []
    #save the original data in the savedirectory
    np.savetxt("%s/raw_data"%savedir, data_set)
    
    #calculate the parameters for binning
    maxvalue = int(np.amax(data_set)/args.spacing) + 1
    minvalue = int(np.amin(data_set)/args.spacing)
    num_bins = maxvalue - minvalue
    ran_size = (minvalue*args.spacing,maxvalue*args.spacing)
    
    # check for existing boot-strap file
    # if exists, load it. If not, make first entry the new data
    if os.path.isfile(bootstrap_file) and args.append:
        boot_strap = np.loadtxt(bootstrap_file)
    else:
        hist, edges = np.histogram(data_set, bins=num_bins, range=ran_size)
        np.savetxt("%s/original_data"%savedir, hist)
        boot_strap = np.array([hist])
        boot_strap = boot_strap.transpose()
    
    centers = (edges[:-1] + edges[1:])*0.5
    save_centers = np.array([centers]).transpose()
    x.append(centers)
    y.append(hist)
    yerr.append((hist**0.5))
    
    for i in range(args.iters):
        print "starting iteration %d" %i
        new_data = get_new_set(data_set)
        hist, edges = np.histogram(new_data, bins=num_bins, range=ran_size)
        boot_strap = np.append(boot_strap, np.array([hist]).transpose(), axis=1)
        
    np.savetxt(bootstrap_file, boot_strap)
    
    ##calculate the mean distribution, the standard deviation 
    N = np.shape(boot_strap)[1]
    mean = np.sum(boot_strap, axis=1)/N
    y.append(mean) ##append the data now, when it's still a one-D array
    mean = (np.array([mean]).transpose())
    std = np.sum(((boot_strap-mean)**2)/N, axis=1)
    std = std ** 0.5
    yerr.append(std) ##append the data now, when it's still a one-D array
    std = np.array([std]).transpose()
    
    np.savetxt(boostrap_results_file, np.append(np.append(save_centers, mean, axis=1), std, axis=1))
    
    x.append(centers) ##append once more for the x-axis

    
    label=["Original Data", "Boot-strapped"]
    
    plt_pkg.plot_error_bar(x, y, label, "Exp Data", "R (nm)", "Counts", yerr=yerr, axis=args.axis)
    
def run_combine(args):
    x = []
    y = []
    label = []
    yerr = []
    
    for i in args.names:
        data = np.loadtxt("%s/%s.dat" % (args.filedir, i))
        yresult = data[:,1]
        factor = np.sum(yresult)*args.spacing
        x.append(data[:,0])
        y.append(data[:,1]/factor)
        yerr.append(data[:,2]/factor)
        label.append(i)
        
    print np.shape(x[0])
    print np.shape(y[0])
    print np.shape(yerr[0])    
    
    print np.shape(x[1])
    print np.shape(y[1])
    print np.shape(yerr[1])
    plt_pkg.plot_error_bar(x, y, label, "Combined Results", "R (nm)", "Normalized Probability", yerr=yerr)
    
  
def get_multi_data_set(args):
    raise IOError("get_multi_data is not ready for use. Aborting.")
    return [0]
    
    
def get_new_set(list_of_values, weights=None):
    if weights==None:
        return np.random.choice(list_of_values, size=np.shape(list_of_values), replace=True)
    else:
        raise IOError("Not configured to use weights yet, aborting")
        return np.random.choice(list_of_values, size=np.shape(list_of_values), replace=True)

def sanitize_args(args):    
    ##sets the save_name to the file's name if it is not specified
    if args.save_name == None:
        args.save_name = args.file.split(".")[0]
        print "No save name specified, defaulting to %s" % args.save_name
        
    if not os.path.isdir(args.savedir):
        os.mkdir(args.savedir)
    
    return args
        
def get_args():
    ##parent parser for shared parameters
    parser = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    parser.add_argument("--savedir", type=str, default=os.getcwd(), help="Directory for saving the outputted data")
    parser.add_argument("--save_name", default="File", type=str, help="Name of file to be saved in")
    parser.add_argument("--iters", type=int, help="Number of iterations of bootstrapping to perform")
    parser.add_argument("--spacing", type=float, default=0.1, help="the spacing of the data")
    parser.add_argument("--no-append", dest="append", action="store_false", default=True, help="Use this flag if you do not want to append to the previously found bootstrap files") 
    parser.add_argument("--axis", default=None, nargs=4, help="Specify axis dimensions") 
    ##The Real Parser
    par = argparse.ArgumentParser(description="Options for finding the pair distance values of a trajectory")
    sub = par.add_subparsers(dest="method")
    
    ##Sub parsers:
    ##gro, for running on a .gro file
    single_sub = sub.add_parser("single", parents=[parser], help="Use option 'single' for analyzing data from a single file")
    single_sub.add_argument("--file", type=str, help="Specify the file to use located inside the filedir")
    single_sub.add_argument("--yshift", default=0, type=float, help="specify the y-shift to the data")
    
    ##xtc, for running on a compressed .xtc file
    multi_sub = sub.add_parser("multi", parents=[parser], help="Use option 'multi' for analyzing data from multiple files")
    multi_sub.add_argument("--dir_layout", type=str, help="Specify the directory layout for the analysis")
    multi_sub.add_argument("--filedir", type=str, default=os.getcwd(), help="location of files for analysis")
    
    combine_sub = sub.add_parser("combine", parents=[parser], help="Use to combine the results form multiple error analysis results")
    combine_sub.add_argument("--names", nargs="+", type=str)
    combine_sub.add_argument("--filedir", type=str, default=os.getcwd(), help="location of files for analysis")
    
    
    args = par.parse_args()
    
    args = sanitize_args(args)
    
    return args

if __name__ == "__main__":
    args = get_args()
    
    if args.method=="single":
        run_analysis(args)
    elif args.method=="combine":
        run_combine(args)
