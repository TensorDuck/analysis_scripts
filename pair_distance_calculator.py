""" Script for calculating the pair distances for a variety of temperatures

Then prints it out and displays the plot
assumes you are in the directory with T_array.txt

by Justin Chen, December 8, 2014
    
"""
import numpy as np
import mdtraj as md
import os
import matplotlib.pyplot as plt
import argparse

def histogram_iterations(pairs,spacing,temperature, fitopts):
    ##assumes you are in the directory with al lthe iterations.
    cwd = os.getcwd()
    print "Beginning histogramming of the directory %s" % cwd
    
    iterations = fitopts["iteration"]+1
    if not os.path.isdir("histanalysis"):
        os.mkdir("histanalysis")
    
    opd = "%s/histanalysis" %cwd
    
    if "y_shift" in fitopts:
	yshift = fitopts["y_shift"]
    else:
        yshift = 0.0
    
    centers_of_bins = [] 
    normalized_valu = []
    labels = []
    
    for i in pairs:## first index, pair, second index, histogram info
        centers_of_bins.append([])
        normalized_valu.append([])
    ##get original data

    #check for all the zeroeth order files
    one_file_found = False
    one_file_notfound = False
    for j in range(np.shape(pairs)[0]):
        if os.path.isfile("T%d_0-pair%d-%d.dat"%(temperature,pairs[j][0]+1, pairs[j][1]+1)):
            one_file_found = True
        else:
            one_file_notfound = True
    #if all zeroeth iteration files found, start the iteration count at 0 for those, and load, otherwise, start others at 0 and continue
    if one_file_found and (not one_file_notfound):
        for j in range(np.shape(pairs)[0]):
            data = np.loadtxt("T%d_0-pair%d-%d.dat"%(temperature,pairs[j][0]+1, pairs[j][1]+1)) 
            centers_of_bins[j].append(data[:,0])
            normalized_valu[j].append(data[:,1])
        labels.append(0)
        start_iteration = 1  
        num_calculated = 1
    else:
        print "Warning, some files missing, skipping and assuming no zeroeth order files exist"
        start_iteration = 0
        num_calculated = 0
        iterations -= 1

    ##start histogram analysis
    for i in range(iterations):
        if (i > iterations-4) or (i%4==0):
            print "Starting histogram analysis"
            labels.append(i+start_iteration)
            num_calculated += 1
            os.chdir("iteration_%d/%d_0"%(i,temperature))
            traj = md.load("traj.xtc", top="Native.pdb")
            compdist = compute_distances(traj, pairs, yshift=yshift)
            for j in range(np.shape(compdist)[1]): ##for each pair
                print "Calculating pair %s" % str(pairs[j]) 
                hist, centers = histogram_data_normalized(compdist[:,j], spacing)
                centers_of_bins[j].append(centers)
                normalized_valu[j].append(hist)
            os.chdir(cwd)
    print "Finished calculating the histograms, begining file writing"
    os.chdir(opd)
    ##start saving the data
    for j in range(np.shape(pairs)[0]):
        for i in range(num_calculated):
            data = np.array([centers_of_bins[j][i],normalized_valu[j][i]])
            data = data.transpose()
            np.savetxt("Iteration%d-pair%d-%d.dat"%(labels[i], pairs[j][0]+1, pairs[j][1]+1), data)
    os.chdir(cwd)  
    print "Completed histogramming and file writing for directory %s" % cwd   
    return centers_of_bins, normalized_valu, labels



def histogram_directory(pairs, spacing, yshift=0.0):
    ##assumes you are in diretory with short_temps.txt
    cwd = os.getcwd()
    print "Beginning histogramming of the directory %s" % cwd
    if not os.path.isfile("short_temps"):
        print "NO short_temps file, need that file to proceed, FAILING"
    
    if not os.path.isdir("histanalysis"):
        os.mkdir("histanalysis")
    
    opd = "%s/histanalysis" %cwd

    temp_directory = np.loadtxt("short_temps", str)
    ##Final matrices. The indices are as follows:
    ##first index = pair
    ##second index = temperature it is taken at
    centers_of_bins = [] 
    normalized_valu = []
    for i in pairs:
        centers_of_bins.append([])
        normalized_valu.append([])
    
    #Begin looping over the temps, then the pairs, to histogram, and then write into files.
    try:
        temp_terms = np.shape(temp_directory)[0]
    except IndexError:
        temp_directory = np.array([temp_directory])
    
    for i in range(np.shape(temp_directory)[0]):
        print "Starting the analysis for directory %s" % temp_directory[i]
        os.chdir(temp_directory[i])
        traj = md.load("traj.xtc", top="Native.pdb")
        compdist = compute_distances(traj, pairs, yshift=yshift)
        for j in range(np.shape(compdist)[1]):
            print "Calculating pair %s" % str(pairs[j])
            hist, centers = histogram_data_normalized(compdist[:,j], spacing)
            centers_of_bins[j].append(centers)
            normalized_valu[j].append(hist)
        os.chdir(cwd)
    print "Finished calculating the histograms, begining file writing"
    
    #write every file
    os.chdir(opd)
    for j in range(np.shape(pairs)[0]):
        for i in range(np.shape(temp_directory)[0]):
            data = np.array([centers_of_bins[j][i],normalized_valu[j][i]])
            data = data.transpose()
            np.savetxt("T%s-pair%d-%d.dat"%(temp_directory[i], pairs[j][0]+1, pairs[j][1]+1), data)
    os.chdir(cwd)  
    
    print "Completed histogramming and file writing for directory %s" % cwd   
    return centers_of_bins, normalized_valu, temp_directory
    
def histogram_multi(args):
    cwd = os.getcwd()
    names = args.names
    spacing = args.spacing
    
    fwd = args.filedir
    swd = args.savedir
    
    fret_data = get_FRET_data(args.fret_data)
    centers_of_bins = [] 
    normalized_valu = []
    for i in args.pairs:
        centers_of_bins.append([])
        normalized_valu.append([])
    
    for name in names:
        if name[-3:] == "dat":
            trace_data = np.loadtxt(name)
            for j in range(np.shape(args.pairs)[0]):
                print "Calculating pair %s" % str(args.pairs[j])
                hist, centers = histogram_data_normalized(trace_data, spacing)
                centers_of_bins[j].append(centers)
                normalized_valu[j].append(hist)
        else:    
            if name[-3:] == "xtc":
                traj_name = md.load("%s/%s"%(fwd, name), top="%s/%s" %(fwd, args.native))
            if name[-3:] == "gro":
                traj_name = md.load("%s/%s"%(fwd, name))
            else:
                raise IOError("Specified file etension not found!")
            compdist = compute_distances(traj_name, args.pairs, yshift=args.y_shift)
            for j in range(np.shape(compdist)[1]):
                print "Calculating pair %s" % str(args.pairs[j])
                hist, centers = histogram_data_normalized(compdist[:,j], spacing)
                centers_of_bins[j].append(centers)
                normalized_valu[j].append(hist)
    
    os.chdir(args.savedir)
    plot_it(centers_of_bins, normalized_valu, args.pairs, args.labels, args.spacing, args.title, fretdata = fret_data)
    os.chdir(cwd)
    

def plot_single(centers_of_bins, normalized_valu, pairs, title, temperature, spacing, fretdata=None):
    plot_it(centers_of_bins, normalized_valu, pairs, "T = %d"%temperature, spacing, title, fretdata=fretdata)   
    
def plot_directory(centers_of_bins, normalized_valu, pairs, temp_directory, spacing):
    for i in np.arange(np.shape(temp_directory)[0]):
        temp_directory[i] = "T=" + temp_directory[i][:-2]
    plot_it(centers_of_bins, normalized_valu, pairs, temp_directory, spacing, "ProbDist")

def plot_iterations(centers_of_bins, normalized_valu, pairs, label, spacing, fit_temp, fretdata=None):
    label_string = []
    for i in label:
        label_string.append("Iter=%d" % i)
        
    plot_it(centers_of_bins, normalized_valu, pairs, label_string, spacing, "T-%d-Iter-%d"%(fit_temp, np.max(label)-1), fretdata=fretdata)
       
    
def plot_it(centers_of_bins, normalized_valu, pairs, label, spacing, title, axis=None, fretdata=None):
    #Plot every file, different graphs for different pairs. Outputs a picture to the current directory
    #centers_of_bins are the bin centers (x), given in [j][i], j=pair, i = label
    #normalized_valu are the normalied histogram values
    #pairs are the list of pairs being plotted
    #label is the legend name for the particular value [i]
    #spacing is the spacing of the histogram, for calculating the FRET data
    #title is the title fo the plot to make
    cwd = os.getcwd()
    print "Beginning plotting of directory %s" %cwd
    
    colors = ["b","g","r","c","m","y","b","g","r","c","m","y","b","g","r","c","m","y","b","g","r","c","m","y"]
    linetype = ["-", "--",":"]
    
    
    if fretdata == None:
        fhist = [0]
        fcenter = [0]
    else:
        fhist, fcenter = histogram_data_normalized(fretdata, spacing)
    for j in range(np.shape(pairs)[0]):
        plt.figure()
        maxvalue = 0.0
        maxcenter = 0.0
        plt.plot(fcenter,fhist, alpha=1, color="k", linewidth=2, marker="o", label="D-FRET Data")
        plt.xlabel("R (nm)", fontsize=20)
        plt.ylabel("Probability",fontsize=20)
        plt.title("R-histogram for pair %d-%d"% (pairs[j][0]+1, pairs[j][1]+1), fontsize=20)
        for i in range(np.shape(label)[0]):
            plt.plot(centers_of_bins[j][i], normalized_valu[j][i], alpha=0.75, linewidth=2, linestyle=linetype[i/6], color=colors[i], label="%s"%label[i], marker="o")
            maxvalue = find_max(maxvalue, np.max(normalized_valu[j][i]))
            maxcenter = find_max(maxcenter, np.max(centers_of_bins[j][i]))
        if axis == None:
            plt.axis([0,int(maxcenter*1.5)+1, 0, maxvalue*1.2],fontsize=20) 
        else:
            plt.axis(axis)
        plt.legend()
        plt.savefig("%s-Pair-%d-%d.png"% (title, pairs[j][0]+1, pairs[j][1]+1))

    print "Finished plotting the directory"

def get_FRET_data(fret_type):
    if fret_type == "obs":
        name = "FRET_trace_obs.dat"
    else:
        name = "FRET_trace.dat"
        
    found = False
    paths = [""]
    fdata = []
    for i in paths:
        if os.path.isfile("%s%s"%(i,name)) and (not found):
            found = True
            print "here is the data"
            print np.loadtxt("%s%s"%(i,name))
             
            fdata = np.loadtxt("%s%s"%(i,name))
            print "FOUND FRET_data"
    
    return fdata

def histogram_data_normalized(y, spacing, wgt=None):
    hist, bincenters = histogram_data(y, spacing, wgt)
    hist = hist / (spacing*np.sum(hist))
    return hist, bincenters


def histogram_data(y, spacing, wgt=None):
    maxstep = int(np.max(y)/spacing) + 1
    minstep = int(np.min(y)/spacing)
    numbins = maxstep - minstep
    ransize = (minstep*spacing,maxstep*spacing)
    hist, edges = np.histogram(y, numbins, ransize, weights=wgt)
    edges = edges + (0.5*spacing)
    bincenters = edges[:-1]
    return hist, bincenters

def compute_distances(traj, prs, yshift=0.0):
    distances = md.compute_distances(traj, prs, periodic=False)
    print distances
    distances += yshift
    print distances
    return distances

def find_max(cmax, nmax):
    if cmax < nmax:
        return nmax
    else:
        return cmax


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
    
    if args.method == "multi":
        if not len(args.labels) == len(args.names):
            if len(args.labels) < len(args.names):
                for i in range(len(args.names)-len(args.labels)):
                     args.labels.append(args.names[len(args.labels)+i][:-4])
                else:
                    print "More labels than names, dropping the ones not corresponding to a file"

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
    parser.add_argument("--fret_data", type=str, default=None, help="specify the type of FRET data using. Either den=Denoised or obs=Observed")
    
    ##The Real Parser
    
    par = argparse.ArgumentParser(description="For deciding which method pair distances will be calculated and made")
    sub = par.add_subparsers(dest="method")
    
    multi_sub = sub.add_parser("multi", parents=[parser], help="for plotting an arbitrary set of files")
    multi_sub.add_argument("--names", type=str, nargs="+", help="Specify names of all the files")
    multi_sub.add_argument("--labels", type=str, nargs="+", help="Specify labels of all the files")
    multi_sub.add_argument("--title", type=str, nargs="+", help="Specify title displayed on the plot")
    multi_sub.add_argument("--native", type=str, help="specify the location of the native.pdb file for use in the filedir")
     
    direc_sub = sub.add_parser("directory", parents=[parser], help="for plotting out a directory from regular MD simulations")
    direc_sub.add_argument("--ran_size", type=float, nargs=2, default=(0,10), help="specify plotting range")
    direc_sub.add_argument("--subdir", type=str, default="1PB7", help="subdir name")
    
    
    args = par.parse_args()

    args = sanitize_args(args)
    
    return args
    

if __name__ == "__main__":
    
    args = get_args()
    
    if args.method == "directory":
        cwd = os.getcwd()
        os.chdir(args.subdir)
        os.chdir("iteration_0")
        centers_of_bins, normalized_valu, temp_directory = histogram_directory(args.pairs, args.spacing, yshift=args.y_shift)
        os.chdir("histanalysis")
        plot_directory(centers_of_bins, normalized_valu, args.pairs, temp_directory, args.spacing)
        os.chdir(cwd)
        
    elif args.method == "multi":
        histogram_multi(args)
    
    
    
    
    
    
    
    
    
    
