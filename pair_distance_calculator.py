""" Script for calculating the pair distances for a variety of temperatures

Then prints it out and displays the plot
assumes you are in the directory with T_array.txt

by Justin Chen, December 8, 2014
    
"""
import numpy as np
import mdtraj as md
import os
import matplotlib.pyplot as plt
import matplotlib

def histogram_directory(pairs, spacing):
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
    
    for i in range(np.shape(temp_directory)[0]):
        print "Starting the analysis for directory %s" % temp_directory[i]
        os.chdir(temp_directory[i])
        traj = md.load("traj.xtc", top="Native.pdb")
        compdist = compute_distances(traj, pairs)
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
    

def plot_directory(centers_of_bins, normalized_valu, pairs, temp_directory, spacing):
    #Plot every file, different graphs for different pairs. Outputs a picture to the current directory
    cwd = os.getcwd()
    print "Beginning plotting of directory %s" %cwd
    
    colors = ["b","g","r","c","m","y","b","g","r","c","m","y","b","g","r","c","m","y","b","g","r","c","m","y"]
    linetype = ["-", "--",":"]
    
    fdata = get_FRET_data()
    fhist, fcenter = histogram_data_normalized(fdata, spacing)
    
    for j in range(np.shape(pairs)[0]):
        plt.figure(j)
        maxvalue = 0.0
        maxcenter = 0.0
        plt.plot(fcenter,fhist, alpha=1, color="k", linewidth=2, marker="o", label="D-FRET Data")
        plt.xlabel("R (nm)", fontsize=20)
        plt.ylabel("Probability",fontsize=20)
        plt.title("R-histogram for pair %d-%d"% (pairs[j][0]+1, pairs[j][1]+1), fontsize=20)
        for i in range(np.shape(temp_directory)[0]):
            plt.plot(centers_of_bins[j][i], normalized_valu[j][i], alpha=0.75, linewidth=2, linestyle=linetype[i/6], color=colors[i], label="T-%s"%temp_directory[i][:-2], marker="o")
            maxvalue = find_max(maxvalue, np.max(normalized_valu[j][i]))
            maxcenter = find_max(maxcenter, np.max(centers_of_bins[j][i]))
        plt.axis([0,int(maxcenter*1.5)+1, 0, maxvalue*1.2],fontisze=20) 
        plt.legend()
        plt.savefig("ProbDist-Pair-%d-%d.png"% (pairs[j][0]+1, pairs[j][1]+1))

    print "Finished plotting the directory"

def get_FRET_data():
    found = False
    paths = [ "/home/jchen/projects/2014/10_00_14-FRET/", "/home/jc49/work/", ""]
    for i in paths:
        if os.path.isfile("%sFRET_trace.dat"%i) and (not found):
            found = True 
            fdata = np.loadtxt("%sFRET_trace.dat"%i)
            print "FOUND FRET_data"
    return fdata

def histogram_data_normalized(y, spacing):
    hist, bincenters = histogram_data(y, spacing)
    hist = hist / (spacing*np.sum(hist))
    return hist, bincenters


def histogram_data(y, spacing):
    maxstep = int(np.max(y)/spacing) + 1
    minstep = int(np.min(y)/spacing)
    numbins = maxstep - minstep
    ransize = (minstep*spacing,maxstep*spacing)
    hist, edges = np.histogram(y, numbins, ransize)
    edges = edges + (0.5*spacing)
    bincenters = edges[:-1]
    return hist, bincenters

def compute_distances(traj, prs):
    return md.compute_distances(traj, prs, periodic=False)

def find_max(cmax, nmax):
    if cmax < nmax:
        return nmax
    else:
        return cmax
    

if __name__ == "__main__":
    pairs = np.array([[114,192]])  ##FRET Probes
    spacing = 0.1   # bin spacings
    ran_size = (0,10)    #range of values for the final distance
    #fdata = np.loadtxt("FRET_trace.dat")
    cwd = os.getcwd()
    os.chdir("1PB7")
    os.chdir("iteration_0")
    centers_of_bins, normalized_valu, temp_directory = histogram_directory(pairs, spacing)
    #os.chdir("histanalysis")
    #np.savetxt("FRET_trace.dat", fdata)
    os.chdir("histanalysis")
    plot_directory(centers_of_bins, normalized_valu, pairs, temp_directory, spacing)
    os.chdir(cwd)
    
    
    
    
    
    
    
    
    
    
    
