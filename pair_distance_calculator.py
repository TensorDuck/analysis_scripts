""" Script for calculating the pair distances for a variety of temperatures

Then prints it out and displays the plot
assumes you are in the directory with T_array.txt

by Justin Chen, December 8, 2014
    
"""
import numpy as np
import mdtraj as md
import os
import matplotlib.pyplot as plt


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
    labels.append(0)
    for j in range(np.shape(pairs)[0]):
        data = np.loadtxt("T%d_0-pair%d-%d.dat"%(temperature,pairs[j][0]+1, pairs[j][1]+1)) 
        centers_of_bins[j].append(data[:,0])
        normalized_valu[j].append(data[:,1])
    ##start histogram analysis
    num_calculated = 1
    for i in range(iterations):
        if (i > iterations-4) or (i%4==0):
            print "Starting histogram analysis"
            labels.append(i+1)
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
    

def plot_directory(centers_of_bins, normalized_valu, pairs, temp_directory, spacing):
    for i in np.arange(np.shape(temp_directory)[0]):
        temp_directory[i] = "T=" + temp_directory[i][:-2]
    plot_it(centers_of_bins, normalized_valu, pairs, temp_directory, spacing, "ProbDist")

def plot_iterations(centers_of_bins, normalized_valu, pairs, label, spacing, fit_temp):
    label_string = []
    for i in label:
        label_string.append("Iter=%d" % i)
        
    plot_it(centers_of_bins, normalized_valu, pairs, label_string, spacing, "T-%d-Iter-%d"%(fit_temp, np.max(label)-1))
       
    
def plot_it(centers_of_bins, normalized_valu, pairs, label, spacing, title, axis=None):
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
    
    fdata = get_FRET_data()
    if len(fdata) == 0:
        fhist = [0]
        fcenter= [0]
    else:
        fhist, fcenter = histogram_data_normalized(fdata, spacing)

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

def get_FRET_data():
    found = False
    paths = [ "/home/jchen/projects/10_00_14-FRET/", "/home/jc49/work/", ""]
    fdata = []
    for i in paths:
        if os.path.isfile("%sFRET_trace.dat"%i) and (not found):
            found = True
            print "here is the data"
            print np.loadtxt("%sFRET_trace.dat"%i)
             
            fdata = np.loadtxt("%sFRET_trace.dat"%i)
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
    

if __name__ == "__main__":
    pairs = np.array([[114,192]])  ##FRET Probes
    spacing = 0.1   # bin spacings
    ran_size = (0,10)    #range of values for the final distance
    yshift = 0.2
    #fdata = np.loadtxt("FRET_trace.dat")
    cwd = os.getcwd()
    os.chdir("1PB7")
    os.chdir("iteration_0")
    centers_of_bins, normalized_valu, temp_directory = histogram_directory(pairs, spacing, yshift=yshift)
    #os.chdir("histanalysis")
    #np.savetxt("FRET_trace.dat", fdata)
    os.chdir("histanalysis")
    plot_directory(centers_of_bins, normalized_valu, pairs, temp_directory, spacing)
    os.chdir(cwd)
    
    
    
    
    
    
    
    
    
    
    
