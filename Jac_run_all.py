"""
This script is for running a jacobian and fitting calculation on every temperature for fitting

"""

import numpy as np
import os
import analysis_scripts.Jac_run_module as jrm
import analysis_scripts.pair_distance_calculator as pdistance

def run_computation(temperature, pairs, spacing, svdt):
    centers_of_bins, normalized_valu, labels, highvalue, lowvalue = jrm.run_main(temperature, pairs, spacing, svdt)
    return centers_of_bins, normalized_valu, labels, highvalue, lowvalue

if __name__ == "__main__":
    spacing = 0.1
    pairs = np.array([[114,192]]) ##for analysis
    temps = np.loadtxt("Temparray.txt")
    cwd = os.getcwd()
    owd = "%s/histograms" % cwd
    
    f = open("fitting.txt","r")
    fsv = open("cutoffs.txt")
    
    method = f.readline()
    svdt = False
    if method[:4] == "tsvd":
        svdt = True
    f.close()
    
    print "svdt is = ", svdt
    if not os.path.isdir(owd):
        os.mkdir(owd)
    print temps
    for t in temps:
        #Analyze and get the new histograms
        print "Starting analysis on temperature %d" % t
        os.chdir("%d"%t)
        centers_of_bins, normalized_valu, labels, highvalue, lowvalue = run_computation(t,pairs, spacing, svdt)
        os.chdir(cwd)
        #Calcualte and Plot the histograms
        os.chdir(owd)
        pdistance.plot_iterations(centers_of_bins, normalized_valu, pairs, labels, spacing, t)
        os.chdir(cwd)
        #write the cutoff into a file for easy access
        fsv.write("The Cutoff occurs between singular values: %.3e and %.3e\n" % (highvalue, lowvalue))
    fsv.close()
    print "Finished all temperatures"
