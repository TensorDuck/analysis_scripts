"""
This script is for running a jacobian and fitting calculation on every temperature for fitting

"""

import numpy as np
import os
import analysis_scripts.Jac_run_module as jrm
import analysis_scripts.pair_distance_calculator as pdistance

def run_computation(temperature, pairs, spacing, svdt):
    centers_of_bins, normalized_valu, labels = jrm.run_main(temperature, pairs, spacing, svdt)
    return centers_of_bins, normalized_valu, labels

if __name__ == "__main__":
    spacing = 0.1
    pairs = np.array([[114,192]]) ##for analysis
    temps = np.loadtxt("Temparray.txt")
    cwd = os.getcwd()
    owd = "%s/histograms" % cwd
    
    f = open("fitting.txt","r")
    method = f.readline()
    svdt = False
    if method == "tsvd":
        svdt = True
    f.close()
    
    if not os.path.isdir(owd):
        os.mkdir(owd)
    print temps
    for t in temps:
        print "Starting temperature %d" % t
        os.chdir("%d"%t)
        centers_of_bins, normalized_valu, labels = run_computation(t,pairs, spacing, svdt)
        os.chdir(cwd)
        os.chdir(owd)
        pdistance.plot_iterations(centers_of_bins, normalized_valu, pairs, labels, spacing, t)
        os.chdir(cwd)
        
    print "Finished all temperatures"
