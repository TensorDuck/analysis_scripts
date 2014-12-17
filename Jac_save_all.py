"""
Script for running a fit for all the solutions at different temperatures
"""

import numpy as np
import os
import analysis_scripts.Jac_save_module as jsm

def run_saving(temperature):
    jsm.run_main(temperature)

if __name__ == "__main__":

    temps = np.loadtxt("Temparray.txt")
    cwd = os.getcwd()
    print temps
    for t in temps:
        print "Starting temperature %d" % t
        os.chdir("%d"%t)
        run_saving(int(t))
        os.chdir(cwd)
    print "Finished all temperatures"
