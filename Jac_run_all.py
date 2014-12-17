"""
This script is for running a jacobian and fitting calculation on every temperature for fitting

"""

import numpy as np
import os
import analysis_scripts.Jac_run_module as jrm

def run_computation(temperature):
    jrm.run_main(temperature)

if __name__ == "__main__":

    temps = np.loadtxt("Temparray.txt")
    cwd = os.getcwd()
    print temps
    for t in temps:
        print "Starting temperature %d" % t
        os.chdir("%d"%t)
        run_computation(t)
        os.chdir(cwd)
    print "Finished all temperatures"
