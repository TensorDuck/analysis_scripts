"""
Script for running a fit for all the solutions at different temperatures
"""

import numpy as np
import os

def run_saving():
    import analysis_scripts.Jac_save_module as jsm
    jsm.run_main()

if __name__ == "__main__":

    temps = np.loadtxt("Temparray.txt")
    cwd = os.getcwd()
    print temps
    for t in temps:
        print "Starting temperature %d" % t
        os.chdir(str(t))
        run_saving()
        os.chdir(cwd)
    print "Finished all temperatures"
