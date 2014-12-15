"""
Script for running a fit for all the solutions at different temperatures
"""

import numpy as np
import os


if __name__ == "__main__":

    temps = np.loadtxt("Temparray.txt")
    cwd = os.getcwd()
    for t in temps:
        print "Starting temperature %d" % t
        os.chdir(t)
        import analysis_scripts.Jac_save_module as jsm
        jsm.run_main()
        os.chdir(cwd)
    print "Finished all temperatures"
