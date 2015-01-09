"""
Script for running a fit for all the solutions at different temperatures
"""

import numpy as np
import os
import analysis_scripts.Jac_save_module as jsm

def run_saving(temperature):
    return jsm.run_main(temperature)
    
if __name__ == "__main__":

    temps = np.loadtxt("Temparray.txt")
    cwd = os.getcwd()
    f = open("Fitting_damping.txt")
    print temps
    for t in temps:
        print "Starting temperature %d" % t
        os.chdir("%d"%t)
        iteration = run_saving(int(t))
        ffit = int(open("1PB7/iteration_%d/newton/fitting_scale").readline())
        if ffit = "1":
            f.write("Temperature %d scaled = True\n"%t)
        else:
            f.write("Temperature %d scaled = False\n"%t)
        os.chdir(cwd)
    f.close()
    print "Finished all temperatures"
