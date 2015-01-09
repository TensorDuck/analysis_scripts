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
    fit_string = ""
    print temps
    for t in temps:
        print "Starting temperature %d" % t
        os.chdir("%d"%t)
        iteration = run_saving(int(t))
        ffit = int(open("1PB7/iteration_%d/newton/fitting_scale"%iteration).readline())
        if ffit == "1":
            fit_string += "Temperature %d scaled = True\n"%t
        else:
            fit_string += "Temperature %d scaled = False\n"%t
        os.chdir(cwd)
    
    f = open("Fitting_damping_%d.txt"%iteration,"w")
    f.write("fit_string")
    f.close()
    print "Finished all temperatures"
