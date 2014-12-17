"""
This script will calculate the new epsilon parametesr for a temperature T_fit

then starts a new simulation at that one temperature.

Meant for the FRET-fitting package in Paramter_fitting
this is  

"""

import numpy as np
#import project_tools.parameter_fitting.FRET.compute_Jacobian as compJ
import model_builder as mdb
import os as os
import project_tools.parameter_fitting as pmfit
from FRET_experiment.MandC2004_hacked import MandC2004hack
from project_tools import simulation
import matplotlib.pyplot as plt

def run_main(T_fit):
    cwd = os.getcwd()
    cwd0 = os.getcwd()
    cwd += "/1PB7"
    log = "%s/modelbuilder.log" % cwd
    model = mdb.check_inputs.load_model(cwd, False)

    rcpmanager = MandC2004hack(os.getcwd())
    append_log = rcpmanager.append_log
    #pmfit.prepare_newtons_method(model,"FRET",rcpmanager.append_log)
    pmfit.save_new_parameters(model,"FRET",append_log)
    
    
    model.iteration += 1
    newdirec = "%s/iteration_%d" % (cwd, model.iteration)
    
    if not os.path.isdir(newdirec):
        os.mkdir(newdirec)
    
    os.chdir(newdirec)
    
    simulation.constant_temp.run_temperature_array(model,T_fit,T_fit,5)
    
    append_log(model.subdir,"Submitting short_temps iteration %d " % model.iteration)
    append_log(model.subdir,"  T_min = %d , T_max = %d , dT = %d" % (T_fit, T_fit, 5))
    append_log(model.subdir,"Starting: Tf_loop_iteration")
    
    os.chdir(cwd0)
    open(model.subdir+"/model.info","w").write(model.get_model_info_string())
    
if __name__ == "__main__":
    T_fit = int(np.loadtxt("fitting_temperature.txt"))
    run_main(T_fit)
    
    print "GOT TO END"



