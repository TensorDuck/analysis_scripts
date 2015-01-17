"""
This is a script for testing the functionality of the compute_Jacobian module for FRET

"""

import numpy as np
import model_builder as mdb
import os as os
import project_tools.parameter_fitting as pmfit
from FRET_experiment.MandC2004_hacked import MandC2004hack
import matplotlib.pyplot as plt
import project_tools.parameter_fitting.FRET.truncated_SVD_FRET as tsvd
import analysis_scripts.pair_distance_calculator as pdistance

def run_main(T_fit, pairs, spacing, svdt, subfolder="1PB7"):
    
    pmfit.FRET.compute_Jacobian.def_temp = T_fit
    cwd = os.getcwd()
    cwd0 = os.getcwd()
    cwd += "/%s" % subfolder
    
    log = "%s/modelbuilder.log" % cwd

    model = mdb.check_inputs.load_model(cwd, False)
    if svdt:
        model.fitting_solver = "TSVD"
    rcpmanager = MandC2004hack(os.getcwd())
    pmfit.prepare_newtons_method(model,"FRET",rcpmanager.append_log)
    '''
    #This will make a histogram of all the different iterations and save the data accordingly
    os.chdir(cwd)
    centers_of_bins, normalized_valu, labels = pdistance.histogram_iterations(pairs,spacing,T_fit)
    os.chdir(cwd0)
    '''
    ##Following will estimate the expected cutoff, along with return the boundaries of the cutoff
    newtondir = "%s/iteration_%d/newton" % (cwd,model.iteration)
    os.chdir(newtondir)
    highvalue, lowvalue = estimate_lambda()
    os.chdir(cwd0)    
    
    return centers_of_bins, normalized_valu, labels, highvalue, lowvalue

def estimate_lambda():
    print "estimating the value of lambda from singular values"    
    svf = np.loadtxt("singular_values.dat") 
    index = 0
    num = np.shape(svf)[0]
    lowvalue = np.min(svf)
    highvalue = np.min(svf)
    for i in range(num-1):
        if svf[i]/svf[i+1] > 1000:
            index = num - 1 - i
            highvalue = svf[i]
            lowvalue = svf[i+1]   
    open("Lambda_index.txt","w").write("%d"%index)
    return highvalue, lowvalue

if __name__ == "__main__":
    T_fit = int(np.loadtxt("fitting_temperature.txt"))
    run_main(T_fit, np.array([[114, 192]]), 0.1, True, "1PBQ")




