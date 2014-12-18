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

def run_main(T_fit, pairs, spacing, svdt):
    
    print pmfit.FRET.compute_Jacobian.def_temp
    pmfit.FRET.compute_Jacobian.def_temp = T_fit
    print pmfit.FRET.compute_Jacobian.def_temp
    cwd = os.getcwd()
    cwd0 = os.getcwd()
    cwd += "/1PB7"

    log = "%s/modelbuilder.log" % cwd

    model = mdb.check_inputs.load_model(cwd, False)
    if svdt:
        model.fitting_solver = "TSVD"
    rcpmanager = MandC2004hack(os.getcwd())
    pmfit.prepare_newtons_method(model,"FRET",rcpmanager.append_log)
    
    os.chdir(cwd)
    centers_of_bins, normalized_valu, labels = pdistance.histogram_iterations(pairs,spacing,T_fit)
    os.chdir(cwd0)
    return centers_of_bins, normalized_valu, labels
    
    
if __name__ == "__main__":
    T_fit = int(np.loadtxt("fitting_temperature.txt"))
    run_main(T_fit, np.array([[114, 192]]), 0.1, True)




