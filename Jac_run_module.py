"""
This is a script for testing the functionality of the compute_Jacobian module for FRET

"""

import numpy as np
import model_builder as mdb
import os as os
import project_tools.parameter_fitting as pmfit
from FRET_experiment.MandC2004_hacked import MandC2004hack
import matplotlib.pyplot as plt

def run_main(T_fit):
    print pmfit.FRET.compute_Jacobian.def_temp
    pmfit.FRET.compute_Jacobian.def_temp = T_fit
    print pmfit.FRET.compute_Jacobian.def_temp
    cwd = os.getcwd()
    cwd0 = os.getcwd()
    cwd += "/1PB7"

    log = "%s/modelbuilder.log" % cwd

    model = mdb.check_inputs.load_model(cwd, False)

    rcpmanager = MandC2004hack(os.getcwd())
    pmfit.prepare_newtons_method(model,"FRET",rcpmanager.append_log)

if __name__ == "__main__":
    T_fit = int(np.loadtxt("fitting_temperature.txt"))
    run_main(T_fit)




