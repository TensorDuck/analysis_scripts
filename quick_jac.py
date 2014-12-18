import numpy as np
import model_builder as mdb
import os as os
import project_tools.parameter_fitting as pmfit
from FRET_experiment.MandC2004_hacked import MandC2004hack

def run_start():
    T_fit = 160
    print pmfit.FRET.compute_Jacobian.def_temp
    pmfit.FRET.compute_Jacobian.def_temp = T_fit
    print pmfit.FRET.compute_Jacobian.def_temp
    cwd = os.getcwd()
    cwd0 = os.getcwd()
    cwd += "/1PB7"

    log = "%s/modelbuilder.log" % cwd

    model = mdb.check_inputs.load_model(cwd, False)
    model.fitting_solver = "TSVD"
    rcpmanager = MandC2004hack(os.getcwd())
    pmfit.prepare_newtons_method(model,"FRET",rcpmanager.append_log)
    
def run_end():
    cwd = os.getcwd()
    cwd0 = os.getcwd()
    cwd += "/1PB7"
    log = "%s/modelbuilder.log" % cwd
    model = mdb.check_inputs.load_model(cwd, False)

    rcpmanager = MandC2004hack(os.getcwd())
    append_log = rcpmanager.append_log
    #pmfit.prepare_newtons_method(model,"FRET",rcpmanager.append_log)
    pmfit.save_new_parameters(model,"FRET",append_log)

if __name__ == "__main__":
    #run_start()
    run_end()

