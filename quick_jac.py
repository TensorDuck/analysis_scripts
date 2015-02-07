import numpy as np
import model_builder as mdb
import os as os
import project_tools.parameter_fitting as pmfit
import argparse
from analysis_scripts.recipe_log_function import log_function
from analysis_scripts.JacRunModule import estimate_lambda

def run_start(args):
    T_fit = args.T_fit
    print pmfit.FRET.compute_Jacobian.def_temp
    pmfit.FRET.compute_Jacobian.def_temp = T_fit
    print pmfit.FRET.compute_Jacobian.def_temp
    cwd = os.getcwd()
    cwd0 = os.getcwd()
    cwd += "/%s" % args.subdir

    log = "%s/modelbuilder.log" % cwd

    model = mdb.inputs.load_model(cwd, False)
    model.fitting_solver = "TSVD"
    rcpmanager = log_function(os.getcwd())
    pmfit.prepare_newtons_method(model,"FRET",rcpmanager.append_log)
    
def run_end(args):
    cwd = os.getcwd()
    cwd0 = os.getcwd()
    cwd += "/%s" % args.subdir
    log = "%s/modelbuilder.log" % cwd
    model = mdb.inputs.load_model(cwd, False)

    rcpmanager = log_function(os.getcwd())
    append_log = rcpmanager.append_log
    #pmfit.prepare_newtons_method(model,"FRET",rcpmanager.append_log)
    pmfit.save_new_parameters(model,"FRET",append_log)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    parser.add_argument("--subdir", default="1PB7", type=str, help="Sub Directory the model.info file is located in")
    parser.add_argument("--T_fit", default=130, type=int)
    
    cwd = os.getcwd()
    
    args = parser.parse_args()
   
    run_start(args)
    
    os.chdir("%s/%s/iteration_0/newton" % (cwd, args.subdir))
    estimate_lambda()
    
    os.chdir(cwd)
    run_end(args)
    
    
    

