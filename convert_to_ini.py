import model_builder.inputs as inp
import os
import argparse 
import numpy as np

def calc_fit_direc(subdir, iters):

    seq = np.arange(130, 171, 10)
    cwd = os.getcwd()


    for i in seq:
        os.chdir("%d"%i)
        model, fitopt = inp.load_model(subdir)
        fitopt["walltime"] = "16:00:00"
        #fitopt["solver"] = "TSVD"
        #fitopt["T_fit"] = i
        #FRET_pairs = [[115, 193]]
        #fitopt["spacing"] = 0.1
        #fitopt["iteration"] = iters
        inp.save_model(model, fitopt)
        os.chdir(cwd)

def calc_reset_direc(subdir, iters, mparams, pparams):

    seq = np.arange(130, 171, 10)
    cwd = os.getcwd()


    for i in seq:
        os.chdir("%d"%i)
        model, fitopt = inp.load_model(subdir)
        fitopt["iteration"] = iters
        fitopt["pairwise_params_file"] = pparams
        fitopt["model_params_file"] = mparams
        inp.save_model(model, fitopt)
        os.chdir(cwd)

def calc_initial_direc(subdir, iters, t_fit, trunc):

    model, fitopt = inp.load_model(subdir)
    #fitopt["data_type"] = "FRET"
    #fitopt["solver"] = "TSVD"
    fitopt["t_fit"] = t_fit
    #fitopt["fret_pairs"] = [[115, 193]]
    #fitopt["spacing"] = 0.1
    #fitopt["iteration"] = iters
    fitopt["truncate_value"] = trunc
    inp.save_model(model, fitopt)
    
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument("--method", default="initial", type=str)
    parser.add_argument("--subdir", default="1PBQ", type=str)
    parser.add_argument("--iteration", default=0, type=int)
    parser.add_argument("--t_fit", default=130, type=int)
    parser.add_argument("--mparams", default="model_params", type=str)
    parser.add_argument("--pparams", default="pairwise_params", type=str)
    parser.add_argument("--trunc", default=0.01, type=float)
    args = parser.parse_args()
    
    if args.method == "initial":
        calc_initial_direc(args.subdir, args.iteration, args.t_fit, args.trunc)
    elif args.method =="fit":
        calc_fit_direc(args.subdir, args.iteration)
    elif args.method == "reset":
        calc_reset_direc(args.subdir, args.iteration, args.mparams, args.pparams)
