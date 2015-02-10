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
        fitopt["data_type"] = "FRET"
        fitopt["solver"] = "TSVD"
        fitopt["T_fit"] = i
        FRET_pairs = [[115, 193]]
        fitopt["spacing"] = 0.1
        fitopt["iteration"] = iters
        inp.save_model(model, fitopt)
        os.chdir(cwd)

def calc_initial_direc(subdir, iters):

    model, fitopt = inp.load_model(subdir)
    fitopt["data_type"] = "FRET"
    fitopt["solver"] = "TSVD"
    fitopt["T_fit"] = 130
    fitopt["FRET_pairs"] = [[115, 193]]
    fitopt["spacing"] = 0.1
    fitopt["iteration"] = iters
    inp.save_model(model, fitopt)
    
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument("--method", default="initial", type=str)
    parser.add_argument("--subdir", default="1PBQ", type=str)
    parser.add_argument("--iteration", default=1, type=int)
    args = parser.parse_args()
    
    if args.method == "initial":
        calc_initial_direc(args.subdir, args.iteration)
    elif args.method =="fit":
        calc_fit_direc(args.subdir, args.iteration)
