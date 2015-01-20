"""
This is a script for calculating the jacobin and FRET bins of a protein

It performs all the steps of: calculating Jacobian and new parameters, saving new parameters for the next model, and submitting a job

"""

import numpy as np
import matplotlib.pyplot as plt
import os as os
import argparse

import model_builder as mdb
import project_tools.parameter_fitting as pmfit
from project_tools import simulation

import analysis_scripts.pair_distance_calculator as pdistance
from analysis_scripts.recipe_log_function import log_function

def run_calc_all(args):
    original_directory = os.getcwd() #starting directory. Not necessarily the cwd
    owd = "%s/histograms" % args.cwd
    
    os.chdir(args.cwd)
    fsv = open("cutoffs.txt", "w") #File for writing if the cutoff SV's
    
    if not os.path.isdir(owd):
        os.mkdir(owd)
        
    for t in args.temps:
        #Analyze and get the new histograms
        print "Starting analysis on temperature %d" % t
        os.chdir("%d"%t)
        centers_of_bins, normalized_valu, labels, highvalue, lowvalue = run_main_calc(t, args)
        #Calcualte and Plot the histograms
        os.chdir(owd)
        pdistance.plot_iterations(centers_of_bins, normalized_valu, args.pairs, labels, args.spacing, t)
        os.chdir(args.cwd)
        #write the cutoff into a file for easy access
        fsv.write("The Cutoff occurs between singular values: %.3e and %.3e\n" % (highvalue, lowvalue))
    fsv.close()
    os.chdir(original_directory)
    print "Finished all temperatures"

def run_main_calc(T_fit, args):
    #calculate for a temperature directory. Loads the arguments into the variables here
    pairs = args.pairs
    spacing = args.spacing
    subfolder = args.subdir
    pmfit.FRET.compute_Jacobian.def_temp = T_fit
    
    #internal directory tree here. Assumings everything is arranged a certain way
    cwd = os.getcwd()
    cwd0 = os.getcwd()
    cwd += "/%s" % subfolder
    
    log = "%s/modelbuilder.log" % cwd
    #load model, and change fitting method if specified
    model = mdb.check_inputs.load_model(cwd, False)
    if not args.fitting_method==None:
        model.fitting_solver = args.fitting_method
        
    rcpmanager = log_function(os.getcwd())
    pmfit.prepare_newtons_method(model,"FRET",rcpmanager.append_log)

    #This will make a histogram of all the different iterations and save the data accordingly
    os.chdir(cwd)
    centers_of_bins, normalized_valu, labels = pdistance.histogram_iterations(pairs,spacing,T_fit)
    os.chdir(cwd0)

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


def run_run_all(args):
    original_directory = os.getcwd() #starting directory. Not necessarily the cwd
    fit_string = ""
    
    #move to the cwd, and begin saving everything
    os.chdir(args.cwd)

    for t in args.temps:
        print "Starting temperature %d" % t
        os.chdir("%d"%t)
        iteration = run_main_run(int(t),args)
        ffit = open("%s/iteration_%d/newton/fitting_scale"%(args.subdir,iteration)).readline().strip()
        if not ffit == "0":
            fit_string += "Temperature %d scaled = True,  by factor = %s\n"%(t,ffit)
        else:
            fit_string += "Temperature %d scaled = False, by factor = %s\n"%(t,ffit)
        os.chdir(args.cwd)
    
    f = open("Fitting_damping_%d.txt"%iteration,"w")
    f.write(fit_string)
    f.close()
    os.chdir(original_directory)
    print "Finished all temperatures"
        
def run_main_run(T_fit, args):
    subfolder = args.subdir
    cwd = os.getcwd()
    cwd0 = os.getcwd()
    cwd += "/%s"%subfolder
    log = "%s/modelbuilder.log" % cwd
    model = mdb.check_inputs.load_model(cwd, False)

    rcpmanager = log_function(os.getcwd())
    append_log = rcpmanager.append_log
    pmfit.save_new_parameters(model,"FRET",append_log)
    
    if args.no_pbs:
        print "Not submiting pbs jobs"
        return model.iteration
    else:
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
        
        return model.iteration-1

def sanitize_args(args):
    original_directory = os.getcwd()
    os.chdir(args.cwd)
    ##If temps was not specified, will attempt to open Temparray.txt, Otherwise prints ERROR
    if args.temps==None and os.path.isfile("Temparray.txt"):
        args.temps = np.loadtxt("Temparray.txt",dtype=int)
    else:
        print "ERROR: No Temperature Directories Specified"
    
    ##set the pairs for fitting in the array format for the calculation
    pairs = np.array([[args.pairs[0], args.pairs[1]]])
    print "Number of pairs is: %d" % len(args.pairs)
    if len(args.pairs)>2:
        for i in np.arange(3, len(args.pairs), 2):
            pairs = np.append(pairs, np.array([[args.pairs[i-1], args.pairs[i]]]), axis=0)
    args.pairs = pairs
    
    #If no fitting method specified, look for the fitting.txt file to re-specify. If none, it will use the model's default method
    if args.fitting_method == None and os.path.isfile("fitting.txt"):
        f = open("fitting.txt","r")
        args.fitting_method = f.readline().strip()
        f.close()
    
    os.chdir(original_directory)
    return args
        
def get_args():
    ##parent parser for shared parameters
    parser = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    parser.add_argument("--subdir", default="1PB7", type=str, help="Sub Directory the model.info file is located in")
    parser.add_argument("--temps", default=None, type=int, nargs="+", help="Temperature Directories for Fitting. Default=load Temarray.txt")
    parser.add_argument("--fitting_method", default=None, type=str, help="Choose either TSVD (def) or Levenberg")
    parser.add_argument("--pairs", nargs="+",type=int, default=[114,192], help="pairs for FRET fitting")
    parser.add_argument("--spacing", type=float, default=0.1, help="spacing for binning the simulated and experimental FRET data")
    parser.add_argument("--cwd", type=str, default=os.getcwd(), help="Directory for fitting")
    
    ##The Real Parser
    par = argparse.ArgumentParser(description="Options for Jac_run_module. Use --cwd for analysis on not the current working directory")
    sub = par.add_subparsers(dest="step")
    
    ##Sub parsers:
    ##Calc = For running a calculation for all the Jacobians
    calc_sub = sub.add_parser("calc", parents=[parser], help="For Computing all the Jacobians for the directory and calculating their fitted value")
    
    ##Run = For calculating the next step and submitting the PBS job.
    run_sub = sub.add_parser("run", parents=[parser], help="For taking the calculated information from Calc and submitting a job with the new parameters")
    run_sub.add_argument("--no_pbs", default=False, action="store_true", help="Use if you do not want to submit the job")
    
    args = par.parse_args()
    
    args = sanitize_args(args)
    
    
    return args
if __name__ == "__main__":
    args = get_args()
    
    if args.step=="calc":
        run_calc_all(args)
    elif args.step=="run":
        run_run_all(args)
    
    print "Finished running %s" % args.step
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




