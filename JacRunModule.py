"""
This is a script for calculating the jacobin and FRET bins of a protein

It performs all the steps of: calculating Jacobian and new parameters, saving new parameters for the next model, and submitting new simulation job

"""

import numpy as np
import matplotlib.pyplot as plt
import os as os
import argparse

import model_builder as mdb
import project_tools.parameter_fitting as pmfit
from project_tools import simulation

import analysis_scripts.pair_distance_calculator as pdistance

def run_calc_all(args):
    original_directory = os.getcwd() #starting directory. Not necessarily the cwd

    os.chdir(args.cwd)
    owd = "%s/histograms" % os.getcwd()
    
    if args.single:
        print "Starting Calc on one temperature"
        model, fitopts = mdb.inputs.load_model(args.subdir, False)
        
        pmfit.prepare_newtons_method(model, fitopts)
        newtondir = "%s/iteration_%d/newton" % (args.subdir,fitopts["iteration"])
        os.chdir(newtondir)
        print "estimating the value of lambda from singular values"
        svf = np.loadtxt("singular_values.dat") 
        index = 0
        num = np.shape(svf)[0]
        if "truncate_value" in fitopts:
            trunc = fitopts["truncate_value"]
        else:
            trunc = 0.1    
        highvalue, lowvalue, lambda_index = estimate_lambda(trunc)    
        os.chdir(args.cwd)
        print "Finished Calc on one temperature"
    else:
        print "Starting Calc on all temperatures"
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
            pdistance.plot_iterations(centers_of_bins, normalized_valu, args.pairs, labels, args.spacing, t, fretdata=args.fret_data)
            os.chdir(args.cwd)
            #write the cutoff into a file for easy access
            fsv.write("The Cutoff occurs between singular values: %.3e and %.3e\n" % (highvalue, lowvalue))
        fsv.close()
        print "Finished Calc on all temperatures"
    os.chdir(original_directory)
    

def run_main_calc(T_fit, args):
    #calculate for a temperature directory. Loads the arguments into the variables here
    pairs = args.pairs
    spacing = args.spacing
    subfolder = args.subdir
    
    #internal directory tree here. Assumings everything is arranged a certain way
    cwd = os.getcwd()
    cwd0 = os.getcwd()
    cwd += "/%s" % subfolder
    
    #load model, and change fitting method if specified
    model, fitopts = mdb.inputs.load_model(args.subdir, False)
    if not args.fitting_method==None:
        fitopts["solver"] = args.fitting_method
        
    pmfit.prepare_newtons_method(model, fitopts)
    
    #This will make a histogram of all the different iterations and save the data accordingly
    os.chdir(cwd)
    centers_of_bins, normalized_valu, labels = pdistance.histogram_iterations(pairs,spacing,T_fit, fitopts)
    print labels
    os.chdir(cwd0)

    ##Following will estimate the expected cutoff, along with return the boundaries of the cutoff
    newtondir = "%s/iteration_%d/newton" % (cwd,fitopts["iteration"])
    os.chdir(newtondir)
    # Use truncate value from fitopts, otherwise use default of 0.01
    if "truncate_value" in fitopts:
        trunc = fitopts["truncate_value"]
    else:
        trunc = 0.1    
    highvalue, lowvalue, lambda_index = estimate_lambda(trunc)
    os.chdir(cwd0)    
    
    fitopts["last_completed_task"] = "Finished: Solving_Newtons_Method"
    
    return centers_of_bins, normalized_valu, labels, highvalue, lowvalue

def estimate_lambda(trunc):
    #will read the singular values file and find where a alrge jump exists
    #assumes singular_values.dat is in the directory this is run in
    print "estimating the value of lambda from singular values"    
    lvalues = np.loadtxt("lambdas.dat") 
    svf = np.loadtxt("singular_values.dat")
    index = 0
    max_search = np.shape(svf)[0]
    lowvalue = np.min(svf)
    highvalue = np.min(svf)
    go = True
    i = 0 
    while (go and i < max_search):
        if test_truncate(lvalues[i], svf, trunc):
            go = False
            lowvalue = svf[i]
            highvalue = svf[i+1]
        i += 1
    open("Lambda_index.txt","w").write("%d"%(i-1))
    return highvalue, lowvalue, index

def test_truncate(lam, svf, trunc):
    if svf[0] > lam:
        return False
    elif svf[np.shape(svf)[0]-1] < lam:
        raise IOError("Larges singular value is smaller than the truncate value, consider changing")
    else:    
        for i in range(np.shape(svf)[0]-1):
            if svf[i] < lam and svf[i+1] >= lamb:
                high = svf[i+1]
                low = svf[i]
        if high >= trunc and low <trunc:
            return True
        else:
            return False
        
def run_save_all(args):
    original_directory = os.getcwd() #starting directory. Not necessarily the cwd
    fit_string = ""
    
    #move to the cwd, and begin saving everything
    os.chdir(args.cwd)
    if args.single:
        print "Starting Save on one temperature"
        iteration = run_main_save(args.temps[0],args)
        print "Finished Save on one temperature"
    else:
        print "Running Save on all temperatures"
        for t in args.temps:
            #Cycle through temperatures, and calculate the ite
            print "Starting temperature %d" % t
            os.chdir("%d"%t)
            iteration = run_main_save(t,args)
            #ffit will write out a file detailing which frames are scaled and which are not
            fit_data_file = open("%s/iteration_%d/newton/fitting_scale"%(args.subdir,iteration))
            ffit = fit_data_file.readline().strip()
            eps_average = fit_data_file.readline().strip()
            if not float(ffit) == 1:
                fit_string += "Temperature %d scaled = True,  by factor = %s"%(t,ffit)
            else:
                fit_string += "Temperature %d scaled = False, by factor = %s"%(t,ffit)
            fit_string += ", eps_average = %s\n" % eps_average
            os.chdir(args.cwd)
        
        f = open("fitting_info_%d.txt"%iteration,"w")
        f.write(fit_string)
        f.close()
        print "Finished Save on all temperatures"
    os.chdir(original_directory)
    
        
def run_main_save(T_fit, args):
    subfolder = args.subdir
    cwd = os.getcwd()
    cwd0 = os.getcwd()
    cwd += "/%s"%subfolder
    
    model, fitopts = mdb.inputs.load_model(args.subdir, False)
    
    #Both parts will fit the Jacobian per Lambda_index.txt
    pmfit.save_new_parameters(model, fitopts)
    
    #args.pbs If true: submit a new job with new parameters, if False, return iteration number
    if args.pbs:
        fitopts["iteration"] += 1
        newdirec = "%s/iteration_%d" % (cwd, fitopts["iteration"])
        
        if not os.path.isdir(newdirec):
            os.mkdir(newdirec)
        
        os.chdir(newdirec)
        
        simulation.constant_temp.run_temperature_array(model,fitopts,T_fit,T_fit,5)
        
        os.chdir(cwd0)
        mdb.inputs.save_model(model, fitopts)
        
        return fitopts["iteration"]-1
    else:
        print "Not submiting pbs jobs"
        return fitopts["iteration"]
        
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
    if args.single:
        if not len(args.temps) == 1:
            raise IOError("There are too many temperatures specified with the single-flag. Either speciy only one temperature in temps, or do not use the single flag")
    
    os.chdir(original_directory)
    return args
        
def get_args():
    ##parent parser for shared parameters
    parser = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    parser.add_argument("--subdir", default="1PB7", type=str, help="Sub Directory, also the name of the .ini file")
    parser.add_argument("--temps", default=None, type=int, nargs="+", help="Temperature Directories for Fitting. Default=load Temarray.txt")
    parser.add_argument("--cwd", type=str, default=os.getcwd(), help="Directory for fitting")
    parser.add_argument("--pairs", nargs="+",type=int, default=[114,192], help="pairs for FRET fitting")
    parser.add_argument("--fitting_method", default=None, type=str, help="Choose either TSVD (def) or Levenberg")
    parser.add_argument("--spacing", type=float, default=0.1, help="spacing for binning the simulated and experimental FRET data")
    parser.add_argument("--fret_data", type=str, default="den", help="specify the type of FRET data using. Either den=Denoised or obs=Observed")
    parser.add_argument("--single", default=False, action="store_true", help="Specify you are working with only one temperature or file")
    parser.add_argument("--title", default="test", type=str, help="specify a title for save files for certain auto-generated files. Generally not needed")
    
    ##The Real Parser
    par = argparse.ArgumentParser(description="Options for Jac_run_module. Use --cwd for analysis on not the current working directory")
    sub = par.add_subparsers(dest="step")
    
    ##Sub parsers:
    ##Calc = For running a calculation for all the Jacobians
    calc_sub = sub.add_parser("calc", parents=[parser], help="For Computing all the Jacobians for the directory and calculating their fitted value")
    
    ##Run = For calculating the next step and submitting the PBS job.
    save_sub = sub.add_parser("save", parents=[parser], help="For taking the calculated information from calc and submitting a job with the new parameters")
    save_sub.add_argument("--pbs", default=False, action="store_true", help="Use if you do not want to submit the job")
    
    #run_sub = sub.add_parser("run", parents=[calc_sub], help="Run both calc and save blind", conflict_handler="resolve")
    
    args = par.parse_args()
    
    args = sanitize_args(args)
    
    
    return args
if __name__ == "__main__":
    args = get_args()
    
    if args.step=="calc":
        run_calc_all(args)
    elif args.step=="save":
        run_save_all(args)
    #elif args.step=="run":
        #run_calc_all(args)
        #run_save_all(args)
    print "Finished running %s" % args.step
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




