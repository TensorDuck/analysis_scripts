import numpy as np
import model_builder as mdb
import os as os
import project_tools.parameter_fitting as pmfit
import argparse
from analysis_scripts.recipe_log_function import log_function
from analysis_scripts.JacRunModule import estimate_lambda

def run_start(args):
    #internal directory tree here. Assumings everything is arranged a certain way
    cwd = os.getcwd()
    cwd0 = os.getcwd()
    cwd += "/%s" % args.subdir
    
    #load model, and change fitting method if specified
    model, fitopts = mdb.inputs.load_model(args.subdir, False)
        
    pmfit.prepare_newtons_method(model, fitopts)
    newtondir = "%s/iteration_%d/newton" % (cwd,fitopts["iteration"])
    os.chdir(newtondir)
    print "estimating the value of lambda from singular values"
    svf = np.loadtxt("singular_values.dat") 
    index = 0
    num = np.shape(svf)[0]
    """
    lowvalue = np.min(svf)
    highvalue = np.min(svf)
    for i in range(num-1):
        if (svf[i]< 0.01 and svf[i]/svf[i+1] > 1000) or svf[i+1]<10**-12:
            index = num - 1 - i
            highvalue = svf[i]
            lowvalue = svf[i+1]   

    for i in range(num):
        if svf[i] <0.01:
            index = num - i
    """
    go = True
    i = 0
    while (go and i <= num):
        if svf[i] < 0.01:
            index = num -i
            go = False
        i += 1
    open("Lambda_index.txt","w").write("%d"%index)    

    os.chdir(cwd0)    
    
def run_end(args):
    subfolder = args.subdir
    cwd = os.getcwd()
    cwd0 = os.getcwd()
    cwd += "/%s"%subfolder
    
    model, fitopts = mdb.inputs.load_model(args.subdir, False)
    
    #Both parts will fit the Jacobian per Lambda_index.txt
    pmfit.save_new_parameters(model, fitopts)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    parser.add_argument("--subdir", default="1PB7", type=str, help="Sub Directory the model.info file is located in")
    parser.add_argument("--T_fit", default=130, type=int)
    
    cwd = os.getcwd()
    
    args = parser.parse_args()
   
    run_start(args)
    
    os.chdir("%s/%s/iteration_0/newton" % (cwd, args.subdir))
    #estimate_lambda()
    
    os.chdir(cwd)
    run_end(args)
    
    
    

