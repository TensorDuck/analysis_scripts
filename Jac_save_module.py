"""
This is a script for testing the functionality of the compute_Jacobian module for FRET

"""

import numpy as np
#import project_tools.parameter_fitting.FRET.compute_Jacobian as compJ
import model_builder as mdb
import os as os
import project_tools.parameter_fitting as pmfit
from FRET_experiment.MandC2004_hacked import MandC2004hack
from project_tools import simulation
import matplotlib.pyplot as plt

'''
cwd = os.getcwd()

cwd += "/1PB7"

model = mdb.check_inputs.load_model(cwd, False)

bs, rs, sp = compJ.get_sim_params(model)
bincenters = compJ.get_sim_centers(model)
compJ.fret_hist_calc(model, bs, rs, sp) 

target, targeterr = compJ.get_target_feature(model)

print target
print compJ.get_sim_array(model)
print targeterr

simf, simferr, Jac, Jacerr = compJ.calculate_average_Jacobian(model)

np.savetxt("Jac_testB.dat", Jac)

'''

def plot_Lcurve_curvature(nrm_resd,nrm_soln,Lambdas):
    x = np.log10(nrm_resd)
    y = np.log10(nrm_soln)
    y2 = np.log10(Lambdas)
    ## Calculate the second derivative of the L-curve.
    x2 = np.array([ 0.5*(x[i] + x[i+1]) for i in range(len(x)-1)])
    x3 = np.array([ 0.5*(x2[i] + x2[i+1]) for i in range(len(x2)-1)])
    
    dy = np.diff(y)
    dx = np.diff(x)
    dx2 = np.diff(x2)
    dydx = dy/dx
    ddy = np.diff(dydx)
    ddyddx = ddy/dx2
    
    plt.plot(x,y2,'r')
    plt.xlabel("$\\log{||J\\delta\\epsilon_{\\lambda} - \\delta f||}$",fontsize=16)
    plt.ylabel("$\\log{\\lambda}$",fontsize=16)
    plt.show()
    
    ## Plot the second derivative of the L-curve to help the user select the
    ## optimum damping parameter, lambda_choice. lambda_choice should be 
    ## the lambda where the L-curve curvature is the most positive.
    fig, ax1 =  plt.subplots()
    ax2 = ax1.twinx()
    line1, = ax1.plot(x3,ddyddx,'b')
    ax1.xaxis.grid(True)
    ax1.yaxis.grid(True)
    line2, = ax2.plot(x,y2,'r')
    ax1.set_xlabel("$\\log{||J\\delta\\epsilon_{\\lambda} - \\delta f||}$",fontsize=16)
    ax1.set_ylabel("$\\frac{d^2}{dx^2}\\log{||\\delta\\epsilon_{\\lambda}}||}$",fontsize=16)
    ax2.set_ylabel("$\\log{\\lambda}$",fontsize=16)
    plt.title("$\\lambda_{choice}$ should be $\\lambda$ of max positive curvature")

    ax1.yaxis.label.set_color(line1.get_color())
    ax2.yaxis.label.set_color(line2.get_color())
    ax1.tick_params(axis='y', colors=line1.get_color())
    ax2.tick_params(axis='y', colors=line2.get_color())
    
def run_main():
    cwd = os.getcwd()
    cwd0 = os.getcwd()
    cwd += "/1PB7"
    log = "%s/modelbuilder.log" % cwd
    T_fit = int(np.loadtxt("fitting_temperature.txt"))
    model = mdb.check_inputs.load_model(cwd, False)

    rcpmanager = MandC2004hack(os.getcwd())
    append_log = rcpmanager.append_log
    #pmfit.prepare_newtons_method(model,"FRET",rcpmanager.append_log)
    pmfit.save_new_parameters(model,"FRET",append_log)
    
    
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
    
if __name__ == "__main__":
    run_main()
    
    print "GOT TO END"



