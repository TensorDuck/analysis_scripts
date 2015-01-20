import numpy as np
from densplot import hist2d
from math import ceil
import os
import matplotlib.pyplot as plt
import argparse
import free_energy_plot_2d as fep

def run_plot_iterations(start, stop, file_location=None, output_location=None, w=False, nbins=50, axisr=None, axisq=None, rmsd_plot=True, Q_plot=True, scatter_only=False):
    print "Starting iteration range from %d to %d" % (start, stop)
    #assumes rc1 = rmsd apo, rc2 = rmsd closed
    cwd = os.getcwd()
    rc1n = "rmsd-apo (nm)"
    rc2n = "rmsd-closed (nm)"
    rcqn = "Q-1PB7"
    
    if file_location == None:
        cfd = cwd
    else:
        cfd = file_location
    if output_location == None:
        if not os.path.isdir("plots"):
            os.mkdir("plots")
        csd = "%s/plots" % cwd
    else:
        csd = output_location
    
    #set final matrices, to append all the data onto
    rc1 = np.array([])
    rc2 = np.array([])
    rcq = np.array([])
    
    for i in np.arange(130, 170, 10):
        for j in np.arange(start, stop+1, 1):
            f1 = np.loadtxt("%s/T%dI%d-rmsd-apo.xvg"%(cfd,i,j), skiprows=13)
            f2 = np.loadtxt("%s/T%dI%d-rmsd-closed.xvg"%(cfd,i,j), skiprows=13)
            fq = np.loadtxt("%s/T%dI%d-Qclosed.out"%(cfd,i,j))
            rc1 = f1[:,1]
            rc2 = f2[:,1]
            rcq = fq
            os.chdir(csd)
            fep.plot_2D_Free_Energy(rcq, rc2, rcqn, rc2n, "T%d-I%d"%(i,j), "scatter", weights=None, nbins=nbins, axis=axisq,temp=i)
            fep.plot_2D_Free_Energy(rc1, rc2, rc1n, rc2n, "T%d-I%d"%(i,j), "scatter", weights=None, nbins=nbins, axis=axisr,temp=i)
            os.chdir(cwd)
            
    for i in np.arange(170, 186, 5):
        for j in np.arange(start, stop+1, 1):
            f1 = np.loadtxt("%s/T%dI%d-rmsd-apo.xvg"%(cfd,i,j), skiprows=13)
            f2 = np.loadtxt("%s/T%dI%d-rmsd-closed.xvg"%(cfd,i,j), skiprows=13)
            fq = np.loadtxt("%s/T%dI%d-Qclosed.out"%(cfd,i,j))
            rc1 = f1[:,1]
            rc2 = f2[:,1]
            rcq = fq
            os.chdir(csd)
            fep.plot_2D_Free_Energy(rcq, rc2, rcqn, rc2n, "T%d-I%d"%(i,j), "scatter", weights=None, nbins=nbins, axis=axisq,temp=i)
            fep.plot_2D_Free_Energy(rc1, rc2, rc1n, rc2n, "T%d-I%d"%(i,j), "scatter", weights=None, nbins=nbins, axis=axisr,temp=i)
            os.chdir(cwd)
    
if __name__ == "__main__":
    run_plot_iterations(0, 16, output_location="plots/same_axis_asDMDMD", axisr=[0, 1, 0 ,1], axisq=[0, 1000, 0, 1])
