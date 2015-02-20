import numpy as np
from densplot import hist2d
from math import ceil
import os
import matplotlib.pyplot as plt
import argparse
import free_energy_plot_2d as fep

def plot_2D_Free_Energy(rc1, rc2, rc1n, rc2n, name, plot_style, weights=None, nbins=50, axis=[0, 1000, 0, 10], temp=300):
    plt.figure()
    plot_style = "scatter"
    x, y, z = hist2d.make(rc1, rc2, nbins, nbins, temperature=temp, weight=weights, plot_style=plot_style, free_energy_plot=True, idx_smoothing=3)
    cp = plt.scatter(x, y, s=10, c=z, marker='o', linewidth=0.)
    #vmax=7
    if not axis == None:
        plt.axis(axis)
    plt.xlabel(rc1n, size=16)
    plt.ylabel(rc2n, size=16)

    cb = plt.colorbar(cp)
    pad = 10
    cb.set_label(r'Free Energy (kT units)',labelpad=pad)

    plt.savefig("%s_%s-%s-%s.png"%(name, rc1n, rc2n, plot_style))
    plt.close()

def run_plot_iterations(start, stop, file_location=None, output_location=None, w=False, nbins=50, axisr=None, axisq=None, axisy=None, rmsd_plot=True, Q_plot=True, scatter_only=False):
    print "Starting iteration range from %d to %d" % (start, stop)
    #assumes rc1 = rmsd apo, rc2 = rmsd closed
    cwd = os.getcwd()
    rc1n = "rmsd-apo (nm)"
    rc2n = "rmsd-closed (nm)"
    rcqn = "Q-1PB7"
    rcyn = "y-115-193(nm)"
    
    if file_location == None:
        cfd = cwd
    else:
        cfd = file_location
    if output_location == None:
        csd = "%s/plots" % cwd
    else:
        csd = output_location
    if not os.path.isdir(csd):
        os.mkdir(csd)
    
    #set final matrices, to append all the data onto
    #rc1 = np.array([])
    #rc2 = np.array([])
    #rcq = np.array([])
    
    for i in np.arange(130, 171, 10):
        for j in np.arange(start, stop+1, 1):
            f1 = np.loadtxt("%s/T%dI%d-rmsd-apo.xvg"%(cfd,i,j), skiprows=13)
            f2 = np.loadtxt("%s/T%dI%d-rmsd-closed.xvg"%(cfd,i,j), skiprows=13)
            fq = np.loadtxt("%s/T%dI%d-Qclosed.out"%(cfd,i,j))
            fy = np.loadtxt("%s/T%dI%d-y114-192.out"%(cfd,i,j))
            rc1 = f1[:,1]
            rc2 = f2[:,1]
            rcq = fq
            rcy = fy
            os.chdir(csd)
            plot_2D_Free_Energy(rcq, rc2, rcqn, rc2n, "T%d-I%d"%(i,j), "scatter", weights=None, nbins=nbins, axis=axisq,temp=i)
            plot_2D_Free_Energy(rc1, rc2, rc1n, rc2n, "T%d-I%d"%(i,j), "scatter", weights=None, nbins=nbins, axis=axisr,temp=i)
            plot_2D_Free_Energy(rcq, rcy, rcqn, rcyn, "T%d-I%d"%(i,j), "scatter", weights=None, nbins=nbins, axis=axisy,temp=i)
            os.chdir(cwd)
    '''       
    for i in np.arange(170, 186, 5):
        for j in np.arange(start, stop+1, 1):
            f1 = np.loadtxt("%s/T%dI%d-rmsd-apo.xvg"%(cfd,i,j), skiprows=13)
            f2 = np.loadtxt("%s/T%dI%d-rmsd-closed.xvg"%(cfd,i,j), skiprows=13)
            fq = np.loadtxt("%s/T%dI%d-Qclosed.out"%(cfd,i,j))
            fy = np.loadtxt("%s/T%dI%d-y114-192.out"%(cfd,i,j))
            rc1 = f1[:,1]
            rc2 = f2[:,1]
            rcq = fq
            rcy = fy
            os.chdir(csd)
            fep.plot_2D_Free_Energy(rcq, rc2, rcqn, rc2n, "T%d-I%d"%(i,j), "scatter", weights=None, nbins=nbins, axis=axisq,temp=i)
            fep.plot_2D_Free_Energy(rc1, rc2, rc1n, rc2n, "T%d-I%d"%(i,j), "scatter", weights=None, nbins=nbins, axis=axisr,temp=i)
            fep.plot_2D_Free_Energy(rcq, rcy, rcqn, rcyn, "T%d-I%d"%(i,j), "scatter", weights=None, nbins=nbins, axis=axisq,temp=i)
            os.chdir(cwd)
    '''
if __name__ == "__main__":
    run_plot_iterations(0, 5, output_location="plots_same", axisr=[0, 1, 0 ,1], axisq=[0, 1000, 0, 12], axisy=[0, 1000, 0, 20])
