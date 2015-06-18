""" Quick script for plotting the time-scales versus time-lag dependence of eigenvalues"""


import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

def run_plotting(args):
    data = np.loadtxt(args.file)
    time_lags = data[:,0]
    eigenvalue = data[:,1]
    ##multiply in time_lags information for the scale
    try:
        time_lags *= args.time_scale
        unit_string = "ns"
    except:
        unit_string = "arb"
    test = np.array([eigenvalue<0])
    if not test.any():
        time_ratio = np.log(eigenvalue)
    else:
        time_ratio = np.log(np.abs(eigenvalue))
        for idx, ev in enumerate(eigenvalue):
            if not ev > 0:
                time_ratio[idx] = -time_ratio[idx]
    time_scale = -1.0 / time_ratio 
    time_scale *= time_lags
    print time_scale
    plt.figure()
    plt.plot(time_lags, time_scale, 'ok')
    plt.xlabel("Time Lag (%s)"%unit_string, fontsize=20)
    plt.ylabel("Time scale (%s)"%unit_string, fontsize=20)
    if not args.axis == None:
        plt.axis(args.axis)
    plt.savefig("%s" %args.title)
    plt.show()

def get_args():
    parser = argparse.ArgumentParser(description="Select a file")
    parser.add_argument("--file", type=str, help="Select a file where x-column is time-lag and y-column is eigenvalue")
    parser.add_argument("--time_scale", type=float, default=0.0005, help="Specify time value for one data step in ns")
    parser.add_argument("--title", type=str, help="Specify time value for one data step in ns")
    parser.add_argument("--axis", nargs=4, type=float, help="Specify time value for one data step in ns")
   
    args = parser.parse_args()
    
    if args.file == None:
        raise IOError("failed to specify a file. required!!")
    return args
    

if __name__ == "__main__":
    args = get_args()
    run_plotting(args)
