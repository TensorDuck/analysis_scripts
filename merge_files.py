"""
Simple method for merging files together for the results from a dmdmd simulation
"""
import argparse
import os
import numpy as np

def merge(start, stop, cfd, csd):
    fout_name = "%s/iter%d-%d.gro"%(csd,start,stop)
    fwout_name = "%s/iter%d-%d.w"%(csd,start,stop)
    
    if os.path.isfile(fout_name) and os.path.isfile(fwout_name):
        print "Merged file already exists, aborting"
    else:
        fout = open(fout_name, "w")
        fwout = open(fwout_name, "w")
        for i in range(start, stop+1, 2):
            fin = open("%s/iter%d.gro"%(cfd,i),"r")
            fwin = open("%s/iter%d.w"%(cfd,i),"r")
            for lines in fin:
                fout.write(lines)
            for lines in fwin:
                fwout.write(lines)
            fin.close()
            fwin.close()
        fout.close()
        fwout.close()
    
if __name__ == "__main__":
    
    par = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    par.add_argument("--file_dir", default=os.getcwd(), type=str)
    par.add_argument("--range", default=[6,100], nargs=2, type=int) 
    par.add_argument("--step", type=int)
    par.add_argument("--save_dir", default=os.getcwd(), type=str)
    
    args = par.parse_args()
    
    if args.file_dir[-1] == "/":
        args.file_dir = file_dir[:-1]
    if args.save_dir[-1] == "/":
        args.save_dir = save_dir[:-1]
    
    if args.range[0] < 6:
        args.range[0] = 6
        print "Fixing range to lower bound of iteration 6"
        
    if args.step == None:
        merge(args.range[0], args.range[1])
    else:
        if args.range[0] < args.step:
            for stop in np.arange(args.step, args.range[1], args.step):
                merge(args.range[0], stop, args.file_dir, args.save_dir)
        else:
            for stop in np.arange(args.range[0]+args.step, args.range[1], args.step):
                merge(args.range[0], stop, args.file_dir, args.save_dir)
    
    
