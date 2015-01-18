"""
Script for determining the frames that fit a certain parameter
Specific to dealing with the output from the dmdmd script
"""

import numpy as np
import analysis_scripts.merge_files as mf
import argparse

def find_frames_2D(start, stop, bounds=[0, 1, 0, 1], fout=None, groupname="group1"):
    frame_num = 1
    if fout == None:
        fout = open("iter%d-%d.ndx"%(start,stop), "w")
        
    fout.write("\n[%s]\n"%groupname)
    for i in np.arange(start, stop+1, 2):

        rc1 = np.loadtxt("iter%d-Q.out"%i)
        rc2 = np.loadtxt("iter%d-rmsd-closed.xvg"%i,  skiprows=13)
        rc2 = rc2[:,1]
        for i in range(np.shape(rc1)[0]):
            if rc1[i]<bounds[1] and rc1[i]>bounds[0] and rc2[i]<bounds[3] and rc2[i]>bounds[2]:
                fout.write("%d\n"%frame_num)
                print "here"
            frame_num += 1
    fout.close()
    
if __name__ == "__main__":
    par = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    par.add_argument("--file_dir", default=os.getcwd(), type=str)
    par.add_argument("--range", default=[6,100], nargs=2, type=int) 
    par.add_argument("--step", type=int)
    par.add_argument("--save_dir", default=os.getcwd(), type=str)
    par.add_argument("--bound", nargs=4, default=[0, 1, 0, 1], type =float) #list a, b, c, d as a < rc1 < b c < rc2 < d
    par.add_argument("--fname", type=str)
    par.add_argument("--append", default=False, type=bool)
    par.add_argument("--group_name", type=str)
    
    args = par.parse_args()
    
    if file_dir[-1] == "/":
        file_dir = file_dir[:-1]
    if save_dir[-1] == "/":
        save_dir = save_dir[:-1]
    
    if args.range[0] < 6:
        args.range[0] = 6
        print "Fixing range to lower bound of iteration 6"
    
    if args.name == None:
        fout = "%s/iter%d-%d-frames.ndx"%(args.save_dir,args.range[0],args.range[1])
    else:
        fout = "%s/%s.ndx"%(args.save_dir,args.fname)
    
    if args.append:
        fsave = open(fout,"a")
    else:
        fsave = open(fout,"w")
    
    if args.step == None:
        merge(args.range[0], args.range[1])
        find_frames_2D(args.range[0], args.range[1], bounds=args.bound, fout=fsave, groupname=args.group_name)
    else:
        if args.range[0] < args.step:
            for stop in np.arange(args.step, args.range[1], args.step):
                mf.merge(args.range[0], stop, args.file_dir, args.save_dir)
                find_frames_2D(args.range[0], stop, bounds=args.bound, fout=fsave, groupname=args.group_name)
        else:
            for stop in np.arange(args.range[0]+args.step, args.range[1], args.step):
                mf.merge(args.range[0], stop, args.file_dir, args.save_dir)     
                find_frames_2D(args.range[0], stop, bounds=args.bound, fout=fsave, groupname=args.group_name) 
    fsave.close()
                
