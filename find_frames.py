"""
Script for determining the frames that fit a certain parameter
Specific to dealing with the output from the dmdmd script
"""

import numpy as np
import analysis_scripts.merge_files as mf
import argparse
import os

import analysis_scripts.free_energy_plot_2d as fep

def find_frames_2D(start, stop, bounds=[0, 1, 0, 1], fout=None, groupname=None, xvg_dir=None, plot_type="QY"):
    print "Begin Frame Loading"
    
    #load the extensions
    names = {"Q":"-Qclosed.out", "A":"-rmsd-apo.xvg", "C":"-rmsd-closed.xvg", "Y":"-y114-192.out", "L":"-comA.xvg"}
    ext1 = names[plot_type[0]]
    ext2 = names[plot_type[1]]
    
    frame_num = 1
    close = False
    if fout == None:
        fout = open("iter%d-%d-%s.ndx"%(start,stop,plot_type), "w")
        close = True
    if groupname == None:
        groupname = "group1"
        
    if xvg_dir == None:
        xvg_str = ""
    else:
        xvg_str = xvg_dir
        
    fout.write("\n[%s]\n"%groupname)
    for i in np.arange(start, stop+1, 2):

        rc1 = fep.get_value("iter%d"%i, ext1, xvg_str) 
        #rc2 = np.loadtxt("%siter%d-rmsd-closed.xvg"%(xvg_str,i),  skiprows=13)
        #rc2 = rc2[:,1]
        rc2 = fep.get_value("iter%d"%i, ext2, xvg_str)
        for i in range(np.shape(rc1)[0]):
            if rc1[i]<bounds[1] and rc1[i]>bounds[0] and rc2[i]<bounds[3] and rc2[i]>bounds[2]:
                fout.write("%d\n"%frame_num)
            frame_num += 1
    if close:
        fout.close()
    print "Finished Frame Loading"
    
if __name__ == "__main__":
    par = argparse.ArgumentParser(description="parent set of parameters", add_help=False)
    par.add_argument("--file_dir", default=os.getcwd(), type=str) #location of the .gro and .w files
    par.add_argument("--xvg_dir", type=str) #location of .xvg files
    par.add_argument("--range", default=[6,100], nargs=2, type=int) 
    par.add_argument("--step", type=int)
    par.add_argument("--save_dir", default=os.getcwd(), type=str)
    par.add_argument("--bound", nargs=4, default=[0, 1, 0, 1], type =float) #list a, b, c, d as a < rc1 < b c < rc2 < d
    par.add_argument("--fname", type=str)
    par.add_argument("--append", default=False, type=bool)
    par.add_argument("--group_name", type=str)
    par.add_argument("--plot_type", type=str, default="QC", help="specify the type of plot in xy format; default QC. C=RMSD-closed, A=RMSD-apo, Q=Q, Y=FRET probe distance, L=Lobes center distance")
    
    args = par.parse_args()
    
    
    
    if args.file_dir[-1] == "/":
        args.file_dir = file_dir[:-1]
    if args.save_dir[-1] == "/":
        args.save_dir = save_dir[:-1]
    
    if not args.xvg_dir == None:
        xvg_file_dir = "%s/%s" % (args.file_dir,args.xvg_dir)
    else:
        xvg_file_dir = None
    
    if args.fname == None:
        fout = "%s/iter%d-%d-%s-frames.ndx"%(args.save_dir,args.range[0],args.range[1], args.plot_type)
        finf = "%s/iter%d-%d-%s-info.txt"%(args.save_dir,args.range[0],args.range[1], ags.plot_type)
    else:
        fout = "%s/%s.ndx"%(args.save_dir,args.fname)
        finf = "%s/%s-info.txt"%(args.save_dir,args.fname)
    
    if args.append:
        print "Appending the File"
        fsave = open(fout,"a")
        finfo = open(finf,"a")
    else:
        print "Overwriting the File"
        fsave = open(fout,"w")
        finfo = open(finf,"w")
    
    if args.step == None:
        mf.merge(args.range[0], args.range[1],args.file_dir, args.save_dir)
        find_frames_2D(args.range[0], args.range[1], bounds=args.bound, fout=fsave, groupname=args.group_name,xvg_dir=xvg_file_dir, plot_type=args.plot_type)
    else:
        if args.range[0] < args.step:
            for stop in np.arange(args.step, args.range[1], args.step):
                mf.merge(args.range[0], stop, args.file_dir, args.save_dir)
                find_frames_2D(args.range[0], stop, bounds=args.bound, fout=fsave, groupname=args.group_name,xvg_dir=xvg_file_dir, plot_type=args.plot_type)
        else:
            for stop in np.arange(args.range[0]+args.step, args.range[1], args.step):
                mf.merge(args.range[0], stop, args.file_dir, args.save_dir)     
                find_frames_2D(args.range[0], stop, bounds=args.bound, fout=fsave, groupname=args.group_name,xvg_dir=xvg_file_dir, plot_type=args.plot_type) 
    #load the extensions
    names = {"Q":"-Qclosed.out", "A":"-rmsd-apo.xvg", "C":"-rmsd-closed.xvg", "Y":"-y114-192.out", "L":"-comA.xvg"}
    ext1 = names[args.plot_type[0]]
    ext2 = names[args.plot_type[1]]
    
    label1, label2 = fep.get_labels(ext1, ext2)
    finfo.write("[%s] is for a range of %5.2f < %s <%5.2f and %5.2f < %s < %5.2f\n"%(args.group_name,args.bound[0],ext1,args.bound[1],args.bound[2],ext2,args.bound[3]))
    finfo.close() 
    fsave.close()
                
                
                
                
                
                
                
                
                
                
                
                
                
