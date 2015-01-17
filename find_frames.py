"""
Script for determining the frames that fit a certain parameter
Specific to dealing with the output from the dmdmd script
"""

import numpy as np

def find_frames_2D(start, stop, dtype="rmsdq", bounds=[0, 1, 0, 1], fout=None, groupname="group1"):
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
    find_frames_2D(6, 8, dtype="rmsdq", bounds=[750, 900, 0.4, 0.6])
            
                
