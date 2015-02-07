"""
Custom script for reading in the dmdmd .gro file outputs and
producing an array in accordance to mdtraj's standard

by Justin Chen
"""
import numpy as np

class traj(object):
    def __init__(self, xyz):
        self.xyz = xyz
        self.n_atoms = np.shape(xyz)[1]

def load_frame(num_atoms, f):
    xyz_frame = np.zeros((1,num_atoms,3))
    for i in range(num_atoms):
        line = f.readline()
        values = line.split()
        xyz_frame[0,i] = np.array([float(values[3]),float(values[4]),float(values[5])])
        
    return xyz_frame

def read(fn):
    f = open(fn,"r")
    f.readline()
    line=f.readline().strip()
    num_atoms = int(line)
    xyz = load_frame(num_atoms,f)
    f.readline()
    f.readline()
    f.readline()
    xyz = np.append(xyz,load_frame(num_atoms,f), axis=0)
    line = f.readline().strip()
    cycle = True
    print "Beging Loading frames"
    while cycle:
        
        if line == "%d"%num_atoms:
            xyz = np.append(xyz,load_frame(num_atoms,f), axis=0)
        elif line == "":
            cycle = False
            print "Finished Reading the Frames, total of %d frames" % np.shape(xyz)[0]
        line = f.readline().strip()
    
    f.close()
    return xyz

if __name__=="__main__":
    xyz = read("iter98.gro")
    print np.shape(xyz)
    np.savetxt("xyz0.dat",xyz[0], fmt="%+.3f")
    np.savetxt("xyz1.dat",xyz[1], fmt="%+.3f")
    np.savetxt("xyz2.dat",xyz[2], fmt="%+.3f")
    np.savetxt("xyzlast.dat",xyz[-1,:,:], fmt="%+.3f")
    
    
    
    
    
