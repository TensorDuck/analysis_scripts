"""script for taking a series of calcualted next steps, and determining at what singular value the last time a solution crosses the x-axis"""
"""assumes you are in the folder containing xp_%d files"""

import os
import numpy as np

if __name__ == "__main__":
    fwrite = open("last_crossing.txt", "w")
    cwd = os.getcwd()
    for i in np.arange(130, 171, 10):
        os.chdir("01_TSVD_all_%d"%i)
        
        last = 0
        current = 0
        go = True

        current_array = np.array([])
        last_array = np.loadtxt("xp_0.dat")

        while go:
            current += 1
            nextfile = "xp_%d.dat" % current
            if os.path.isfile(nextfile):
                current_array = np.loadtxt(nextfile)
                if np.min(current_array*last_array) < 0:
                    last = current
                last_array = current_array
            else:
                go = False
        
        singular_values = np.loadtxt("singular_values.dat")
        fwrite.write("last crossing happens for file xp_%d.dat, for singular value %f\n" % (last, singular_values[len(singular_values)-last-1]))
        os.chdir(cwd)
        


