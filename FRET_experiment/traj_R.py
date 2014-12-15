""" A class that contains a histogram of R pair distances
and functions for analysing such data.


"""

import numpy as np

class traj_R(object):
    def __init__(self, p, bs, rs):
        print "There are four lights"
        self.bin_size = bs    #number of bins
        self.ran_size = rs    #range of values for the final distance
        self.pairs = p
        step = (float(self.ran_size[1]-self.ran_size[0])) / (2*self.bin_size)
        start = self.ran_size[0] + step
        end = self.ran_size[1] - step 
        self.xspace = np.linspace(start, end, self.bin_size)
        print "Shut Up Wesley."
        self.hist_data=[]
        self.print_once = 0
    
    def erase_all(self):
        self.hist_data = []
        
    def normalize_all(self):
        if not self.hist_data == []:
            for w in self.hist_data:
                for q in w:
                    if self.print_once == 0:
                        print q
                        self.print_once = 1
                    N = np.sum(q[:,1])
                    q[:,1] = q[:,1] / N        
        else:
            print "There is no histogram data loaded"
