""" Set of functions for working with FRET data

description:

remove_cut(cut,data): Remove all values in data greater than value of cut
convert_r(eff,R): Converts the efficiency values to real distancde values

Warning: R is a constant specific to the fluorescent tags used. Make sure
you use the correct constant.

reference: 
C.F. Landes et al. Nature: Chemical Biology. Vol 7, pp168-173, March 2011

"""

import numpy as np

def remove_cut(cut,data):
    """ Deletes all entries from data greater than value of cut """
    index = 0
    rmData=np.array([])
    for i in data:
        if i  > cut:
            rmData = np.append(rmData,index)
        index += 1
    data = np.delete(data,rmData)
    return data

def convert_r(eff,R):
    """Convert an efficiency value to a real distance"""
    rnew = R*(((1.0/eff)-1.0)**(1.0/6.0))
    return rnew
   
def double_gaussian(x, Ampa, ka, ya, Ampb, kb, yb):
    """Double 1-D gaussian with options
    x = independent variable
    Amp = Amplitude
    k = Exponent constant
    y = x-shift    """
    return Ampa*np.exp((-1)*ka*((x-ya)**2)) + Ampb*np.exp((-1)*kb*((x-yb)**2))
    


def chisq_reduced(obs, cal):
    """Calculated a reduced chi squared = assume sigma on observed is sqrt(N)
    Assume obs and cal are both 1-d and same size.    """
    diff = cal - obs
    count = 0
    total = 0.0
    for i in diff:
        if obs[count] == 0:
            count += 1
        else:
            total += (i**2) / obs[count]
            count += 1
    return total
    
def single_gaussian(x, Amp, k, y):
    """Single 1-D gaussian with options
    x = independent variable
    Amp = Amplitude
    k = Exponent constant
    y = x-shift    """
    return Amp*np.exp((-1)*k*((x-y)**2))
    

