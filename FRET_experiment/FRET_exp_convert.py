""" A recipe for converting the experimental Efficiency(time) trace to R(time)
Data is from Landes group, should be formatted as a trace.

Efficiency will be converted directly to R, and then a histogram will be plotted

A cutoff efficiency Ecut will be used for efficiency values near 1.

Units are in nanometers

Assumes these files are present in the current directory:
FRET_obsf.txt : Observed FRET efficinecies time trace file
FRET_denf.txt : De-noised FRET efficiencies time trace file



"""
import numpy as np
import matplotlib.pyplot as plt
import FRET_experiment.FRET_exp_pack as fep
from scipy.stats import norm

Ro = 5.1   #5.1nm constant Ro, dependent on FRET probe pairs   
Ecut = 1   #Set to 1, anything larger will give an imaginary solution.
bin_size = 100.0    #number of bins
ran_size = (0,10)    #range of values for the final distance
step = (float(ran_size[1]-ran_size[0])) / bin_size
start = ran_size[0] + (step/2)
end = ran_size[1] - (step/2)

##variables for outputting a file with the final gaussian data 
xspace = np.linspace(start, end, bin_size)
print xspace
dist_data = []

Obs_trace = np.loadtxt("FRET_obsf.txt", delimiter=",")
Den_trace = np.loadtxt("FRET_denf.txt", delimiter=",")

Obs_trace = fep.remove_cut(Ecut, Obs_trace)
Den_trace = fep.remove_cut(Ecut, Den_trace)

Obs_r = np.array([])
Den_r = np.array([])

##Convert efficincies to distance R
for i in Obs_trace:
    Obs_r = np.append(Obs_r, fep.convert_r(i,Ro))

for i in Den_trace:
    Den_r = np.append(Den_r, fep.convert_r(i,Ro))
    
np.savetxt("FRET_trace.dat", Den_r)

##Begin Plotting the Data for niceties

plt.figure(1)
n, bins, patches = plt.hist(Obs_r, bins=bin_size, range=ran_size, facecolor='green', alpha=0.8)
plt.xlabel('R (nm)',fontsize=20)
plt.ylabel('Total Found',fontsize=20)
plt.title('Observed FRET: R distribution',fontsize=20)
plt.axis([0,8,0,8000])
plt.savefig("Observed-R-distribution.png")

dist_data = np.transpose(np.array([xspace,n]))
np.savetxt('Observed-R-distribution.txt', dist_data)



plt.figure(2)
n, bins, patches = plt.hist(Den_r, bins=bin_size, range=ran_size, facecolor='green', alpha=0.8)
plt.xlabel('R (nm)',fontsize=20)
plt.ylabel('Total Found',fontsize=20)
plt.title('De-noised FRET: R distribution',fontsize=20)
plt.axis([0,8,0,8000])
plt.savefig("Denoised-R-distribution.png")

dist_data = np.transpose(np.array([xspace,n]))
np.savetxt('Denoised-R-distribution.txt', dist_data)


plt.figure(3)
n, bins, patches = plt.hist(Obs_trace, bins=100, range=(0,1), facecolor='green', alpha=0.9)
plt.xlabel('R (nm)',fontsize=20)
plt.ylabel('Total Found',fontsize=20)
plt.title('Observed FRET: R distribution',fontsize=20)
plt.axis([0,8,0,8000])
plt.savefig("Observed-Eff-distribution.pdf")


plt.figure(4)
n, bins, patches = plt.hist(Den_trace, bins=100, range=(0,1), facecolor='green', alpha=0.9)
plt.xlabel('R (nm)',fontsize=20)
plt.ylabel('Total Found',fontsize=20)
plt.title('Observed FRET: R distribution',fontsize=20)
plt.axis([0,8,0,8000])
plt.savefig("Denoised-Eff-distribution.pdf")


plt.show()
