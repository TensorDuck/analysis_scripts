""" Method for fitting an outputted distribution of R values
from the FRET_exp_convert.py script

Attempts to fit multiple gaussians, or a singular gaussian, then plot

"""

import numpy as np
import FRET_experiment.FRET_exp_pack as fep
import scipy.optimize as scopt
import matplotlib.pyplot as plt
import scipy.stats as stat

##load files, formatted: x = [0,:], y = [1,:]
Obs_dist = np.loadtxt("Observed-R-distribution.txt")
Den_dist = np.loadtxt("Denoised-R-distribution.txt")

##Initial guesses, need to do these manually
Obs_initial = [5000.0, 0.5, 4.0]
Den_initial = [7000.0, 0.5, 4.0]

##Store final optimized coordinates
Obs_optimal = []
Den_optiaml = []

##variables for plotting the final fitted curves
xspace = np.linspace(0,10,1001)
print xspace
yobs = []
yden = []

##Perform the curve fitting
Obs_optimal, obscov = scopt.curve_fit(fep.single_gaussian, Obs_dist[:,0], Obs_dist[:,1], Obs_initial)
Den_optimal, dencov = scopt.curve_fit(fep.single_gaussian, Den_dist[:,0], Den_dist[:,1], Den_initial)

##Calculated a Chi-Squared Test on the fitted gaussians
for n in Obs_dist[:,0]:
    yobs = np.append(yobs, fep.single_gaussian(n, Obs_optimal[0], Obs_optimal[1], Obs_optimal[2]))
for n in Den_dist[:,0]:
    yden = np.append(yden, fep.single_gaussian(n, Den_optimal[0], Den_optimal[1], Den_optimal[2]))

print "The chisquare from scipy is = ", stat.chisquare(Obs_dist[:,1], yobs)
obs_chi = fep.chisq_reduced(Obs_dist[:,1], yobs)
den_chi = fep.chisq_reduced(Den_dist[:,1], yden)

yobs = []
yden = []

##Plot the stuff
for n in xspace:
    yobs = np.append(yobs, fep.single_gaussian(n, Obs_optimal[0], Obs_optimal[1], Obs_optimal[2]))
    yden = np.append(yden, fep.single_gaussian(n, Den_optimal[0], Den_optimal[1], Den_optimal[2]))

#observed
plt.figure(1)
plt.bar(Obs_dist[:,0], Obs_dist[:,1], width=0.1, color="green", alpha=0.5, linewidth=0, align="center")
plt.plot(xspace, yobs, color='k', linewidth=2)
plt.xlabel('R (nm)')
plt.ylabel('Total Found')
plt.title("Observed FRET: Chi2 = %f" % obs_chi)
plt.savefig("Fitted-Observed-R-distribution.pdf")

#Denoised
plt.figure(2)
plt.bar(Den_dist[:,0], Den_dist[:,1], width=0.1, color="green", alpha=0.5, linewidth=0, align="center")
plt.plot(xspace, yden, color='k', linewidth=2)
plt.xlabel('R (nm)')
plt.ylabel('Total Found')
plt.title('Denoised FRET: Chi2 = %f' %den_chi)
plt.savefig("Fitted-Denoised-R-distribution.pdf")


plt.show()