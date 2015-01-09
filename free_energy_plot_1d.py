"""
Script for plotting a 1-D free energy plot
"""

import numpy as np
import math
import matplotlib.pyplot as plt

kb = 6.022141*1.380650/(4.184*1000.0)

def make(x, nbins, T):
    
    hist, edges = np.histogram(x, bins=nbins, normed=True)
    coord_center = edges+((edges[1]-edges[0])/2.0)
    coord_center = coord_center[:-1]
    energy = np.log(hist)
    energy *= -1
    energy *= kb
    energy *= T
    return energy, coord_center
    
def plot_Q(energies, coord_centers, label, title): 
    plt.figure()

    colors = ["b","g","r","c","m","y","b","g","r","c","m","y","b","g","r","c","m","y","b","g","r","c","m","y"]
    linetype = ["-", "--",":"]
    
    plt.xlabel("Q", fontsize=20)
    plt.ylabel("F",fontsize=20)
    plt.title(title)
    for i in range(np.shape(label)[0]):
        plt.plot(coord_centers[i], energies[i], alpha=0.75, linewidth=2, linestyle=linetype[i/6], color=colors[i], label="%s"%label[i], marker="o")
    plt.legend()
    plt.savefig("%s.png"%title)

if __name__ == "__main__":
    labels = ["iter_0", "iter_10"]
    temperature = 185
    energies = []
    coord_centers = []
    nbins = 70
    title = "Q-plot-%d"%temperature
    for i in range(2):
        x = np.loadtxt("Q%d_%s.out"%(temperature,labels[i]))
        energy, coord_center = make(x, nbins, temperature)
        energies.append(energy)
        coord_centers.append(coord_center)
    plot_Q(energies, coord_centers, labels, title)
