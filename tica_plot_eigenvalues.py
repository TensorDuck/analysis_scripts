import numpy as np
import pyemma
import pyemma.coordinates as coor
import mdtraj as md
import time
import analysis_scripts.plot_package as pltpkg


def run_plot(lag_times, collected_eigenvalues, eigen_number):
    title = "lag-plot-eigenvalue-%d" % eigen_number
    pltpkg.plot_simple(lag_times, collected_eigenvalues, title, "lag time", "eigenvalue") 
    
    ##save a txt file of all eigenvalues it's found, as well as the figure
    np.savetxt("lag-sequence-eigen-%d.dat"%eigen_number, np.array(lag_times, collected_eigenvalues).transpose())


topology = "Native.pdb"
ticadim = 10

feat = coor.featurizer(topology)

pair = []
cutoff = 5
start = 5
stop = 287
for i in np.arange(start, stop, 1):
    for j in np.arange(i+1, stop, cutoff):
        pair.append([i, j])
    
print np.shape(pair)

pairs = np.array(pair)

feat.add_distances(pairs)
#feat.add_distances_ca()
X1 = coor.load("traj.xtc", feat, stride=1)

#traj = md.load("traj.xtc", top="Native.pdb")
#X1 = md.compute_distances(traj, [[115, 192]], periodic=False)

print np.shape(X1)
possible_times = np.logspace(0,3,50)
possible_times = possible_times.astype(int)
lag_times = []
for i in possible_times:
    if i not in lag_times:
        lag_times.append(i)

print lag_times

collected_eigenvalues=[]
for i in range(ticadim):
    collected_eigenvalues.append([])
for i in lag_times:
    tica_obj = coor.tica(X1, stride=1, lag=i, dim=ticadim)
    outputs = tica_obj.get_output()[0]
    eigen = tica_obj.eigenvalues
    np.savetxt("output_L%d.dat"%i, outputs)
    np.savetxt("eigenvalues_L%d.dat"%i, eigen)
    for j in range(ticadim):
        collected_eigenvalues[j].append(eigen[j])
        run_plot(lag_times[:len(collected_eigenvalues[j])], collected_eigenvalues[j], j)
        
        
    
    

