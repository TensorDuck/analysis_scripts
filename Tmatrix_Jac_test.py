"""
Test code for transition matrix Jacobian

Last updated: June 15, 2015
"""

import numpy as np
import scipy.stats as stats

# CHECK TO SEE HOW LONG IT TAKES TO RUN
# import time
# start_time = time.time()

def find_FRET_bins(FRETr):
	# Histograms the trace data into macrostate bins, analogous to find_sim_bins from compute_Jacobian script
	
	weights = np.ones(np.shape(FRETr)[0])
	spacing = 0.1
	
	# Taken from find_sim_bins, included for consistency
	maxvalue = int(np.amax(FRETr)/spacing) + 1
	minvalue = int(np.amin(FRETr)/spacing)
	num_bins = maxvalue - minvalue
	ran_size = (minvalue*spacing,maxvalue*spacing)
	
	hist, bins, slices = stats.binned_statistic(FRETr, weights, statistic="sum", range=[ran_size], bins=num_bins)
	
	# Call get_transition_bins to make transition bin indices
	F_indices, t_indices = get_transition_bins(slices, num_bins)
	
	return hist, F_indices, t_indices
	
def get_transition_bins(slices, num_bins):
	# Possible new function, or could add into find_sim_bins
	# Returns linear indices of "transition bins," should match with numpy.ndarray.flatten() results
	# 'slices' has dimension Nx1, where N=number of frames
	
	# Get i and j slices, i=macrostate before transition, j=macrostate after transition
	# Subtract one to convert to zero-based indices
	F_indices = slices - 1
	state_i = F_indices[:-1]
	state_j = F_indices[1:]
	
	# Compute bin indices for transitions
	t_indices = state_i*num_bins + state_j
	
	return F_indices, t_indices
	
def Jacobian_test(qij, hist, F_indices, t_indices):
	# Would-be modified compute_Jacobian_basic
	# qij has dimension NxM
	# hist = histogram of FRET distances
	# F_indices = FRET bin indices
	# t_indices has dimension (N-1)x1
	# N = number of frames, M = number of contact pairs
	
	# number of FRET bins
	nstates = np.size(hist)
	
	# number of TRANSITION bins
	nbins = nstates**2
		
	(N_total_traj, npairs) = np.shape(qij)
	
	qi = np.zeros((nstates, npairs))
	qi_count = np.zeros(nstates)
	##PRINT THINGS FOR DEBUGGING
	print F_indices
	print np.shape(F_indices)
	print np.shape(qi)
	print np.shape(qij)
	##
	
	# Q values for all frames starting in a particular bin
	for idx, F_bin_location in enumerate(F_indices[:-1]):
		qi[F_bin_location,:] += (qij[idx,:] + qij[idx+1,:])
		qi_count[F_bin_location] += 1
	
	#normalize the average value for the pair sum starting in state i
	for i in range(np.shape(qi)[0]):
	    qi[i,:] /= float(qi_count[i])
	    
	Jacobian = np.zeros((nbins, npairs))
	for idx, t_bin_location in enumerate(t_indices):
		# Add q values for specific transition
		Jacobian[t_bin_location, :] += (qij[idx,:] + qij[idx+1,:])	
	
	# Normalize by i bin count
	for i in range(np.shape(Jacobian)[0]):
		Jacobian[i,:] /= (2*qi_count[np.floor(i/nbins)])
		
	for idx, t_bin_location in enumerate(t_indices):
		# Index for q value of all transitions starting at state i
		state_i_idx = np.floor(t_bin_location/nbins)
		
		# Subtract q values for all starting at state i
		Jacobian[t_bin_location, :] -= qi[state_i_idx,:]
				
	return Jacobian
	
if __name__=="__main__":

	# Number of frames and contacts for testing purposes
	frames = 100000
	contacts = 980
	
	# FRET trace (all frames) in actual implementation
	trace = np.linspace(0, 0.999999, frames)
	trace = trace*10
	
	# Q values, would be 'qij' in compute_Jacobian (in real life, this comes from simulation data)
	qij = np.ones((frames, contacts))
	qij = -qij
	
	hist, F_indices, t_indices = find_FRET_bins(trace) # would normally go in compute_Jacobian_basic
	
	Jacobian = Jacobian_test(qij, hist, F_indices, t_indices)
	np.savetxt("Jacobian.dat", Jacobian)
	
	# PRINT STUFF	
#	print hist
#	print np.shape(hist)
	
	print F_indices
	print np.shape(F_indices)
	
#	print qij
#	print np.shape(qij)
	
	print t_indices
	print np.shape(t_indices)
	print np.max(t_indices)	
	
#	print Jacobian
	print np.shape(Jacobian)
	
	# ELAPSED TIME TO COMPUTE JACOBIAN
#	print("--- %s seconds ---" % (time.time() - start_time))
	
