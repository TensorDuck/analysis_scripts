"""
Test code for transition matrix Jacobian

Last updated: June 15, 2015
"""

import numpy as np
import scipy.stats as stats

# CHECK TO SEE HOW LONG IT TAKES TO RUN
# import time
# start_time = time.time()

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

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
	
	return hist, F_indices, t_indices, num_bins
	
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
	
def get_T_matrix(FRET_trace):
	# Calculate flattened transition matrix for a given FRET trace
	# Based on fret_analysis/compute_transitions.py
	
	# Get FRET bins
	hist, F_indices, t_indices, num_bins = find_FRET_bins(FRET_trace)
	
	T_matrix = np.zeros((num_bins, num_bins))	
	
	# Add ones to transition bins in square transition matrix
	for i in range(np.shape(F_indices)[0] - 1):
		T_matrix[F_indices[i], F_indices[i+1]] += 1
	
	# Mask zeros to avoid divide-by-zero in normalization
	T_masked = np.ma.masked_where(T_matrix == 0, T_matrix)
	
	# Normalize each row
	for i in range(np.shape(T_matrix)[0]):
		T_masked[i,:] /= np.sum(T_masked[i,:])
		
	# Reshape to column vector
	T_matrix_flat = np.ndarray.flatten(T_masked)
	T_matrix_flat = np.transpose(T_matrix_flat)
	
	print np.shape(T_matrix_flat)
	
	return T_matrix_flat
	
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
	for idx, F_bin_location in enumerate(F_indices):
		qi[F_bin_location,:] += (qij[idx,:])
		qi_count[F_bin_location] += 1
	
	##DEBUGGING
	print qi_count
		
	#normalize the average value for the pair sum starting in state i
	for i in range(np.shape(qi)[0]):
		if qi_count[i] != [0]:
			qi[i,:] /= float(qi_count[i])
	
	##DEBUGGING
	print qi
	    
	Jacobian = np.zeros((nbins, npairs))
	for idx, t_bin_location in enumerate(t_indices):
		# Add q values for specific transition
		Jacobian[t_bin_location, :] += (qij[idx,:] + qij[idx+1,:])	
	
	# Normalize by i bin count
	for i in range(np.shape(Jacobian)[0]):
		state_i_idx = np.floor(i/nstates)
		if qi_count[state_i_idx] != 0:
			Jacobian[i,:] /= (2*qi_count[state_i_idx])
		
	for idx, t_bin_location in enumerate(t_indices):
		# Index for q value of all transitions starting at state i
		state_i_idx = np.floor(t_bin_location/nstates)
		
		# Subtract q values for all starting at state i
		Jacobian[t_bin_location, :] -= qi[state_i_idx,:]
	
#	Jacobian /= -GAS_CONSTANT_KJ_MOL*170
	
	return Jacobian
	
if __name__=="__main__":

	trace = np.loadtxt("FRETr.dat")
	print np.min(trace)
	print np.max(trace)
	
	T_vec = get_T_matrix(trace)
	print np.shape(T_vec)
	np.savetxt("T_vector.dat", T_vec)
	
	# Takes a really long time
	qij = np.loadtxt("qij.dat")
		
	hist, F_indices, t_indices, num_bins = find_FRET_bins(trace) # would normally go in compute_Jacobian_basic
	
	Jacobian = Jacobian_test(qij, hist, F_indices, t_indices)
	np.savetxt("Jacobian.dat", Jacobian)
	
	# PRINT STUFF	
#	print hist
	print np.shape(hist)
	
#	print F_indices
	print np.shape(F_indices)
	
#	print qij
	print np.shape(qij)
	
#	print t_indices
	print np.shape(t_indices)
	print np.max(t_indices)
	
#	print Jacobian
	print np.shape(Jacobian)
	
	# ELAPSED TIME TO COMPUTE JACOBIAN
#	print("--- %s seconds ---" % (time.time() - start_time))
	
