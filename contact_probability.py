import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

class solution_container(object):
    def __init__(self, contacts_list, contact_probability_list, distance_pairs, hamming_matrix, factor_distances, shift_distances):
        self.contacts_list = contacts_list
        self.contact_probability_list = contact_probability_list
        self.distance_pairs = distance_pairs
        self.hamming_matrix = hamming_matrix
        self.factor_distances = factor_distances #distances as a factor of r0
        self.shift_distances = shift_distances #distances as a difference of r0


def plot_probability(pairs, n_residues, save_file, native_pairs=None, probability=None, native_probability=None):
    """ Plot the contact probability in a grid

    Native pairs in the lower-right half.

    Parameters
    ----------
    pairs : numpy.ndarray
        A Nx2 array of pairs. If its a list, will attempt to turn into array.
    n_residues : int
        Number of total residues in the protein system.
    save_file : str
        Location to save plot.
    native_pairs : numpy.ndarray, Default = None
        A Mx2 array of pairs. If its a list, will attempt to turn into array.
        If None, native_pairs is equal to pairs.
    probability : list-like, Default = None
        List of size N. Probability of contact between residues specified in
        pairs. If None, defaults to 1 for all contacts.
    native_probability : list-like
        List of size M. Probability of contact between residues specified in
        native_pairs. If None, defaults to 1 for all contacts.
    """

    if isinstance(pairs, list):
        pairs = np.array(pairs)
    if isinstance(native_pairs, list):
        native_pairs = np.array(native_pairs)

    if native_pairs is None:
        native_pairs = pairs
    n_pairs = np.shape(pairs)[0]
    n_native_pairs = np.shape(native_pairs)[0]

    if probability is None:
        probability = np.ones(n_pairs)
    if native_probability is None:
        native_probability = np.ones(n_native_pairs)

    assert np.shape(probability)[0] == n_pairs
    assert np.shape(native_probability)[0] == n_native_pairs

    contact_matrix = np.zeros((n_residues,n_residues))
    for idx in range(n_native_pairs):
        i = native_pairs[idx, 0]
        j = native_pairs[idx, 1]
        contact_matrix[i,j] = native_probability[idx]
    for idx in range(n_pairs):
        i = pairs[idx, 0]
        j = pairs[idx, 1]
        if contact_matrix[i,j] != 0:
            contact_matrix[j,i] = probability[idx]
        else:
            contact_matrix[j,i] = -probability[idx]

    edges = np.arange(0.5, n_residues+1, 1)
    plt.figure()
    qmesh = plt.pcolormesh(edges, edges, contact_matrix, vmin=-1, vmax=1, cmap="bwr_r")
    plt.plot([0, n_residues + 1], [0, n_residues + 1], color="k", linewidth=2)
    cb = plt.colorbar(qmesh)
    cb.set_label("contact", size=16)
    plt.xlabel("i", size=18)
    plt.ylabel("j", size=18)
    plt.axis([0.5, n_residues+0.5, 0.5, n_residues+0.5])
    plt.savefig(save_file)

def compare_contacts(pairs_1, pairs_2, probability_1=None, probability_2=None, hamming=False, no_abs=False):
    """ Compare two sets of contacts and returns a distance_matrix.
    Euclidean Distance : hamming = False
        Computes a euclidean distance between the two list of contacts, by
        taking the difference in their respective probabilities.
    Hamming distance : hamming = True
        Computes a hamming distance between the two list of contacts.
        Intuitively, this is the number of contacts present in only one of the
        contacts list but not the other.

    """
    try:
        max1 = np.max(pairs_1)
        max2 = np.max(pairs_2)
        num_make = int(np.max([max1, max2]) + 1)
    except:
        np.savetxt("DEBUG_pair_1.dat", pairs_1)
        np.savetxt("DEBUG_pair_2.dat", pairs_2)
        raise IOError("Something is wrong with the inputted pairs. See DEBUG_pair_1.dat and DEBUG_pair_2.dat")

    if not hamming:
        if probability_1 is None:
            probability_1 = np.ones()
        if probability_2 is None:
            probability_2 = np.ones()

    compare_matrix = np.zeros((num_make,num_make))
    for idx in range(np.shape(pairs_1)[0]): #for pairs_1
        i = pairs_1[idx,0]
        j = pairs_1[idx,1]
        assert i < j
        if hamming:
            compare_matrix[i,j] = 1
        else:
            compare_matrix[i,j] = probability_1[idx]

    for idx in range(np.shape(pairs_2)[0]): #for pairs_2
        try:
            i = pairs_2[idx,0]
            j = pairs_2[idx,1]
        except:
            print "ERROR"
            print pairs_2
            print type(pairs_2)
            print idx
            raise
        assert i < j
        if hamming:
            compare_matrix[j,i] = 1
        else:
            compare_matrix[j,i] = probability_2[idx]

    total = 0. #treat like euclidean distance
    for i in range(num_make):
        for j in range(i+1, num_make):
            diff =  compare_matrix[j,i] - compare_matrix[i,j]
            if no_abs:
                total += diff
            else:
                total += diff ** 2

    if not hamming:
        total = np.sqrt(total)

    return total


def get_traj_sliced(traj, hmm_memberships, state_idx, dtraj_all, state_cutoff=0.95):
    """ Return trajectory object sliced for specific states"""
    selected_states = np.where(hmm_memberships[:,state_idx] > state_cutoff)[0]
    for count,index in enumerate(selected_states):
        if count == 0:
            rows = np.where(dtraj_all == index)[0]
        else:
            rows = np.append(rows, np.where(dtraj_all == index)[0])

    return traj[rows]

def get_rows_sliced(hmm_memberships, state_idx, dtraj_all, state_cutoff=0.95):
    """ Return trajectory object sliced for specific states"""
    selected_states = np.where(hmm_memberships[:,state_idx] > state_cutoff)[0]
    for count,index in enumerate(selected_states):
        if count == 0:
            rows = np.where(dtraj_all == index)[0]
        else:
            rows = np.append(rows, np.where(dtraj_all == index)[0])

    return rows

def compute_aa_contacts_from_array(selected_distances, distance_pairs, distance_cutoff=0.35, contact_probability_cutoff=0.3):
    """ Compute Contacts by finding closest heavy atoms within a cutoff """

    contacts_list = []
    contact_probability_list = []

    num_distances = np.shape(distance_pairs)[0]
    num_rows = np.shape(selected_distances)[0]
    hamming_matrix = np.zeros((num_rows, num_distances)) # 1 when in contact
    for distance_index in range(num_distances):
        one_d_hamming = np.zeros(num_rows)
        contact_indices = np.where(selected_distances[:,distance_index]<=distance_cutoff)
        one_d_hamming[contact_indices] = 1
        hamming_matrix[:,distance_index] = one_d_hamming
        in_contact = contact_indices[0]
        num_in_contact = np.shape(in_contact)[0]
        contact_probability = float(num_in_contact) / float(num_rows)
        if contact_probability > contact_probability_cutoff:
            #print "contact formed between %d-%d, Probability: %f" % (distance_pairs[distance_index,0], distance_pairs[distance_index,1], contact_probability)
            contacts_list.append(distance_pairs[distance_index])
            contact_probability_list.append(contact_probability)

    return contacts_list, contact_probability_list, distance_pairs, hamming_matrix


def compute_aa_contacts(traj, distance_cutoff=0.35, contact_probability_cutoff=0.3):
    """ Compute Contacts by finding closest heavy atoms within a cutoff """

    contacts_list = []
    contact_probability_list = []

    selected_distances, distance_pairs = md.compute_contacts(traj)
    num_distances = np.shape(distance_pairs)[0]
    num_rows = np.shape(selected_distances)[0]
    hamming_matrix = np.zeros((num_rows, num_distances)) # 1 when in contact
    for distance_index in range(num_distances):
        one_d_hamming = np.zeros(num_rows)
        contact_indices = np.where(selected_distances[:,distance_index]<=distance_cutoff)
        one_d_hamming[contact_indices] = 1
        hamming_matrix[:,distance_index] = one_d_hamming
        in_contact = contact_indices[0]
        num_in_contact = np.shape(in_contact)[0]
        contact_probability = float(num_in_contact) / float(num_rows)
        if contact_probability > contact_probability_cutoff:
            #print "contact formed between %d-%d, Probability: %f" % (distance_pairs[distance_index,0], distance_pairs[distance_index,1], contact_probability)
            contacts_list.append(distance_pairs[distance_index])
            contact_probability_list.append(contact_probability)


    return solution_container(contacts_list, contact_probability_list, distance_pairs, hamming_matrix, None, None)

def get_cg_contact_params(data_file="contacts_by_residue_ww_domain.dat"):
    f = open(data_file, "r")

    residue_pairs = []
    atom_pairs = []
    r0 = []
    for line in f:
        if line[0] == "#":
            pass
        else:
            data = line.strip().split()
            residue_pairs.append([int(data[0]), int(data[1])])
            atom_pairs.append([int(data[2]), int(data[3])])
            r0.append(float(data[4]))

    atom_pairs = np.array(atom_pairs)
    r0 = np.array(r0)

    max_count = len(residue_pairs)
    f.close()

    return residue_pairs, atom_pairs, r0

def run_cg_computations(traj, compute_pairs, r0, r_original, num_total, hamming_matrix=None, shift_distances=None, factor_distances=None):
    md_distances = md.compute_distances(traj, compute_pairs, periodic=False)
    num_total = np.shape(md_distances)[0]
    data_shape = np.shape(md_distances)
    distances = np.zeros(data_shape)
    shift_dist = np.zeros(data_shape)
    factor_dist = np.zeros(data_shape)
    for idx in range(num_total):
        distances[idx,:] = md_distances[idx,:] - r0
        shift_dist[idx,:] = md_distances[idx,:] - r_original
        factor_dist[idx,:] = md_distances[idx,:] / r_original

    this_factor_distance = np.amin(factor_dist, axis=1)
    this_shift_distance = np.amin(shift_dist, axis=1)
    assert np.shape(this_factor_distance)[0] == num_total
    assert np.shape(this_shift_distance)[0] == num_total

    contact = np.where(distances <= 0)
    not_contact = np.where(distances > 0)
    distances[not_contact] = 0
    distances[contact] = 1
    contact_distances = np.sum(distances, axis=1)
    assert np.shape(contact_distances)[0] == num_total
    contact_indices = np.where(contact_distances > 0)
    #compute probability
    num_in_contact = np.shape(contact_indices[0])[0]
    contact_probability = float(num_in_contact) / float(num_total)
    # compute and append hamming matrix
    one_d_hamming = np.zeros(num_total)
    one_d_hamming[contact_indices] = 1

    if hamming_matrix is None:
        hamming_matrix = [one_d_hamming]
    else:
        hamming_matrix.append(one_d_hamming)

    if shift_distances is None:
        shift_distances = [this_shift_distance]
    else:
        shift_distances.append(this_shift_distance)

    if factor_distances is None:
        factor_distances = [this_factor_distance]
    else:
        factor_distances.append(this_factor_distance)


    return hamming_matrix, factor_distances, shift_distances, contact_probability

def compute_cg_contacts(traj, residue_pairs, atom_pairs, r0, factor=None, shift=None, probability_cutoff=0.3):
    max_count = len(residue_pairs)
    r_start = np.copy(r0)
    if not factor is None:
        r0 = r0*factor
    if not shift is None:
        r0 = r0 + shift

    contacts_list = []
    contact_probability_list = []

    go = True
    count = 0
    current_residue_pair = residue_pairs[count]
    compute_idxs = []

    distance_pairs = []
    hamming_matrix = None
    shift_distances = None
    factor_distances = None
    num_total = traj.n_frames
    while go:
        if current_residue_pair[0] ==  residue_pairs[count][0] and current_residue_pair[1] == residue_pairs[count][1]:
            # get all the atoms for a residue-interaction
            compute_idxs.append(count)
        else:
            # once you have all the residue stuff, compute the previous residue stuff
            assert (len(compute_idxs)==2) or (len(compute_idxs)==4) or (len(compute_idxs)==1)

            distance_pairs.append(current_residue_pair)

            compute_pairs = atom_pairs[compute_idxs, :]
            r0_comparison = r0[compute_idxs]
            r_original = r_start[compute_idxs]

            hamming_matrix, factor_distances, shift_distances, contact_probability = run_cg_computations(traj, compute_pairs, r0_comparison, r_original, num_total, hamming_matrix=hamming_matrix, factor_distances=factor_distances, shift_distances=shift_distances)

            if contact_probability > probability_cutoff:
                contacts_list.append(current_residue_pair)
                contact_probability_list.append(contact_probability)

            current_residue_pair = residue_pairs[count]
            compute_idxs = [count]
        count += 1
        if count >= max_count:
            go = False

    assert (len(compute_idxs)==2) or (len(compute_idxs)==4) or (len(compute_idxs)==1)

    distance_pairs.append(current_residue_pair)

    compute_pairs = atom_pairs[compute_idxs, :]
    r0_comparison = r0[compute_idxs]
    r_original = r_start[compute_idxs]

    hamming_matrix, factor_distances, shift_distances, contact_probability = run_cg_computations(traj, compute_pairs, r0_comparison, r_original, num_total, hamming_matrix=hamming_matrix, factor_distances=factor_distances, shift_distances=shift_distances)

    if contact_probability > probability_cutoff:
        contacts_list.append(current_residue_pair)
        contact_probability_list.append(contact_probability)

    hamming_matrix = np.array(hamming_matrix).transpose()
    factor_distances = np.array(factor_distances).transpose()
    shift_distances = np.array(shift_distances).transpose()

    try:
        assert np.shape(hamming_matrix)[0] == num_total
    except:
        print np.shape(hamming_matrix)
        print num_total
        raise

    return solution_container(contacts_list, contact_probability_list, distance_pairs, hamming_matrix, factor_distances, shift_distances)
