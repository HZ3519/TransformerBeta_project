import torch
import numpy as np
import numba


def hamming_distance_training(seq_pred, seq_true):
	"""Return the Hamming distancefor training between 2 equal lists of equal-length sequences(list of number)."""
	
	hd_list = [(x!=y).sum() for x, y in zip(seq_pred, seq_true)]
	return torch.sum(torch.stack(hd_list))



def hamming_distance_sum(seq_pred, seq_true):
	"""Return the Hamming distance sum between 2 equal list of equal-length sequences(list of amino sequences)."""

	hd_list = [sum(np.array(list(x))!=np.array(list(y))) for x, y in zip(seq_pred, seq_true)]
	return sum(hd_list)



def hamming_distance_list(seq_pred, seq_true):
	"""Return the Hamming distance list between 2 equal list of equal-length sequences(list of amino sequences)."""

	hd_list = [sum(np.array(list(x))!=np.array(list(y))) for x, y in zip(seq_pred, seq_true)]
	return hd_list



@numba.jit(nopython=True)
def hamming_distance(s1, s2):
	"""
	Compute the Hamming distance between two sequences.
	"""
	distance = 0
	for i in range(len(s1)):
		if s1[i] != s2[i]:
			distance += 1
	return distance



@numba.jit(nopython=True)
def levenshtein_distance(s1, s2):
	"""
	Compute the Levenshtein distance between two sequences.
	"""

	if len(s1) < len(s2):
		return levenshtein_distance(s2, s1)
	if len(s2) == 0:
		return len(s1)

	previous_row = list(range(len(s2) + 1))
	for i, c1 in enumerate(s1):
		current_row = [i + 1]
		for j, c2 in enumerate(s2):
			insertions = previous_row[j + 1] + 1
			deletions = current_row[j] + 1
			substitutions = previous_row[j] + (c1 != c2)
			current_row.append(min(insertions, deletions, substitutions))
		previous_row = current_row
	return previous_row[-1]



@numba.jit(nopython=True, parallel=True)
def compute_MSA_weights_hamming(MSA, threshold = 0.1, verbose=False):
	"""
	Compute the MSA sequence weights as the inverse of the
	number of neighbouring sequences.
	A neighbour of a sequence is another sequence within a given
	threshold of Hamming distance (given as a fraction of the total
	sequence length).
	"""

	B = len(MSA)
	num_neighbours = np.zeros(B)
	for a in range(0, B):
		num_neighbours[a] += 1
		t = MSA[a]
		if (a % 10000 == 0) and verbose:
			print(a)
		num_neighbours_a_count = 0
		for b in numba.prange(a + 1, B):
			dist = hamming_distance(t, MSA[b])/len(t)
			if dist < threshold:
				num_neighbours_a_count += 1
				num_neighbours[b] += 1
		num_neighbours[a] += num_neighbours_a_count # avoid race condition
	return 1.0/num_neighbours



@numba.jit(nopython=True, parallel=True)
def compute_MSA_weights_Levenshtein(MSA, threshold = 0.1, verbose=False):
	"""
	Compute the MSA sequence weights as the inverse of the
	number of neighbouring sequences.
	A neighbour of a sequence is another sequence within a given
	threshold of Levenshtein distance (given as a fraction of the total
	sequence length).MSA string sequences do not need to be converted into
	number sequences before passing it to compute_MSA_weights.
	"""

	B = len(MSA)
	num_neighbours = np.zeros(B)
	for a in range(0, B):
		num_neighbours[a] += 1
		t = MSA[a]
		if (a % 10000 == 0) and verbose:
			print(a)
		num_neighbours_a_count = 0
		for b in numba.prange(a + 1, B):
			dist = levenshtein_distance(t, MSA[b])/max(len(t), len(MSA[b]))
			if dist < threshold:
				num_neighbours_a_count += 1
				num_neighbours[b] += 1
		num_neighbours[a] += num_neighbours_a_count # avoid race condition
	return 1.0/num_neighbours
