import torch
import numpy as np


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
