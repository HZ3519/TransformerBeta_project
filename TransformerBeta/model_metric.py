import torch
import numpy as np


def hamming_distance_training(seq_pred, seq_true):
    """Return the Hamming distancefor training between 2 equal lists of equal-length sequences(list of number)."""

    hd_list = [(x != y).sum() for x, y in zip(seq_pred, seq_true)]
    return torch.sum(torch.stack(hd_list))


def hamming_distance_sum(seq_pred, seq_true):
    """Return the Hamming distance sum between 2 equal list of equal-length sequences(list of amino sequences)."""

    hd_list = [
        sum(np.array(list(x)) != np.array(list(y))) for x, y in zip(seq_pred, seq_true)
    ]
    return sum(hd_list)


def hamming_distance_list(seq_pred, seq_true):
    """Return the Hamming distance list between 2 equal list of equal-length sequences(list of amino sequences)."""

    hd_list = [
        sum(np.array(list(x)) != np.array(list(y))) for x, y in zip(seq_pred, seq_true)
    ]
    return hd_list


def compute_MSA_weights(MSA, threshold=0.1, verbose=False):
    """
    Compute the MSA sequence weights as the inverse of the
    number of neighbouring sequences.
    A neighbour of a sequence is another sequence within a given
    threshold of Hamming distance (given as a fraction of the total
    sequence length). MSA string sequences must be converted into
    number sequences with seq2num function before passing it to
    compute_MSA_weights.
    """

    B = MSA.shape[0]
    N = MSA.shape[1]
    num_neighbours = np.zeros(B)
    for a in range(0, B):
        num_neighbours[a] += 1
        t = MSA[a]
        if (a % 10000 == 0) and verbose:
            print(a)
        for b in range(a + 1, B):
            dist = np.mean(t != MSA[b])
            if dist < threshold:
                num_neighbours[a] += 1
                num_neighbours[b] += 1
    return 1.0 / num_neighbours
