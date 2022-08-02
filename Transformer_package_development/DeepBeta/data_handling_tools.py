def search_target(seq_list, seq):
	"""Search whether one target are in the dataset or not."""

	return seq_list[seq_list[:, 0]==seq]



def window_divide(sequence, windown_size):
	"""Divide the sequence into windows."""

	letter_list = list(sequence)

	sequence_list = []

	num_windows = len(letter_list) - windown_size + 1

	for i in range(num_windows):

		letter_list_i = letter_list[i:i+windown_size]
		sequence_i = "".join(letter_list_i)
		sequence_list.append(sequence_i)

	return sequence_list