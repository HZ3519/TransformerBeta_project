import numpy as np

def search_target(seq_list, seq):
	"""Search whether one target are in the dataset or not."""

	return seq_list[seq_list[:, 0]==seq]



def search_target_similar(seq_list, seq, insimilarity_num = 0):
	"""Search whether similar targets are in the dataset or not."""

	similar_target_list = []

	for i in range(len(seq_list)):
		count = 0
		for j in range(len(seq_list[i, 0])):
			if seq_list[i, 0][j] != seq[j]:
				count += 1
		if count <= insimilarity_num:
			list_with_count = seq_list[i].tolist() + [count]
			similar_target_list.append(list_with_count)

	return np.array(similar_target_list)



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



def screen(validation_list, screen_list, window_size, return_target_condition = False, return_target_sequence = False, return_screen_condition = False, return_screen_sequence = False):

	if return_target_condition == True or return_target_sequence == True:
		target_condition_list = []
		count = 0

		for target in validation_list:
			sub_target = window_divide(target, window_size)
			target_condition = False

			for sub in sub_target:
				if sub in screen_list:
					target_condition = True
					break
			target_condition_list.append(target_condition)

			if count%10000 == 0:
				print(count)
			count += 1
		target_condition_array = np.array(target_condition_list)
		
		if return_target_condition:
			return target_condition_array
		if return_target_sequence:
			return validation_list[target_condition_array]
		else:
			raise ValueError("Please specify return_target_condition or return_target_sequence")

	if return_screen_condition == True or return_screen_sequence == True:
		screen_condition_list = []
		count = 0

		target_sequence_list = []
		for target in validation_list:
			sub_target = window_divide(target, window_size)
			for sub in sub_target:
				target_sequence_list += sub_target
		
		for screen in screen_list:
			if screen in target_sequence_list:
				screen_condition_list.append(True)
			else:
				screen_condition_list.append(False)
		
			if count%10000 == 0:
				print(count)
			count += 1
		screen_condition_array = np.array(screen_condition_list)

		if return_screen_condition:
			return screen_condition_array
		if return_screen_sequence:
			return screen_list[screen_condition_array]
		else:
			raise ValueError("Please specify return_screen_condition or return_screen_sequence")