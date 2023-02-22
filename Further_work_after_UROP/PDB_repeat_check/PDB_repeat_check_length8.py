import pickle
import numpy as np
import matplotlib.pyplot as plt
import time 

from TransformerBeta import *

"""---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
"""Data preprocessing"""

amino_dict = {
		'<bos>': 0, 
		'<eos>': 1, 
		'<pad>': 2, 
		'<unk>': 3,
		'A': 4,
		'C': 5,
		'D': 6, 
		'E': 7,
		'F': 8, 
		'G': 9, 
		'H': 10,
		'I': 11, 
		'K': 12, 
		'L': 13, 
		'M': 14, 
		'N': 15, 
		'P': 16, 
		'Q': 17, 
		'R': 18, 
		'S': 19, 
		'T': 20, 
		'V': 21, 
		'W': 22, 
		'Y': 23 
		}

data_list_antiparallel= []
data_list_parallel= []

for i in range(1, 9):
	with open('BSn_libraries/BSn_libraries_copy/anti_frag_dic_{}.pkl'.format(i), 'rb') as f:
		data = pickle.load(f, encoding='latin1')
		data_list_antiparallel.append(data)

	with open('BSn_libraries/BSn_libraries_copy/para_frag_dic_{}.pkl'.format(i), 'rb') as f:
		data = pickle.load(f, encoding='latin1')
		data_list_parallel.append(data)

# target, complementary_seq, counts, promiscuity, length, working_score, hb_pattern, para/anti, freq
BSn_data = []
least_length = 3

for frag_i_data in data_list_parallel[least_length-1:]:
	for keys in frag_i_data.keys():

		length = len(keys)
		freq = len(frag_i_data[keys])
		for element in frag_i_data[keys]:

			working_score = length**2 * element.count_score - 0.01 * length * element.promiscuity_score
			list_i = [keys, element.complementary_sequence, element.count_score, element.promiscuity_score, length, working_score, element.hb_pattern, 0, freq]
			BSn_data.append(list_i)

for frag_i_data in data_list_antiparallel[least_length-1:]:
	for keys in frag_i_data.keys():

		length = len(keys)
		freq = len(frag_i_data[keys])
		for element in frag_i_data[keys]:

			working_score = length**2 * element.count_score - 0.01 * length * element.promiscuity_score
			list_i = [keys, element.complementary_sequence, element.count_score, element.promiscuity_score, length, working_score, element.hb_pattern, 1, freq]
			BSn_data.append(list_i)

# target, complementary_seq, counts, promiscuity, length, working_score, hb_pattern, para/anti, freq

BSn_data_dataset_sequence = np.array(BSn_data, dtype=object)
BSn_data_dataset1 = np.array(BSn_data_dataset_sequence[BSn_data_dataset_sequence[:, 5] >= 0])
BSn_data_dataset2 = np.array(BSn_data_dataset_sequence[BSn_data_dataset_sequence[:, 5] >= 0])


target_indices = np.arange(BSn_data_dataset2.shape[0]).reshape(-1, 1)
BSn_data_dataset2_indices = np.hstack([BSn_data_dataset2, target_indices])

condition1 = np.nonzero(np.array([len(sequence)==8 for sequence in BSn_data_dataset2_indices[:, 0]]))
BSn_data_dataset2_indices_length8 = BSn_data_dataset2_indices[condition1]

condition2 = np.nonzero(np.array([freq==1 for freq in BSn_data_dataset2_indices_length8[:, -2]]))
BSn_data_dataset2_indices_length8_unique = BSn_data_dataset2_indices_length8[condition2]

# set seed
np.random.seed(0)
validation_indices = np.random.choice(BSn_data_dataset2_indices_length8_unique[:, -1], size=5000, replace=False).astype(np.int32)

X_train_letter = np.delete(BSn_data_dataset2[:, 0], validation_indices, axis=0)
Y_train_letter = np.delete(BSn_data_dataset2[:, 1], validation_indices, axis=0)
X_validation_letter = BSn_data_dataset2[validation_indices, 0]
Y_validation_letter = BSn_data_dataset2[validation_indices, 1]

# keep length 8 training data
condition3 = np.nonzero(np.array([len(sequence)>=8 for sequence in X_train_letter]))
X_train_letter = X_train_letter[condition3]
Y_train_letter = Y_train_letter[condition3]


"""---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
"""repeat check"""
train_letter = [i+j for i, j in zip(X_train_letter, Y_train_letter)]
validation_letter = [i+j for i, j in zip(X_validation_letter, Y_validation_letter)]

print("train letter repeat check: ", len(train_letter) == len(set(train_letter)))
# pick out the repeated training sequence
# save the index of the repeated sequence
for i in train_letter:
	if train_letter.count(i) > 1:
		print(i)
		print(train_letter.index(i))
		print(train_letter.count(i))

print("validation letter repeat check: ", len(validation_letter) == len(set(validation_letter)))
# pick out the repeated validation sequence
# save the index of the repeated sequence
for i in validation_letter:
	if validation_letter.count(i) > 1:
		print(i)
		print(validation_letter.index(i))
		print(validation_letter.count(i))

print("train letter in validation letter check: ", len(set(train_letter).intersection(set(validation_letter))) == 0)
# pick out the training sequence hat is also in the validation set
for i in train_letter:
	if i in validation_letter:
		print(i)
		print(train_letter.index(i))
		print(validation_letter.index(i))

