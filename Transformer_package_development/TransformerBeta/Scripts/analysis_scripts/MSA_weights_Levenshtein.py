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
"""Similarity weight histogram"""
import numba
# please spefify:
numba.set_num_threads(48)
threshold = 0.1

# compute Levenshtein distance
train_letter = [i+j for i, j in zip(X_train_letter, Y_train_letter)]
train_letter = numba.typed.List(train_letter)
MSA_weights_Levenshtein = compute_MSA_weights_Levenshtein(train_letter, threshold)

file_time = time.strftime('%y%b%d_%I%M%p', time.gmtime())

# save MSA weights
# save MSA weights histogram
plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
plt.rcParams.update({'font.size': 15})
plt.hist(MSA_weights_Levenshtein, bins=100, range=(0, 1), color='blue', edgecolor='black', linewidth=1.2)
plt.xlabel('MSA weights')
plt.ylabel('Frequency')
plt.title('MSA weights histogram_{}%Identity'.format((1-threshold)*100))
plt.grid()
plt.savefig('MSA_weights_Levenshtein_{}.png'.format(file_time))

np.savetxt('MSA_weights_Levenshtein_{}.csv'.format(file_time), MSA_weights_Levenshtein, delimiter=',', fmt='%s')	