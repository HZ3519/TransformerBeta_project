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

# reverse screen
# please specify:
window_size_list = [3, 4, 5, 6, 7] # fragments length for reversing screening in training 
length = 7  # fragment length above to use for training 

target_indices = np.arange(BSn_data_dataset2.shape[0]).reshape(-1, 1)
BSn_data_dataset2_indices = np.hstack([BSn_data_dataset2, target_indices])

condition1 = np.nonzero(np.array([length==8 for length in BSn_data_dataset2_indices[:, 4]]))
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

sequence_remove_list = []
X_train_letter_screening_sum= 0 

for window_size in window_size_list:

	condition3 = np.nonzero(np.array([len(sequence)==window_size for sequence in X_train_letter]))
	X_train_letter_length_i = X_train_letter[condition3]
	X_train_letter_screening_sum += X_train_letter_length_i.shape[0]
	print("number of length {} samples for screening: {}".format(window_size, X_train_letter_length_i.shape))
	screen_sequence_to_remove = screen(X_validation_letter, X_train_letter_length_i, window_size, return_screen_sequence=True)
	sequence_remove_list.extend(screen_sequence_to_remove)

sequence_remove_dict = {}
for sequence in sequence_remove_list:
	sequence_remove_dict[sequence] = 1

condition_list = []
for sequence in X_train_letter:
	try:
		sequence_remove_dict[sequence]
		condition_list.append(False)
	except:
		condition_list.append(True)

X_train_letter = X_train_letter[np.array(condition_list)]
Y_train_letter = Y_train_letter[np.array(condition_list)]

condition3 = np.nonzero(np.array([len(sequence)>=length for sequence in X_train_letter]))
X_train_letter = X_train_letter[condition3]
Y_train_letter = Y_train_letter[condition3]

print("------------validation set info---------")
print("number of length 8 samples: ", BSn_data_dataset2_indices_length8.shape)
print("number of unique length 8 sample: ", BSn_data_dataset2_indices_length8_unique.shape)
print("number of validation samples: ", X_validation_letter.shape)
print("------------training set info---------")
print("number of all training samples going though reverse screening", X_train_letter_screening_sum)
print("number of removed traininng sequences: ", len(sequence_remove_list))
print("number of all training samples after reverse screening: ", len(X_train_letter))
print("Above may not be consistent if screen and length used is not the same")
print("------------ data below---------")
print(X_train_letter.shape)
print(Y_train_letter.shape)
print(X_validation_letter.shape)
print(Y_validation_letter.shape)



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
plt.title('MSA weights histogram_{}%Identity_length78'.format((1-threshold)*100))
plt.grid()
plt.savefig('MSA_weights_Levenshtein_{}.png'.format(file_time))

np.savetxt('MSA_weights_Levenshtein_{}.csv'.format(file_time), MSA_weights_Levenshtein, delimiter=',', fmt='%s')	