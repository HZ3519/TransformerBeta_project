import numpy as np
from TransformerBeta import amino_dict
import os


# record start time
import time
start_time = time.time()

# load validation data
validation_dict = np.load('validation_l7_anti/validation_dict_fold0.npy', allow_pickle=True) # <------------------------------------------------------------------------change
validation_dict = validation_dict.tolist()
validation_list = []
for target, value_dict in validation_dict.items():
    for comp, count in value_dict.items():
        validation_list.append([target, comp, count])
validation_array = np.array(validation_list)

# load training data
training_dict = np.load('train_l7_anti/train_dict_fold0.npy', allow_pickle=True) # <------------------------------------------------------------------------change
training_dict = training_dict.tolist()
training_list = []
for target, value_dict in training_dict.items():
	for comp, count in value_dict.items():
		training_list.append([target, comp, count])
training_array = np.array(training_list)

# record time 1 for loading data
time_1 = time.time()
print('Time for loading data: ', time_1 - start_time)

# each row is a data point which is a list of [target, comp, count]
# target is a string of length 7, comp is a string of length 7, count is an integer
# combine train target and train comp into a list of length 14
training_combined = []
for i in range(len(training_array)):
	training_combined.append(training_array[i, 0] + training_array[i, 1])

# combine validation target and validation comp into a list of length 14
validation_combined = []
for i in range(len(validation_array)):
	validation_combined.append(validation_array[i, 0] + validation_array[i, 1])
training_combined = np.array(training_combined, dtype=object)
validation_combined = np.array(validation_combined, dtype=object)

# comvert the list of strings into arrays of integers
def seq2num(seq, amino_dict):
	seq_num = []
	for i in range(len(seq)):
		seq_num.append(amino_dict[seq[i]])
	return seq_num

training_combined_num = []
for i in range(len(training_combined)):
	training_combined_num.append(seq2num(training_combined[i], amino_dict))

validation_combined_num = []
for i in range(len(validation_combined)):
	validation_combined_num.append(seq2num(validation_combined[i], amino_dict))

# convert the list of arrays of integers into a 2D array
training_combined_num = np.array(training_combined_num)
validation_combined_num = np.array(validation_combined_num)

# record time 2 for converting data
time_2 = time.time()
print('Time for converting data: ', time_2 - time_1)
# print statistics of training data and validation data
print('training_combined_num.shape: ', training_combined_num.shape)
print('validation_combined_num.shape: ', validation_combined_num.shape)

validation_combined_num_split = np.array_split(validation_combined_num, 30)
# only keep 2 indices 10, 15 
validation_combined_num_split = [validation_combined_num_split[i] for i in [10, 15]]

print(len(validation_combined_num_split))

# HPC setting
num_resources = 2
counter = 0
array_index = int(os.environ['PBS_ARRAY_INDEX'])
array_index_list = list(range(1, num_resources+1))

# for each validation data point, find the closest hamming distance data point in the training set, record the distance
hamming_dist = []
for save_index, fold in enumerate(validation_combined_num_split):
	if array_index_list[counter%num_resources] == array_index:
		for i in range(len(fold)): # <------------------------------------------------------------------------change
			validation_data = fold[i] # <------------------------------------------------------------------------change
			hamming_dist_validation = np.min(np.sum(validation_data != training_combined_num, axis=1))
			hamming_dist.append(hamming_dist_validation)
		hamming_dist = np.array(hamming_dist)
		if save_index == 0:
			np.save('validation_l7_anti/closest_dist_fold0_' + str(10) + '.npy', hamming_dist) # <------------------------------------------------------------------------change
		elif save_index == 1:
			np.save('validation_l7_anti/closest_dist_fold0_' + str(15) + '.npy', hamming_dist)
	counter += 1

# record time 3 for finding hamming distance
time_3 = time.time()
print('Time for finding hamming distance: ', time_3 - time_2)