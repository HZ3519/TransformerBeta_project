import pickle
import numpy as np
import torch
import torch.nn as nn
from d2l import torch as d2l
import os
import matplotlib.pyplot as plt

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

# 1. load the stored data
# length of interest >= 7
length_of_interest = []
type_of_interest = ['antiparallel'] # <------------------------------------------------------------------------change 
for i in range(7, 8): # <------------------------------------------------------------------------change 
	length_of_interest.append(str(i))
# length_of_interest.append('20more') # <------------------------------------------------------------------------change 

beta_strand_data = {}
# create subdictionary for antiparallel and parallel
beta_strand_data['antiparallel'] = {}
# beta_strand_data['parallel'] = {} # <------------------------------------------------------------------------change 
# create subdictionary of subdictionary for each length, string 3-20 and 20more
for i in range(7, 8): # <-------------------------------------------------------------------------change 
	beta_strand_data['antiparallel'][str(i)] = {}
	# beta_strand_data['parallel'][str(i)] = {} # <------------------------------------------------------------------------change 
# beta_strand_data['antiparallel']['20more'] = {} # <------------------------------------------------------------------------change 
# beta_strand_data['parallel']['20more'] = {} # <------------------------------------------------------------------------change 

total_data_stored = 0
for type in ['antiparallel', 'parallel']:
	if type not in type_of_interest:
		continue
	length_list = os.listdir("AF2_beta_strand_database/" + type)
	length_list.sort()
	for length in length_list:
		length_name = length.split('_')[1]
		if length_name not in length_of_interest:
			continue
		num_load = 0
		# get all subfolders
		length_subfolders = os.listdir("AF2_beta_strand_database/" + type + "/" + length)
		for subfolder in length_subfolders:
			# load .npy file
			data = np.load("AF2_beta_strand_database/" + type + "/" + length + "/" + subfolder, allow_pickle=True)
			data = data.tolist()
			# data is a dictionary
			# add data to beta_strand_data
			for key in data:
				beta_strand_data[type][length_name][key] = data[key]
				total_data_stored += len(data[key])
				num_load += len(data[key])
		# print number of data loaded in this length
		print('number of ' + type + ' beta strands with length ' + length_name + ' loaded: ' + str(num_load))
print("Total data stored: " + str(total_data_stored))

# 2. create the dataset
# append key, comp, count, freq to the dataset
AF_beta_strand_dataset = []
for type in beta_strand_data:
	for length in beta_strand_data[type]:
		for key in beta_strand_data[type][length]:
			freq = len(beta_strand_data[type][length][key])
			for comp in beta_strand_data[type][length][key]:
				AF_beta_strand_dataset.append([key, comp, beta_strand_data[type][length][key][comp], freq])
AF_beta_strand_dataset = np.array(AF_beta_strand_dataset, dtype=object)

dataset_indices = list(range(len(AF_beta_strand_dataset)))
dataset_indices = np.array(dataset_indices)

# filter training indices from the dataset: condition 1: count >= 5
condition_count = np.nonzero(np.array([count >=5 for count in AF_beta_strand_dataset[:, -2]]))
dataset_indices_high_count = dataset_indices[condition_count]
AF_beta_strand_dataset = AF_beta_strand_dataset[dataset_indices_high_count, :]
dataset_indices = list(range(len(AF_beta_strand_dataset)))
dataset_indices = np.array(dataset_indices)

# 4. build validation dataset of 10 fold of dataset_indices
# random shuffle the dataset_indices
np.random.seed(0)
np.random.shuffle(dataset_indices)
validation_indices_folds = []
for i in range(10):
	validation_indices_folds.append(dataset_indices[i*len(dataset_indices)//10:(i+1)*len(dataset_indices)//10])
# fold 0
validation_indices = validation_indices_folds[0]


X_train = np.delete(AF_beta_strand_dataset[:, 0], validation_indices, axis=0)
Y_train = np.delete(AF_beta_strand_dataset[:, 1], validation_indices, axis=0)
X_validation = AF_beta_strand_dataset[validation_indices, 0]
Y_validation = AF_beta_strand_dataset[validation_indices, 1]
train_count = np.delete(AF_beta_strand_dataset[:, 2], validation_indices, axis=0)
validation_count = AF_beta_strand_dataset[validation_indices, 2]

# combine X_train and Y_train into the dictionary
# combine X_validation and Y_validation into the dictionary
# record count as value
train_dict = {}
validation_dict = {}
for i in range(len(X_train)):
	if X_train[i] not in train_dict:
		train_dict[X_train[i]] = {}
	train_dict[X_train[i]][Y_train[i]] = train_count[i]
for i in range(len(X_validation)):
	if X_validation[i] not in validation_dict:
		validation_dict[X_validation[i]] = {}
	validation_dict[X_validation[i]][Y_validation[i]] = validation_count[i]


# save train_dict and validation_dict
if not os.path.exists("train_l7_anti"):
	os.mkdir("train_l7_anti")
if not os.path.exists("validation_l7_anti"):
	os.mkdir("validation_l7_anti")
np.save("train_l7_anti/train_dict_fold0.npy", train_dict)
np.save("validation_l7_anti/validation_dict_fold0.npy", validation_dict)

# print statistics
print("Number of training data (count>=5): " + str(X_train.shape[0]))
print("Number of validation data: " + str(X_validation.shape[0]))
