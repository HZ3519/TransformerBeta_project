import pickle
import numpy as np
import torch
import torch.nn as nn
from d2l import torch as d2l
import os
from transformers import T5Tokenizer, T5Model
import re

from TransformerBeta import *

"""---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
"""Data preprocessing"""

length_want = '10'
type_want = "parallel"

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
type_of_interest = []
length_of_interest.append(length_want)
type_of_interest.append(type_want)

beta_strand_data = {}
beta_strand_data[type_want] = {}
beta_strand_data[type_want][length_want] = {}

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
data_total_len = len(AF_beta_strand_dataset)
dataset_indices = list(range(len(AF_beta_strand_dataset)))
# convert to numpy array
AF_beta_strand_dataset = np.array(AF_beta_strand_dataset, dtype=object)

# build a dataset of training that is similar to beta_strand_data so that comp and key are switched
beta_strand_data_reverse = {}
beta_strand_data_reverse[type_want] = {}
beta_strand_data_reverse[type_want][length_want] = {}

for data in AF_beta_strand_dataset:
	key = data[0]
	comp = data[1]
	count = data[2]
	length = str(len(key))
	if length not in length_of_interest:
		continue
	if comp not in beta_strand_data_reverse[type_want][length]:
		beta_strand_data_reverse[type_want][length][comp] = {}
		beta_strand_data_reverse[type_want][length][comp][key] = count
	else:
		beta_strand_data_reverse[type_want][length][comp][key] = count


"""---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
"""Calculate hamming weights"""

def hamming_weight_AF2(AF_beta_strand_dataset, beta_strand_data, beta_strand_data_reverse):
	hamming_weight_vector = np.ones(len(AF_beta_strand_dataset), dtype=int)

	# calculate hamming weights for AF_beta_strand_dataset
	for index, data in enumerate(AF_beta_strand_dataset):
		this_key = data[0]
		this_comp = data[1]
		other_comp = beta_strand_data[type_want][str(len(this_key))][this_key]
		for comp in other_comp:
			if comp == this_comp:
				continue
			# calculate hamming weight, if equals to 1
			hamming_weight = 0
			for i in range(len(this_key)):
				if this_comp[i] != comp[i]:
					hamming_weight += 1
			if hamming_weight == 1:
				hamming_weight_vector[index] += 1
		other_key = beta_strand_data_reverse[type_want][str(len(this_key))][this_comp]
		for key in other_key:
			if key == this_key:
				continue
			# calculate hamming weight, if equals to 1
			hamming_weight = 0
			for i in range(len(this_key)):
				if this_key[i] != key[i]:
					hamming_weight += 1
			if hamming_weight == 1:
				hamming_weight_vector[index] += 1

	# take inverse of hamming weight vector
	hamming_weight_vector = 1 / hamming_weight_vector

	return hamming_weight_vector

import time
start_time = time.time()

AF_beta_strand_dataset = AF_beta_strand_dataset
hamming_weight_vector = hamming_weight_AF2(AF_beta_strand_dataset, beta_strand_data, beta_strand_data_reverse)
print("data runtime: " + str(time.time() - start_time))

# get distribution of hamming weights
hamming_weight_distribution = {}
for weight in hamming_weight_vector:
	if weight not in hamming_weight_distribution:
		hamming_weight_distribution[weight] = 1
	else:
		hamming_weight_distribution[weight] += 1

# put the dict in order of key from largest to smallest
hamming_weight_distribution = dict(sorted(hamming_weight_distribution.items(), key=lambda item: item[0], reverse=True))

# save hamming_weight_vector as one .npy file
print("hamming_weight_vector" + "_" + type_want + "_" + length_want + ": ", hamming_weight_vector)
np.save("AF2_beta_strand_hamming_weight_vector" + "_" + type_want + "_" + length_want, hamming_weight_vector)

# print and save hamming_weight_distribution
print("hamming_weight_distribution" + "_" + type_want + "_" + length_want + ": ", hamming_weight_distribution)
with open("AF2_beta_strand_hamming_weight_distribution" + "_" + type_want + "_" + length_want + ".txt", "w") as f:
	# write total data 
	f.write("total data " + type_want + " " + length_want + ": " + str(data_total_len) + "\n")
	for key in hamming_weight_distribution:
		f.write(str(key) + " " + str(hamming_weight_distribution[key]) + "\n")
