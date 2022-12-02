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

# plot the distribution of "freq == 1" count as bar chart, the interval is at 1, 10, 100, 1000, 10000, 100000
count_distribution = np.zeros(6)
for i in range(len(AF_beta_strand_dataset)):
	if AF_beta_strand_dataset[i][3] == 1:
		count = AF_beta_strand_dataset[i][2]
		if count < 10:
			count_distribution[0] += 1
		elif count < 100:
			count_distribution[1] += 1
		elif count < 1000:
			count_distribution[2] += 1
		elif count < 10000:
			count_distribution[3] += 1
		elif count < 100000:
			count_distribution[4] += 1
		else:
			count_distribution[5] += 1
print("count distribution: " + str(count_distribution))

plt.figure(figsize=(24, 16), facecolor=(1, 1, 1))
plt.rcParams.update({'font.size': 25})
plt.bar(['1-9', '10-99', '100-999', '1000-9999', '10000-99999', '100000+'], count_distribution, color='blue')
plt.xlabel('Unique Count')
plt.ylabel('Number of sequences')
plt.title('Count distribution of sequences for unique anti-parallel beta strands with length 7')
plt.savefig('count_unique_distribution_l7_anti.png', dpi=300)
plt.close()
