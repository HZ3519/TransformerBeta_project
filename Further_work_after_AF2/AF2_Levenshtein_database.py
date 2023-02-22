import numpy as np
import matplotlib.pyplot as plt
import time 
import os

from TransformerBeta import *

# 1. load the stored data
# length of interest >= 8
length_of_interest = []
type_of_interest = ['antiparallel'] # <------------------------------------------------------------------------change 
for i in range(8, 9): # <------------------------------------------------------------------------change 
	length_of_interest.append(str(i))
# length_of_interest.append('20more') # <------------------------------------------------------------------------change 

beta_strand_data = {}
# create subdictionary for antiparallel and parallel
beta_strand_data['antiparallel'] = {}
# beta_strand_data['parallel'] = {} # <------------------------------------------------------------------------change 
# create subdictionary of subdictionary for each length, string 3-20 and 20more
for i in range(8, 9): # <-------------------------------------------------------------------------change 
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

# 2. define the Levenshtein distance
def levenshtein_distance_AF2(s1, s2):
	"""
	Compute the Levenshtein distance between two sequences.
	"""

	if len(s1) < len(s2):
		return levenshtein_distance_AF2(s2, s1)

	# len(s1) >= len(s2)
	if len(s2) == 0:
		return len(s1)

	previous_row = list(range(len(s2) + 1))
	for i, c1 in enumerate(s1):
		current_row = [i + 1]
		for j, c2 in enumerate(s2):
			insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
			deletions = current_row[j] + 1       # than s2
			substitutions = previous_row[j] + (c1 != c2)
			current_row.append(min(insertions, deletions, substitutions))
		previous_row = current_row

	return previous_row[-1]

# 3. calculate similarity weighting at 90% identity using Levenshtein distance
num_test = [10, 100] 
for num in num_test:
	t1 = time.time()
	for type in beta_strand_data:
		for length in beta_strand_data[type]:
			list_same_length = []
			for key in beta_strand_data[type][length]:
				for comp in beta_strand_data[type][length][key]:
					list_same_length.append(key+comp)
			list_same_length = list_same_length[0:num] # <-------------------------------------------------------------- test
			# calculate the similarity weighting
			B = len(list_same_length)
			num_neighbours = np.zeros(B)
			for a in range(0, B):
				num_neighbours[a] += 1
				t = list_same_length[a]
				if (a % 10000 == 0):
					print(a)
				for b in range(a + 1, B):
					dist = levenshtein_distance_AF2(t, list_same_length[b])/max(len(t), len(list_same_length[b]))
					if dist <= 0.1:
						num_neighbours[a] += 1
						num_neighbours[b] += 1
			# store the similarity weighting
			for a in range(0, B):
				t = list_same_length[a]
				if length != '20more':
					count = beta_strand_data[type][length][t[0: int(length)]][t[int(length):]]
					beta_strand_data[type][length][t[0: int(length)]][t[int(length):]] = []
					beta_strand_data[type][length][t[0: int(length)]][t[int(length):]].append(count)
					beta_strand_data[type][length][t[0: int(length)]][t[int(length):]].append(num_neighbours[a])
				else:
					current_length = len(t//2)
					count = beta_strand_data[type]['20more'][t[0: current_length]][t[current_length:]]
					beta_strand_data[type]['20more'][t[0: current_length]][t[current_length:]] = []
					beta_strand_data[type]['20more'][t[0: current_length]][t[current_length:]].append(count)
					beta_strand_data[type]['20more'][t[0: current_length]][t[current_length:]].append(num_neighbours[a])
	t2 = time.time()
	print("Time taken for " + str(num) + " test for length " + length + " for type " + type + ": " + str(t2-t1) + "s")
	example_type = type_of_interest[0]
	example_length = length_of_interest[0]
	print("example of first 20 data: ", list(beta_strand_data[example_type][example_length].items())[0:20])



