import pickle
import numpy as np
import torch
import torch.nn as nn
from d2l import torch as d2l
import os

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
type_of_interest = ['antiparallel', 'parallel'] # <------------------------------------------------------------------------change 
for i in range(7, 8): # <------------------------------------------------------------------------change 
	length_of_interest.append(str(i))
# length_of_interest.append('20more') # <------------------------------------------------------------------------change 

beta_strand_data = {}
# create subdictionary for antiparallel and parallel
beta_strand_data['antiparallel'] = {}
beta_strand_data['parallel'] = {} # <------------------------------------------------------------------------change 
# create subdictionary of subdictionary for each length, string 3-20 and 20more
for i in range(7, 8): # <-------------------------------------------------------------------------change 
	beta_strand_data['antiparallel'][str(i)] = {}
	beta_strand_data['parallel'][str(i)] = {} # <------------------------------------------------------------------------change 
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

# keep the training data with top 5 million count
AF_beta_strand_dataset = AF_beta_strand_dataset[AF_beta_strand_dataset[:, 2].argsort()[::-1]]
AF_beta_strand_dataset = AF_beta_strand_dataset[:2000000, :]

dataset_indices = list(range(len(AF_beta_strand_dataset)))
dataset_indices = np.array(dataset_indices)

# filter validation indices from the dataset
condition_freq = np.nonzero(np.array([freq==1 for freq in AF_beta_strand_dataset[:, -1]]))
dataset_indices_unique = dataset_indices[condition_freq]

# set seed
np.random.seed(0)
validation_indices = np.random.choice(dataset_indices_unique, size=10000, replace=False).astype(np.int32)

X_train = np.delete(AF_beta_strand_dataset[:, 0], validation_indices, axis=0)
Y_train = np.delete(AF_beta_strand_dataset[:, 1], validation_indices, axis=0)
X_validation = AF_beta_strand_dataset[validation_indices, 0]
Y_validation = AF_beta_strand_dataset[validation_indices, 1]

# split the data in training and validation
num_steps_training = 9 # <------------------------------------------------------------------------change 

X_train, X_valid_len, Y_train, Y_valid_len, X_validation, X_validation_valid_len, Y_validation, Y_validation_valid_len = preprocess_train(X_train, Y_train, amino_dict, num_steps_training, X_validation_letter=X_validation, Y_validation_letter=Y_validation)

working_score_tensor = torch.ones(X_train.shape[0], dtype=torch.float32) # equal weight for all training data

# print statistics
print("Number of training data (after filtering top 2M): " + str(X_train.shape[0]))
print("Number of unique training data (after filtering top 2M): ", len(dataset_indices_unique))
print("Number of validation data (after filtering top 2M): " + str(X_validation.shape[0]))





"""---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
""""Model Training"""
# please specify:
# 1. training_steps: 300k
# 2. model_name:
# 3. warmup_steps: 8k

model_name = 'AF2_transformer_standard_validation_58M_DP_2Mdata'
if not os.path.exists(model_name):
	os.makedirs(model_name)
with open(model_name + '/print_message.txt', 'w') as f:
	f.write("Number of training data (after filtering top 2M): " + str(X_train.shape[0]) + "\n")
	f.write("Number of unique training data (after filtering top 2M): " + str(len(dataset_indices_unique)) + "\n")
	f.write("Number of validation data (after filtering top 2M): " + str(X_validation.shape[0]) + "\n")

query_size, key_size, value_size, num_hiddens = 512, 512, 512, 512
num_layers, dropout = 8, 0.1
lr, training_steps, batch_size, label_smoothing = 0.0004, 300000, 4096, 0.1
ffn_num_input, ffn_num_hiddens, num_heads = 512, 2048, 8

norm_shape = [512] # 512 corresponds to the dim of such number to normalize
device = d2l.try_gpu()

encoder_standard = TransformerEncoder(
	len(amino_dict), key_size, query_size, value_size, num_hiddens, 
	norm_shape, ffn_num_input, ffn_num_hiddens, num_heads,
	num_layers, dropout)
decoder_standard = TransformerDecoder(
	len(amino_dict), key_size, query_size, value_size, num_hiddens, 
	norm_shape, ffn_num_input, ffn_num_hiddens, num_heads,
	num_layers, dropout, shared_embedding=encoder_standard.embedding)
model_standard = EncoderDecoder(encoder_standard, decoder_standard)


model_standard_total_params = sum(p.numel() for p in model_standard.parameters())
model_standard_total_trainable_params = sum(p.numel() for p in model_standard.parameters() if p.requires_grad)

print('Standard model: total number of parameters: {}'.format(model_standard_total_params))
print('Standard model: total number of trainable parameters: {}'.format(model_standard_total_trainable_params))

with open(model_name + '/print_message.txt', 'a') as f:
	f.write('Standard model: total number of parameters: {}'.format(model_standard_total_params) + "\n")
	f.write('Standard model: total number of trainable parameters: {}'.format(model_standard_total_trainable_params) + "\n")


optimizer = torch.optim.Adam(model_standard.parameters(), lr=lr, betas=(0.9, 0.98), eps = 1.0e-9)
warmup = 8000
scheduler = WarmupCosineSchedule(optimizer, warmup, t_total=training_steps)

if torch.cuda.device_count() > 1:
  print("Let's use", torch.cuda.device_count(), "GPUs!")
  model_standard = nn.DataParallel(model_standard)
  with open(model_name + '/print_message.txt', 'a') as f:
  	f.write("Let's use " + str(torch.cuda.device_count()) + " GPUs!" + "\n")

train_seq2seq_training_steps_DP(model_standard, X_train, X_valid_len, Y_train, Y_valid_len, working_score_tensor, lr, training_steps, batch_size, label_smoothing, amino_dict, device, model_name=model_name, warmup=scheduler, optimizer=optimizer, X_validation=X_validation, Y_validation=Y_validation, X_validation_valid_len=X_validation_valid_len, Y_validation_valid_len=Y_validation_valid_len)