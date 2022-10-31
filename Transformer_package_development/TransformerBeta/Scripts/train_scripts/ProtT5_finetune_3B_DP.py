import pickle
import numpy as np
import torch
import torch.nn as nn
from d2l import torch as d2l
import os
from transformers import T5Tokenizer, T5Model
import re
import torch_optimizer as optim

from TransformerBeta import *

# set environment variable set 'PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:512'
# to avoid CUDA out of memory error
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:512'

"""---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
"""Data preprocessing"""

tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_bfd', do_lower_case=False)
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
#amino_T5_dict = {}
#for key, value in amino_dict.items():
#	amino_T5_dict[key] = tokenizer.convert_tokens_to_ids(key)
#amino_T5_dict['<eos>'] = 1
#amino_T5_dict['<bos>'] = 1000
#amino_T5_dict['<unk>'] = 2

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
#beta_strand_data['parallel'] = {} # <------------------------------------------------------------------------change 
# create subdictionary of subdictionary for each length, string 3-20 and 20more
for i in range(7, 8): # <-------------------------------------------------------------------------change 
	beta_strand_data['antiparallel'][str(i)] = {}
	#beta_strand_data['parallel'][str(i)] = {} # <------------------------------------------------------------------------change 
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
dataset_indices = list(range(len(AF_beta_strand_dataset)))
# convert to numpy array
AF_beta_strand_dataset = np.array(AF_beta_strand_dataset, dtype=object)
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

X_train, X_valid_len, Y_train, Y_valid_len, X_validation, X_validation_valid_len, Y_validation, Y_validation_valid_len = preprocess_train_T5(X_train, Y_train, amino_dict, num_steps_training, X_validation_letter=X_validation, Y_validation_letter=Y_validation)

working_score_tensor = torch.ones(X_train.shape[0], dtype=torch.float32) # equal weight for all training data

# print statistics
print("Number of training data: " + str(X_train.shape[0]))
print("Number of unique training data: ", len(dataset_indices_unique))
print("Number of validation data: " + str(X_validation.shape[0]))





"""---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
""""Model Training"""
# please specify:
# 1. training_steps: 200k
# 2. model_name:

model_name = 'ProtT5_finetune_3B_DP'
if not os.path.exists(model_name):
	os.makedirs(model_name)
with open(model_name + '/print_message.txt', 'w') as f:
	f.write("Number of training data: " + str(X_train.shape[0]) + "\n")
	f.write("Number of unique training data: " + str(len(dataset_indices_unique)) + "\n")
	f.write("Number of validation data: " + str(X_validation.shape[0]) + "\n")

lr, training_steps, batch_size, label_smoothing = 0.0004, 10000, 4, 0.1
device = d2l.try_gpu()

ProtT5 = T5Model.from_pretrained("Rostlab/prot_t5_xl_bfd")
ProtT5_finetune_model = ProtT5_finetune(len(amino_dict), ProtT5)

ProtT5_finetune_model_total_params = sum(p.numel() for p in ProtT5_finetune_model.parameters())
ProtT5_finetune_model_total_trainable_params = sum(p.numel() for p in ProtT5_finetune_model.parameters() if p.requires_grad)

print('ProtT5_finetune_model: total number of parameters: {}'.format(ProtT5_finetune_model_total_params))
print('ProtT5_finetune_model: total number of trainable parameters: {}'.format(ProtT5_finetune_model_total_trainable_params))

with open(model_name + '/print_message.txt', 'a') as f:
	f.write('ProtT5_finetune_model: total number of parameters: {}'.format(ProtT5_finetune_model_total_params) + "\n")
	f.write('ProtT5_finetune_model: total number of trainable parameters: {}'.format(ProtT5_finetune_model_total_trainable_params) + "\n")


# optimizer = torch.optim.Adam(ProtT5_finetune_model.parameters(), lr=lr, betas=(0.9, 0.98), eps = 1.0e-9)
optimizer = optim.Adafactor(ProtT5_finetune_model.parameters(), lr=lr)
warmup = 4000
scheduler = WarmupCosineSchedule(optimizer, warmup, t_total=training_steps)

if torch.cuda.device_count() > 1:
  print("Let's use", torch.cuda.device_count(), "GPUs!")
  ProtT5_finetune_model = nn.DataParallel(ProtT5_finetune_model)
  with open(model_name + '/print_message.txt', 'a') as f:
  	f.write("Let's use " + str(torch.cuda.device_count()) + " GPUs!" + "\n")

train_prott5_DP(ProtT5_finetune_model, X_train, X_valid_len, Y_train, Y_valid_len, working_score_tensor, lr, training_steps, batch_size, label_smoothing, amino_dict, device, model_name=model_name, warmup=scheduler, optimizer=optimizer, X_validation=X_validation, Y_validation=Y_validation, X_validation_valid_len=X_validation_valid_len, Y_validation_valid_len=Y_validation_valid_len)