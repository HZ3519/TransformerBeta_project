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


train_dict = np.load('train_l7_anti/train_dict_fold0.npy', allow_pickle=True)
train_dict = train_dict.tolist()
train_list = []
for target, value_dict in train_dict.items():
    for comp, count in value_dict.items():
        train_list.append([target, comp, count])
train_array = np.array(train_list)

X_train = train_array[:, 0]
Y_train = train_array[:, 1]
# count_array = train_array[:, 2]


# 5. split the data in training and validation
num_steps_training = 9 # <------------------------------------------------------------------------change 

X_train, X_valid_len, Y_train, Y_valid_len = preprocess_train(X_train, Y_train, amino_dict, num_steps_training)

working_score_tensor = torch.ones(X_train.shape[0], dtype=torch.float32) # equal weight for all training data

# print statistics
print("Number of training data: " + str(X_train.shape[0]))





"""---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
""""Model Training"""
# please specify:
# 1. training_steps: 100k
# 2. model_name:
# 3. warmup_steps: 20k

model_name = 'AF2_transformer_l6h512_validation_anti_l7'
if not os.path.exists(model_name):
	os.makedirs(model_name)
with open(model_name + '/print_message.txt', 'w') as f:
	f.write("Number of training data: " + str(X_train.shape[0]) + "\n")

query_size, key_size, value_size, num_hiddens = 512, 512, 512, 512
num_layers, dropout = 6, 0.1
lr, training_steps, batch_size, label_smoothing = 0.0004, 100000, 4096, 0.1
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
warmup = 20000
scheduler = WarmupCosineSchedule(optimizer, warmup, t_total=training_steps)

if torch.cuda.device_count() > 1:
  print("Let's use", torch.cuda.device_count(), "GPUs!")
  model_standard = nn.DataParallel(model_standard)
  with open(model_name + '/print_message.txt', 'a') as f:
  	f.write("Let's use " + str(torch.cuda.device_count()) + " GPUs!" + "\n")

train_seq2seq_training_steps_DP(model_standard, X_train, X_valid_len, Y_train, Y_valid_len, working_score_tensor, lr, training_steps, batch_size, label_smoothing, amino_dict, device, model_name=model_name, warmup=scheduler, optimizer=optimizer)