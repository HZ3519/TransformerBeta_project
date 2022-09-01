import pickle
import numpy as np
import torch
from d2l import torch as d2l

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


# split the data in training and validation
num_steps_training = 10 

X_train, X_valid_len, Y_train, Y_valid_len, X_validation, X_validation_valid_len, Y_validation, Y_validation_valid_len = preprocess_train(X_train_letter, Y_train_letter, amino_dict, num_steps_training, X_validation_letter=X_validation_letter, Y_validation_letter=Y_validation_letter)

working_score_tensor = torch.tensor(list(np.delete(BSn_data_dataset2[:, 5], validation_indices, axis=0)))






"""---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
""""Model Training"""

query_size, key_size, value_size, num_hiddens = 512, 512, 512, 512
num_layers, dropout = 4, 0.1
lr, num_epochs, batch_size, label_smoothing = 0.0004, 500, 6000, 0.1
ffn_num_input, ffn_num_hiddens, num_heads = 512, 2048, 8

norm_shape = [512] # 32 corresponds to the dim of such number to normalize
device = d2l.try_gpu()

encoder_base = TransformerEncoder(
	len(amino_dict), key_size, query_size, value_size, num_hiddens, 
	norm_shape, ffn_num_input, ffn_num_hiddens, num_heads,
	num_layers, dropout)
decoder_base = TransformerDecoder(
	len(amino_dict), key_size, query_size, value_size, num_hiddens, 
	norm_shape, ffn_num_input, ffn_num_hiddens, num_heads,
	num_layers, dropout)
model_wide = EncoderDecoder(encoder_base, decoder_base)


model_wide_total_params = sum(p.numel() for p in model_wide.parameters())
model_wide_total_trainable_params = sum(p.numel() for p in model_wide.parameters() if p.requires_grad)

print('Wide model: total number of parameters: {}'.format(model_wide_total_params))
print('Wide model: total number of trainable parameters: {}'.format(model_wide_total_trainable_params))


train_seq2seq(model_wide, X_train, X_valid_len, Y_train, Y_valid_len, working_score_tensor, lr, num_epochs, batch_size, label_smoothing, amino_dict, device, model_name='model_wide_validation', warmup=35000, X_validation=X_validation, Y_validation=Y_validation, X_validation_valid_len=X_validation_valid_len, Y_validation_valid_len=Y_validation_valid_len)