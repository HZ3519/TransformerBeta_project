import pickle
import numpy as np
import torch
from d2l import torch as d2l

from TransformerBeta import *

"""---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
""""Data preprocessing"""

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


# target, complementary_seq, counts, promiscuity, length, working_score, hb_pattern, para/anti(0/1)
BSn_data = []
least_length = 3

for frag_i_data in data_list_parallel[least_length-1:]:
	for keys in frag_i_data.keys():

		length = len(keys)
		for element in frag_i_data[keys]:

			working_score = length**2 * element.count_score - 0.01 * length * element.promiscuity_score
			list_i = [keys, element.complementary_sequence, element.count_score, element.promiscuity_score, length, working_score, element.hb_pattern, 0]
			BSn_data.append(list_i)

for frag_i_data in data_list_antiparallel[least_length-1:]:
	for keys in frag_i_data.keys():
		length = len(keys)
		for element in frag_i_data[keys]:

			working_score = length**2 * element.count_score - 0.01 * length * element.promiscuity_score
			list_i = [keys, element.complementary_sequence, element.count_score, element.promiscuity_score, length, working_score, element.hb_pattern, 1]
			BSn_data.append(list_i)
# target, complementary_seq, counts, promiscuity, length, working_score, hb_pattern, para/anti

BSn_data_dataset_sequence = np.array(BSn_data, dtype=object)[:, 0:2]
scores_array = np.array(BSn_data, dtype=object)[:, 5].reshape(-1, 1)
BSn_data_dataset_scores = np.hstack([BSn_data_dataset_sequence, scores_array])
BSn_data_dataset1 = BSn_data_dataset_scores[BSn_data_dataset_scores[:, 2] >= 0]
working_score_tensor = torch.tensor(list(BSn_data_dataset1[:, 2]))


num_steps_training = 10
X_train, X_valid_len, Y_train, Y_valid_len = preprocess_train(BSn_data_dataset1[:, 0], BSn_data_dataset1[:, 1], amino_dict, num_steps_training)






"""---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
""""Model Training"""

query_size, key_size, value_size, num_hiddens = 512, 512, 512, 512
num_layers, dropout = 4, 0.1
lr, num_epochs, batch_size, label_smoothing = 0.0004, 500, 6000, 0.1
ffn_num_input, ffn_num_hiddens, num_heads = 512, 2048, 8

norm_shape = [512] 
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

print('Base model: total number of parameters: {}'.format(model_wide_total_params))
print('Base model: total number of trainable parameters: {}'.format(model_wide_total_trainable_params))


train_seq2seq(model_wide, X_train, X_valid_len, Y_train, Y_valid_len, working_score_tensor, lr, num_epochs, batch_size, label_smoothing, amino_dict, device, model_name='model_wide', warmup=35000)