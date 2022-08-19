import numpy as np
import torch
from d2l import torch as d2l
import time 

from TransformerBeta import *

"""---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
""""Build the model -- use corresponding settings"""
query_size, key_size, value_size, num_hiddens = 512, 512, 512, 512
num_layers, dropout = 4, 0.1
lr, num_epochs, batch_size, label_smoothing = 0.0004, 500, 6000, 0.1
ffn_num_input, ffn_num_hiddens, num_heads = 512, 2048, 8

norm_shape = [512] # 32 corresponds to the dim of such number to normalize
device = d2l.try_gpu()


encoder_wide = TransformerEncoder(
	len(amino_dict), key_size, query_size, value_size, num_hiddens, 
	norm_shape, ffn_num_input, ffn_num_hiddens, num_heads,
	num_layers, dropout)
decoder_wide = TransformerDecoder(
	len(amino_dict), key_size, query_size, value_size, num_hiddens, 
	norm_shape, ffn_num_input, ffn_num_hiddens, num_heads,
	num_layers, dropout)
model_wide = EncoderDecoder(encoder_wide, decoder_wide)


model_wide_total_params = sum(p.numel() for p in model_wide.parameters())
model_wide_total_trainable_params = sum(p.numel() for p in model_wide.parameters() if p.requires_grad)

print('Wide model: total number of parameters: {}'.format(model_wide_total_params))
print('Wide model: total number of trainable parameters: {}'.format(model_wide_total_trainable_params))


model_wide.load_state_dict(torch.load("model_wide_22Jul16_1011AM", map_location = ('cpu')))


"""---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
""""sampling"""
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

task_target = 'EQVTNVGG' 
model_use = model_wide
model_name = 'model_wide'
file_time = time.strftime('%y%b%d_%I%M%p', time.gmtime())

prediction_length = len(task_target)
num_candidates = 2000
max_iter = 20

peptide_candidates = sample_candidates(model_use, task_target, num_candidates, amino_dict, prediction_length + 2, device, max_iter=max_iter)
np.savetxt('{}_{}_{}.csv'.format(task_target, model_name, file_time), peptide_candidates, delimiter=',', fmt='%s')