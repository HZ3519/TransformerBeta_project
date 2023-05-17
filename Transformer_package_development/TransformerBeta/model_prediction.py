import torch
from torch import nn
from d2l import torch as d2l
import numpy as np
from torch.distributions.categorical import Categorical
import os


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



def get_key(val, my_dict):
	"""Return the key of a given value in a dictionary."""

	for key, value in my_dict.items():
			if val == value:
				return key
	return "key doesn't exist"



def unique(x, dim=None):
    """Unique elements of x and indices of those unique elements."""

    unique, inverse = torch.unique(
        x, sorted=True, return_inverse=True, dim=dim)
    perm = torch.arange(inverse.size(0), dtype=inverse.dtype,
                        device=inverse.device)
    inverse, perm = inverse.flip([0]), perm.flip([0])
    return unique, inverse.new_empty(unique.size(0)).scatter_(0, inverse, perm)



def predict_greedy_single(net, target_sequence_raw, amino_dict, num_steps, device, save_attention_weights=False, print_info=False):
	"""Predict for sequence to sequence."""

	softmax_layer = nn.Softmax(dim=2)
	target_sequence = list(target_sequence_raw) 
	
	net.eval()
	target_sequence_unpad = [amino_dict[letter] for letter in target_sequence] + [amino_dict['<eos>']]
	target_valid_len = torch.tensor([len(target_sequence_unpad)], dtype=torch.long, device=device)
	target_sequence = d2l.truncate_pad(target_sequence_unpad, num_steps, amino_dict['<pad>'])

	target_sequence_batch = torch.unsqueeze(torch.tensor(target_sequence, dtype=torch.long, device=device), dim=0)

	enc_outputs = net.encoder(target_sequence_batch, target_valid_len)
	dec_state = net.decoder.init_state(enc_outputs, target_valid_len)

	dec_X = torch.unsqueeze(torch.tensor([amino_dict['<bos>']], dtype=torch.long, device=device), dim=0)
	output_seq, attention_weight_seq = [], []

	prob = 1

	for i in range(num_steps):

		Y_raw, dec_state = net.decoder(dec_X, dec_state)
		Y = softmax_layer(Y_raw)
		dec_X = Y.argmax(dim=2)
		pred = dec_X.squeeze(dim=0).type(torch.int32).item()

		if save_attention_weights:
			attention_weight_seq.append(net.decoder.attention_weights())

		if pred == amino_dict['<eos>']:
			break

		prob_i = torch.max(Y, dim=2).values.squeeze(dim=0).type(torch.float32).item()
		prob *= prob_i

		if print_info == True:
			print("Conditional probability at position {} is {}".format(i+1, prob_i))

		output_seq.append(pred)

	comple_peptide_pred = "".join([get_key(number, amino_dict) for number in output_seq])

	if print_info == True:
		print('Input target sequence is {}, predicted complementary peptide is {}'.format(target_sequence_raw, comple_peptide_pred))
		print("Condition on input, predicted probability is {}".format(prob))

	return comple_peptide_pred, prob, attention_weight_seq



def predict_greedy_batch(net, target_sequence_list, amino_dict, num_steps, device):
	"""Predict for lists of target sequence to lists of complementary peptides sequence."""

	net.eval()
	softmax_layer = nn.Softmax(dim=2)
  
	target_sequence_list_split = [list(sequence) for sequence in target_sequence_list]
	target_sequence_unpad= [[amino_dict[letter] for letter in sequence] + [amino_dict['<eos>']] for sequence in target_sequence_list_split]
	target_valid_len = torch.tensor([len(sequence) for sequence in target_sequence_unpad], dtype=torch.long, device=device)
	target_sequence_batch = torch.tensor([d2l.truncate_pad(sequence, num_steps, amino_dict['<pad>']) for sequence in target_sequence_unpad], dtype=torch.long, device=device)

	enc_outputs = net.encoder(target_sequence_batch, target_valid_len)
	dec_state = net.decoder.init_state(enc_outputs, target_valid_len)

	dec_X = torch.tensor([amino_dict['<bos>']] * target_sequence_batch.shape[0], device=device).reshape(-1, 1)

	prob = torch.tensor([1] * target_sequence_batch.shape[0], device=device).reshape(-1, 1)

	for i in range(num_steps-2):

		Y_raw, dec_state = net.decoder(dec_X, dec_state)
		Y = softmax_layer(Y_raw)
		dec_X = Y.argmax(dim=2)
		Y_pred = dec_X.type(torch.int32) 

		prob_i = torch.max(Y, dim=2).values.squeeze(dim=0)
		prob = torch.mul(prob, prob_i)

	Y_pred = torch.cat((net.decoder.seqX, dec_X), dim=1).type(torch.int32)

	comple_peptide_pred = ["".join([get_key(number, amino_dict) for number in output_seq]) for output_seq in Y_pred[:, 1:]]

	final = zip(target_sequence_list, comple_peptide_pred, prob.reshape(-1).tolist())
	return np.array(list(final))



def predict_greedy_minibatch(net, target_sequence_list, amino_dict, num_steps, device, num_minibatch, print_info=False, output_file=None):

	peptide_pred_list = []

	if len(target_sequence_list)%num_minibatch == 0:
		num_iter = int(len(target_sequence_list)//num_minibatch)
	else:
		num_iter = int(len(target_sequence_list)//num_minibatch + 1)

	for i in range(num_iter):

		if i%10 == 0 and print_info == True:
			print("Iteration: {}".format(i))
			print("Number of peptide predicted: {}".format(num_minibatch*i))
			if output_file != None:
				if os.path.exists(output_file):
					with open(output_file, 'a') as f:
						f.write("Iteration: {}\n".format(i))
						f.write("Number of peptide predicted: {}\n".format(num_minibatch*i))

		if num_minibatch*i+num_minibatch > len(target_sequence_list):
			peptide_pred = predict_greedy_batch(net, target_sequence_list[num_minibatch*i:], amino_dict, num_steps, device)
			peptide_pred_list.append(peptide_pred)
			break

		peptide_pred = predict_greedy_batch(net, target_sequence_list[num_minibatch*i:num_minibatch*i+num_minibatch], amino_dict, num_steps, device)
		peptide_pred_list.append(peptide_pred)

	peptide_pred_array = np.vstack(peptide_pred_list)
	return peptide_pred_array



def evaluate_single(net, target_sequence_raw, peptide_sequence_raw,amino_dict, num_steps, device, save_attention_weights=False, print_info=False):
	"""Predict for sequence to sequence."""

	net.eval()
	softmax_layer = nn.Softmax(dim=2)

	target_sequence = list(target_sequence_raw)
	peptide_sequence = list(peptide_sequence_raw)
	
	target_sequence_unpad = [amino_dict[letter] for letter in target_sequence] + [amino_dict['<eos>']]
	target_valid_len = torch.tensor([len(target_sequence_unpad)], dtype=torch.long, device=device)
	target_sequence = d2l.truncate_pad(target_sequence_unpad, num_steps, amino_dict['<pad>'])
	target_sequence_batch = torch.unsqueeze(torch.tensor(target_sequence, dtype=torch.long, device=device), dim=0)

	peptide_sequence_unpad = [amino_dict[letter] for letter in peptide_sequence] + [amino_dict['<eos>']]
	peptide_valid_len = torch.tensor([len(peptide_sequence_unpad)], dtype=torch.long, device=device)
	peptide_sequence = d2l.truncate_pad(peptide_sequence_unpad, num_steps, amino_dict['<pad>'])
	peptide_sequence_batch = torch.unsqueeze(torch.tensor(peptide_sequence, dtype=torch.long, device=device), dim=0)


	enc_outputs = net.encoder(target_sequence_batch, target_valid_len)
	dec_state = net.decoder.init_state(enc_outputs, target_valid_len)
	dec_X = torch.unsqueeze(torch.tensor([amino_dict['<bos>']], dtype=torch.long, device=device), dim=0)
	output_seq, attention_weight_seq = [], []

	prob = 1

	for i in range(num_steps):
	
		Y_raw, dec_state = net.decoder(dec_X, dec_state)
		Y = softmax_layer(Y_raw)
		dec_X = peptide_sequence_batch[:, i].reshape(1,-1)
		if save_attention_weights:
			attention_weight_seq.append(net.decoder.attention_weights())

		if peptide_sequence_unpad[i] == amino_dict['<eos>']:
			break

		index_i =  peptide_sequence_batch[:, i].squeeze(dim=0).item()
		prob_i = Y[:, :, index_i].squeeze(dim=0).squeeze(dim=0).type(torch.float32).item()
		prob *= prob_i

		if print_info == True:
			print("Conditional probability at position {} is {}".format(i+1, prob_i))

	if print_info == True:
		print('Input target sequence is {}, complementary peptide is {}'.format(target_sequence_raw, peptide_sequence_raw))
		print('Evaluated probability is {}'.format(prob))

	return prob, attention_weight_seq



def evaluate_batch(net, target_sequence_raw, peptide_sequence_raw, amino_dict, num_steps, device):

	net.eval()
	softmax_layer = nn.Softmax(dim=2)

	target_sequence = list(target_sequence_raw)
	peptide_sequence = list(peptide_sequence_raw)
	
	target_sequence_list_split = [list(sequence) for sequence in target_sequence]
	target_sequence_unpad= [[amino_dict[letter] for letter in sequence] + [amino_dict['<eos>']] for sequence in target_sequence_list_split]
	target_valid_len = torch.tensor([len(sequence) for sequence in target_sequence_unpad], dtype=torch.long, device=device)
	target_sequence_batch = torch.tensor([d2l.truncate_pad(sequence, num_steps, amino_dict['<pad>']) for sequence in target_sequence_unpad], dtype=torch.long, device=device)

	peptide_sequence_list_split = [list(sequence) for sequence in peptide_sequence]
	peptide_sequence_unpad= [[amino_dict[letter] for letter in sequence] + [amino_dict['<eos>']] for sequence in peptide_sequence_list_split]
	peptide_valid_len = torch.tensor([len(sequence) for sequence in peptide_sequence_unpad], dtype=torch.long, device=device)
	peptide_sequence_batch = torch.tensor([d2l.truncate_pad(sequence, num_steps, amino_dict['<pad>']) for sequence in peptide_sequence_unpad], dtype=torch.long, device=device)

	enc_outputs = net.encoder(target_sequence_batch, target_valid_len)
	dec_state = net.decoder.init_state(enc_outputs, target_valid_len)

	dec_X = torch.tensor([amino_dict['<bos>']] * target_sequence_batch.shape[0], device=device).reshape(-1, 1) 
	
	prob = torch.tensor([1] * target_sequence_batch.shape[0], device=device).reshape(-1, 1)

	for i in range(num_steps-2): 
		
		Y_raw, dec_state = net.decoder(dec_X, dec_state)
		Y = softmax_layer(Y_raw)
		dec_X = peptide_sequence_batch[:, i].reshape(-1,1)

		index = dec_X.type(torch.int64)
		prob_i = torch.gather(Y, dim=2, index = index.unsqueeze(dim=2))
		prob_i = prob_i.squeeze(dim=2).squeeze(dim=0)
		prob = torch.mul(prob, prob_i)
	
	final = zip(target_sequence_raw, peptide_sequence_raw, prob.reshape(-1).tolist())
	return np.array(list(final)) 



def evaluate_minibatch(net, target_sequence_list, peptide_sequence_list, amino_dict, num_steps, device, num_minibatch, print_info = False, output_file = None):

	target_peptide_list = []

	if len(target_sequence_list)%num_minibatch == 0:
		num_iter = int(len(target_sequence_list)//num_minibatch)
	else:
		num_iter = int(len(target_sequence_list)//num_minibatch + 1)

	for i in range(num_iter):
		if i%10 == 0 and print_info == True:
			print("Iteration: {}".format(i))
			print("Number of sequences evaluated: {}".format(num_minibatch*i))
			if output_file != None:
				if os.path.exists(output_file):
					with open(output_file, 'a') as f:
						f.write("Iteration: {}\n".format(i))
						f.write("Number of sequences evaluated: {}\n".format(num_minibatch*i))

		if num_minibatch*i+num_minibatch > len(target_sequence_list):
			
			target_peptide_eval = evaluate_batch(net, target_sequence_list[num_minibatch*i:], peptide_sequence_list[num_minibatch*i:], amino_dict, num_steps, device)
			target_peptide_list.append(target_peptide_eval)
			break

		target_peptide_eval = evaluate_batch(net, target_sequence_list[num_minibatch*i:num_minibatch*i+num_minibatch], peptide_sequence_list[num_minibatch*i:num_minibatch*i+num_minibatch], amino_dict, num_steps, device)
		target_peptide_list.append(target_peptide_eval)

	target_peptide_array = np.vstack(target_peptide_list)
	return target_peptide_array



def sample_candidates(net, target_sequence_raw, num_candidates, amino_dict, num_steps, device, max_iter=100):
	"""Given one target sequence, sample a list of unique probable candidates."""

	bos = amino_dict['<bos>']
	eos = amino_dict['<eos>']
	pad = amino_dict['<pad>']
	unk = amino_dict['<unk>']
	
	net.eval()
	softmax_layer = nn.Softmax(dim=2)
	target_sequence = list(target_sequence_raw)

	num_remaining = num_candidates

	Y_pred_track = torch.tensor([1] *(num_steps-2), device=device, dtype=torch.int32).reshape(1, -1) 
	prob_track = torch.tensor([1], device=device, dtype=torch.long).reshape(1, 1)

	target_sequence_unpad = [amino_dict[letter] for letter in target_sequence] + [amino_dict['<eos>']]

	samples_total = 0

	for j in range(max_iter):

		target_valid_len = torch.tensor([len(target_sequence_unpad)], dtype=torch.long, device=device).repeat_interleave(num_remaining*2)
		target_sequence = torch.tensor(d2l.truncate_pad(target_sequence_unpad, num_steps, amino_dict['<pad>']), dtype=torch.long, device=device)
		target_sequence_batch = target_sequence.repeat(num_remaining*2, 1)

		enc_outputs = net.encoder(target_sequence_batch, target_valid_len)
		dec_state = net.decoder.init_state(enc_outputs, target_valid_len)

		dec_X = torch.tensor([amino_dict['<bos>']] * target_sequence_batch.shape[0], device=device).reshape(-1, 1) 
		
		prob = torch.tensor([1] * target_sequence_batch.shape[0], dtype = torch.long, device=device).reshape(-1, 1)


		for i in range(num_steps-2): 

			Y_raw, dec_state = net.decoder(dec_X, dec_state)
			Y = softmax_layer(Y_raw)
			
			m = Categorical(probs=Y)
			dec_X = m.sample()
			Y_pred = dec_X.type(torch.int32) 
			
			index = Y_pred.type(torch.int64)
	
			prob_i = torch.gather(Y, dim=2, index = index.unsqueeze(dim=2))
			prob_i = prob_i.squeeze(dim=2).squeeze(dim=0)
			prob = torch.mul(prob, prob_i)
		
		Y_pred = torch.cat((net.decoder.seqX, dec_X), dim=1).type(torch.int32)
		samples_total += Y_pred.shape[0]

		#Y_raw, dec_state = net.decoder(dec_X, dec_state)
		#Y = softmax_layer(Y_raw)
		#prob_i = torch.max(Y, dim=2).values.squeeze(dim=0)
		#prob = torch.mul(prob, prob_i)

		Y_pred = Y_pred[:, 1:] 
		Y_pred_track = torch.cat((Y_pred_track, Y_pred), dim=0)
		prob_track = torch.cat((prob_track, prob), dim=0)

		if j==0:
			Y_pred_track = Y_pred_track[1:, :]
			prob_track = prob_track[1:, :]

		Y_pred_track, unique_indices = unique(Y_pred_track, dim=0)
		prob_track = prob_track[unique_indices] 
		
		prob_track = prob_track[~(Y_pred_track == bos).any(1), :]
		Y_pred_track = Y_pred_track[~(Y_pred_track == bos).any(1), :]

		prob_track = prob_track[~(Y_pred_track == eos).any(1), :]
		Y_pred_track = Y_pred_track[~(Y_pred_track == eos).any(1), :]

		prob_track = prob_track[~(Y_pred_track == pad).any(1), :]
		Y_pred_track = Y_pred_track[~(Y_pred_track == pad).any(1), :]

		prob_track = prob_track[~(Y_pred_track == unk).any(1), :]
		Y_pred_track = Y_pred_track[~(Y_pred_track == unk).any(1), :]		

		if Y_pred_track.shape[0] >= num_candidates:

			comple_peptide_pred = ["".join([get_key(number, amino_dict) for number in output_seq]) for output_seq in Y_pred_track] 

			prob_track = prob_track.reshape(-1).tolist()

			dtype = [('peptide', '<U21'), ('prob', float)]
			values = list(zip(comple_peptide_pred, prob_track))

			a = np.array(values, dtype=dtype)

			final = list(np.sort(a, order='prob')[-1::-1])
			final = final[0:num_candidates]
			final_array = np.array([[str, prob] for str, prob in final])

			print("number of total candidates sampled: {}".format(samples_total))
			print("number of unique top candidates successfully sampled: {}".format(num_candidates))
			return final_array
			
		num_remaining = num_candidates - Y_pred_track.shape[0]

	comple_peptide_pred = ["".join([get_key(number, amino_dict) for number in output_seq]) for output_seq in Y_pred_track[:, 1:-1]]

	prob_track = prob_track.reshape(-1).tolist()

	dtype = [('peptide', '<U21'), ('prob', float)]
	values = list(zip(comple_peptide_pred, prob_track))

	a = np.array(values, dtype=dtype)
	final = list(np.sort(a, order='prob')[-1::-1])
	final_array = np.array([[str, prob] for str, prob in final])

	print("number of total candidates sampled: {}".format(samples_total))
	print("number of unique candidates successfully sampled: {}".format(len(final)))
	return final_array



def sample_single_candidate(net, target_sequence_raw, amino_dict, num_steps, device, max_iter=100):
	"""Given one target sequence, sample a list of unique probable candidates."""

	bos = amino_dict['<bos>']
	eos = amino_dict['<eos>']
	pad = amino_dict['<pad>']
	unk = amino_dict['<unk>']
	
	net.eval()
	softmax_layer = nn.Softmax(dim=2)
	target_sequence = list(target_sequence_raw)

	num_candidates = 1
	num_remaining = num_candidates

	Y_pred_track = torch.tensor([1] *(num_steps-2), device=device, dtype=torch.int32).reshape(1, -1) 
	prob_track = torch.tensor([1], device=device, dtype=torch.long).reshape(1, 1)

	target_sequence_unpad = [amino_dict[letter] for letter in target_sequence] + [amino_dict['<eos>']]

	samples_total = 0

	for j in range(max_iter):

		target_valid_len = torch.tensor([len(target_sequence_unpad)], dtype=torch.long, device=device).repeat_interleave(num_remaining)
		target_sequence = torch.tensor(d2l.truncate_pad(target_sequence_unpad, num_steps, amino_dict['<pad>']), dtype=torch.long, device=device)
		target_sequence_batch = target_sequence.repeat(num_remaining, 1)

		enc_outputs = net.encoder(target_sequence_batch, target_valid_len)
		dec_state = net.decoder.init_state(enc_outputs, target_valid_len)

		dec_X = torch.tensor([amino_dict['<bos>']] * target_sequence_batch.shape[0], device=device).reshape(-1, 1) 
		
		prob = torch.tensor([1] * target_sequence_batch.shape[0], dtype = torch.long, device=device).reshape(-1, 1)


		for i in range(num_steps-2): 

			Y_raw, dec_state = net.decoder(dec_X, dec_state)
			Y = softmax_layer(Y_raw)
			
			m = Categorical(probs=Y)
			dec_X = m.sample()
			Y_pred = dec_X.type(torch.int32) 
			
			index = Y_pred.type(torch.int64)
	
			prob_i = torch.gather(Y, dim=2, index = index.unsqueeze(dim=2))
			prob_i = prob_i.squeeze(dim=2).squeeze(dim=0)
			prob = torch.mul(prob, prob_i)
		
		Y_pred = torch.cat((net.decoder.seqX, dec_X), dim=1).type(torch.int32)
		samples_total += Y_pred.shape[0]

		#Y_raw, dec_state = net.decoder(dec_X, dec_state)
		#Y = softmax_layer(Y_raw)
		#prob_i = torch.max(Y, dim=2).values.squeeze(dim=0)
		#prob = torch.mul(prob, prob_i)

		Y_pred = Y_pred[:, 1:] 
		Y_pred_track = torch.cat((Y_pred_track, Y_pred), dim=0)
		prob_track = torch.cat((prob_track, prob), dim=0)

		if j==0:
			Y_pred_track = Y_pred_track[1:, :]
			prob_track = prob_track[1:, :]

		Y_pred_track, unique_indices = unique(Y_pred_track, dim=0)
		prob_track = prob_track[unique_indices] 
		
		prob_track = prob_track[~(Y_pred_track == bos).any(1), :]
		Y_pred_track = Y_pred_track[~(Y_pred_track == bos).any(1), :]

		prob_track = prob_track[~(Y_pred_track == eos).any(1), :]
		Y_pred_track = Y_pred_track[~(Y_pred_track == eos).any(1), :]

		prob_track = prob_track[~(Y_pred_track == pad).any(1), :]
		Y_pred_track = Y_pred_track[~(Y_pred_track == pad).any(1), :]

		prob_track = prob_track[~(Y_pred_track == unk).any(1), :]
		Y_pred_track = Y_pred_track[~(Y_pred_track == unk).any(1), :]		

		if Y_pred_track.shape[0] == num_candidates:

			comple_peptide_pred = ["".join([get_key(number, amino_dict) for number in output_seq]) for output_seq in Y_pred_track] 

			prob_track = prob_track.reshape(-1).tolist()

			dtype = [('peptide', '<U21'), ('prob', float)]
			values = list(zip(comple_peptide_pred, prob_track))

			a = np.array(values, dtype=dtype)

			final = list(np.sort(a, order='prob')[-1::-1])
			final = final[0:num_candidates]
			final_array = np.array([[str, prob] for str, prob in final])

			print("number of total candidates sampled: {}".format(samples_total))
			print("number of unique top candidates successfully sampled: {}".format(num_candidates))
			return final_array
			
		num_remaining = num_candidates - Y_pred_track.shape[0]

	comple_peptide_pred = ["".join([get_key(number, amino_dict) for number in output_seq]) for output_seq in Y_pred_track[:, 1:-1]]

	prob_track = prob_track.reshape(-1).tolist()

	dtype = [('peptide', '<U21'), ('prob', float)]
	values = list(zip(comple_peptide_pred, prob_track))

	a = np.array(values, dtype=dtype)
	final = list(np.sort(a, order='prob')[-1::-1])
	final_array = np.array([[str, prob] for str, prob in final])

	print("number of total candidates sampled: {}".format(samples_total))
	print("number of unique candidates successfully sampled: {}".format(len(final)))
	return final_array