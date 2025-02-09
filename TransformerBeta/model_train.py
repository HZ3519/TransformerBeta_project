import torch
from torch import nn
from d2l import torch as d2l
from TransformerBeta.model_architecture import sequence_mask
from TransformerBeta.model_metric import hamming_distance_training
import time 
import matplotlib.pyplot as plt 
import numpy as np
from torch.optim.lr_scheduler import LambdaLR
import math



class WarmupCosineSchedule(LambdaLR):
    """ Linear warmup and then cosine decay.
        Linearly increases learning rate from 0 to 1 over `warmup_steps` training steps.
        Decreases learning rate from 1. to 0. over remaining `t_total - warmup_steps` steps following a cosine curve.
        If `cycles` (default=0.5) is different from default, learning rate follows cosine function after warmup.
    """
    def __init__(self, optimizer, warmup_steps, t_total, cycles=.5, last_epoch=-1):
        self.warmup_steps = warmup_steps
        self.t_total = t_total
        self.cycles = cycles
        super(WarmupCosineSchedule, self).__init__(optimizer, self.lr_lambda, last_epoch=last_epoch)

    def lr_lambda(self, step):
        if step < self.warmup_steps:
            return float(step) / float(max(1.0, self.warmup_steps))
        # progress after warmup
        progress = float(step - self.warmup_steps) / float(max(1, self.t_total - self.warmup_steps))
        return max(0.0, 0.5 * (1. + math.cos(math.pi * float(self.cycles) * 2.0 * progress)))



class MaskedSoftmaxCELoss(nn.CrossEntropyLoss):
	"""The softmax cross-entropy loss with masks."""

	def forward(self, pred, label, valid_len, sample_weights, label_smoothing):
		masks = torch.ones_like(label)
		masks = sequence_mask(masks, valid_len)
		self.reduction = 'none'
		self.label_smoothing = label_smoothing
		unmaskeded_loss = super(MaskedSoftmaxCELoss, self).forward(pred.permute(0, 2, 1), label)
		masked_loss = (unmaskeded_loss * masks).mean(dim=1) # mean of each batch along num_steps
		
		if sample_weights is None:
			return masked_loss.mean()
		weighted_masked_loss = masked_loss * sample_weights
		weighted_average_masked_loss = weighted_masked_loss.sum() / sample_weights.sum()  # weighted average across batch
		return weighted_average_masked_loss



def preprocess_train(X_train_letter, Y_train_letter, amino_dict, num_steps, X_validation_letter=None, Y_validation_letter=None):

	X_train_letter_split = [list(sequence) for sequence in X_train_letter]
	Y_train_letter_split = [list(sequence) for sequence in Y_train_letter]

	X_train_unpad= [[amino_dict[letter] for letter in sequence] + [amino_dict['<eos>']] for sequence in X_train_letter_split]
	Y_train_unpad = [[amino_dict[letter] for letter in sequence] + [amino_dict['<eos>']] for sequence in Y_train_letter_split]

	X_valid_len = torch.tensor([len(sequence) for sequence in X_train_unpad])
	Y_valid_len = torch.tensor([len(sequence) for sequence in Y_train_unpad])

	X_train = torch.tensor([d2l.truncate_pad(sequence, num_steps, amino_dict['<pad>']) for sequence in X_train_unpad])
	Y_train = torch.tensor([d2l.truncate_pad(sequence, num_steps, amino_dict['<pad>']) for sequence in Y_train_unpad])

	if X_validation_letter is not None:

		num_steps = len(X_validation_letter[0]) + 2

		X_validation_letter_split = [list(sequence) for sequence in X_validation_letter]
		Y_validation_letter_split = [list(sequence) for sequence in Y_validation_letter]

		X_validation_unpad= [[amino_dict[letter] for letter in sequence] + [amino_dict['<eos>']] for sequence in X_validation_letter_split]
		Y_validation_unpad = [[amino_dict[letter] for letter in sequence] + [amino_dict['<eos>']] for sequence in Y_validation_letter_split]

		X_validation_valid_len = torch.tensor([len(sequence) for sequence in X_validation_unpad])
		Y_validation_valid_len = torch.tensor([len(sequence) for sequence in Y_validation_unpad])

		X_validation = torch.tensor([d2l.truncate_pad(sequence, num_steps, amino_dict['<pad>']) for sequence in X_validation_unpad])
		Y_validation = torch.tensor([d2l.truncate_pad(sequence, num_steps, amino_dict['<pad>']) for sequence in Y_validation_unpad])
		
		return X_train, X_valid_len, Y_train, Y_valid_len, X_validation, X_validation_valid_len, Y_validation, Y_validation_valid_len

	return X_train, X_valid_len, Y_train, Y_valid_len



def seq2num(X_train_letter, Y_train_letter, amino_dict, X_validation_letter=None, Y_validation_letter=None):

	X_train_letter_split = [list(sequence) for sequence in X_train_letter]
	Y_train_letter_split = [list(sequence) for sequence in Y_train_letter]

	X_train= np.array([[amino_dict[letter] for letter in sequence] for sequence in X_train_letter_split])
	Y_train = np.array([[amino_dict[letter] for letter in sequence] for sequence in Y_train_letter_split])

	if X_validation_letter is not None:

		X_validation_letter_split = [list(sequence) for sequence in X_validation_letter]
		Y_validation_letter_split = [list(sequence) for sequence in Y_validation_letter]

		X_validation = np.array([[amino_dict[letter] for letter in sequence] for sequence in X_validation_letter_split])
		Y_validation = np.array([[amino_dict[letter] for letter in sequence] for sequence in Y_validation_letter_split])

		return X_train, Y_train, X_validation, Y_validation

	return X_train, Y_train



def preprocess_train_T5(X_train_letter, Y_train_letter, amino_dict, num_steps, X_validation_letter=None, Y_validation_letter=None):

	X_train_letter_split = [list(sequence) for sequence in X_train_letter]
	Y_train_letter_split = [list(sequence) for sequence in Y_train_letter]

	X_train_unpad= [[amino_dict[letter] for letter in sequence] + [amino_dict['<eos>']] for sequence in X_train_letter_split]
	Y_train_unpad = [[amino_dict[letter] for letter in sequence] + [amino_dict['<eos>']] for sequence in Y_train_letter_split]

	X_valid_len = torch.tensor([[1 for i in range(len(sequence))] + [0 for i in range(num_steps-len(sequence))] for sequence in X_train_unpad])
	Y_valid_len = torch.tensor([len(sequence) for sequence in Y_train_unpad])

	X_train = torch.tensor([d2l.truncate_pad(sequence, num_steps, amino_dict['<pad>']) for sequence in X_train_unpad])
	Y_train = torch.tensor([d2l.truncate_pad(sequence, num_steps, amino_dict['<pad>']) for sequence in Y_train_unpad])

	if X_validation_letter is not None:

		num_steps = len(X_validation_letter[0]) + 2

		X_validation_letter_split = [list(sequence) for sequence in X_validation_letter]
		Y_validation_letter_split = [list(sequence) for sequence in Y_validation_letter]

		X_validation_unpad = [[amino_dict[letter] for letter in sequence] + [amino_dict['<eos>']] for sequence in X_validation_letter_split]
		Y_validation_unpad = [[amino_dict[letter] for letter in sequence] + [amino_dict['<eos>']] for sequence in Y_validation_letter_split]

		X_validation_valid_len = torch.tensor([[1 for i in range(len(sequence))] + [0 for i in range(num_steps-len(sequence))] for sequence in X_validation_unpad])
		Y_validation_valid_len = torch.tensor([len(sequence) for sequence in Y_validation_unpad])
		
		X_validation = torch.tensor([d2l.truncate_pad(sequence, num_steps, amino_dict['<pad>']) for sequence in X_validation_unpad])
		Y_validation = torch.tensor([d2l.truncate_pad(sequence, num_steps, amino_dict['<pad>']) for sequence in Y_validation_unpad])
		
		return X_train, X_valid_len, Y_train, Y_valid_len, X_validation, X_validation_valid_len, Y_validation, Y_validation_valid_len

	return X_train, X_valid_len, Y_train, Y_valid_len

	

def train_seq2seq(net, X_train, X_valid_len, Y_train, Y_valid_len, sample_weights, lr, num_epochs, batch_size, label_smoothing, amino_dict, device, optimizer=None, warmup=1, model_name = 'model_demo', X_validation=None, X_validation_valid_len=None, Y_validation=None, Y_validation_valid_len=None):
	"""Train a model for sequence to sequence."""

	def init_weights(module):
		if isinstance(module, (nn.Linear)):
			nn.init.normal_(module.weight, mean=0.0, std=0.01) 
		if isinstance(module, (nn.Embedding)):
			nn.init.normal_(module.weight, mean=0.0, std=1e-5) 

	net.apply(init_weights)
	net.to(device)
	if optimizer is None:
		optimizer = torch.optim.Adam(net.parameters(), lr=lr, betas=(0.9, 0.98), eps = 1.0e-9)
	else:
		optimizer = optimizer
	if isinstance(warmup, int):
		scheduler = torch.optim.lr_scheduler.LinearLR(optimizer, total_iters=warmup, start_factor=1/warmup)
	else:
		scheduler = warmup
	loss = MaskedSoftmaxCELoss()
	net.train()

	animator = d2l.Animator(xlabel='epoch', ylabel='loss')
	training_loss = []
	validation_score_hamming = []
	validation_score_sequence_accuracy = []
	validation_score_average_log_prob = []


	for epoch in range(num_epochs):
		
		net.train()

		timer = d2l.Timer()
		metric = d2l.Accumulator(2)  

		index = torch.randperm(X_train.shape[0])
		X_train_shuffled = X_train[index]
		X_valid_len_shuffled = X_valid_len[index]
		Y_train_shuffled = Y_train[index]
		Y_valid_len_shuffled = Y_valid_len[index]
		sample_weights_shuffled = sample_weights[index]

		X_train_batch = torch.split(X_train_shuffled, batch_size)
		X_valid_len_batch = torch.split(X_valid_len_shuffled, batch_size)
		Y_train_batch = torch.split(Y_train_shuffled, batch_size)
		Y_valid_len_batch = torch.split(Y_valid_len_shuffled, batch_size)
		sample_weights_batch = torch.split(sample_weights_shuffled, batch_size)

		for batch in zip(X_train_batch, X_valid_len_batch, Y_train_batch, Y_valid_len_batch, sample_weights_batch):
			
			optimizer.zero_grad()
			X_train_minibatch, X_valid_len_minibatch, Y_train_minibatch, Y_valid_len_minibatch, sample_weights_minibatch = [x.to(device) for x in batch]
			bos = torch.tensor([amino_dict['<bos>']] * Y_train_minibatch.shape[0],
								device=device).reshape(-1, 1)
			dec_input = torch.cat([bos, Y_train_minibatch[:, :-1]], 1)  
			
			Y_hat, _ = net(X_train_minibatch, dec_input, X_valid_len_minibatch)
			l = loss(Y_hat, Y_train_minibatch, Y_valid_len_minibatch, sample_weights_minibatch, label_smoothing)
			l.backward()  
			d2l.grad_clipping(net, 1)
			num_tokens = Y_valid_len_minibatch.sum()
			optimizer.step()
			scheduler.step()
			with torch.no_grad():
				metric.add(l*num_tokens, num_tokens) 

		animator.add(epoch + 1, (metric[0] / metric[1],))
		training_loss.append(metric[0] / metric[1])
		print("epoch {}, loss: {}".format(epoch+1, metric[0] / metric[1]))
		# save the print message to a a txt file in model name folder
		with open(model_name + '/print_message.txt', 'a') as f:
			f.write("epoch {}, loss: {}".format(epoch+1, metric[0] / metric[1]) + '\n')

		if X_validation is not None:
			net.eval()

			with torch.no_grad():

				softmax_layer = nn.Softmax(dim=2)
				X_validation = X_validation.to(device)
				X_validation_valid_len = X_validation_valid_len.to(device)
				Y_validation = Y_validation.to(device)
				Y_validation_valid_len = Y_validation_valid_len.to(device)

				bos = torch.tensor([amino_dict['<bos>']] * Y_validation.shape[0],
									device=device).reshape(-1, 1)
				Y_true = torch.cat([bos, Y_validation[:, :-1]], 1)  
				Y_true = Y_true.type(torch.int32)

				enc_outputs = net.encoder(X_validation, X_validation_valid_len)
				dec_state = net.decoder.init_state(enc_outputs, X_validation_valid_len)

				dec_X = torch.tensor([amino_dict['<bos>']] * X_validation.shape[0],
									device=device).reshape(-1, 1)

				for i in range(Y_validation.shape[1]-1): 

						Y_raw, dec_state = net.decoder(dec_X, dec_state)
						Y = softmax_layer(Y_raw)
						dec_X = Y.argmax(dim=2)
						Y_pred = dec_X.type(torch.int32)

				Y_pred = torch.cat((net.decoder.seqX, dec_X), dim=1).type(torch.int32)
				
				hamming_scores = hamming_distance_training(Y_pred, Y_true)
				validation_score_hamming.append(hamming_scores.cpu().numpy().item())
				print("epoch {}, hamming distance: {}".format(epoch+1,hamming_scores))
				# save the print message to a a txt file in model name folder
				with open(model_name + '/print_message.txt', 'a') as f:
					f.write("epoch {}, hamming distance: {}".format(epoch+1,hamming_scores) + '\n')

				sequence_accuracy_list = [torch.equal(a, b) for a, b in zip(Y_pred, Y_true)]
				sequence_accuracy = sum(sequence_accuracy_list) / Y_true.shape[0]
				validation_score_sequence_accuracy.append(sequence_accuracy)
				print("epoch {}, sequence accuracy: {}".format(epoch+1,sequence_accuracy))
				# save the print message to a a txt file in model name folder
				with open(model_name + '/print_message.txt', 'a') as f:
					f.write("epoch {}, sequence accuracy: {}".format(epoch+1,sequence_accuracy) + '\n')

				enc_outputs = net.encoder(X_validation, X_validation_valid_len)
				dec_state = net.decoder.init_state(enc_outputs, X_validation_valid_len)
				bos = torch.tensor([amino_dict['<bos>']] * Y_validation.shape[0], device=device).reshape(-1, 1)
				dec_X = torch.tensor([amino_dict['<bos>']] * X_validation.shape[0],
									device=device).reshape(-1, 1)
				prob = torch.tensor([1] * Y_validation.shape[0], device=device).reshape(-1, 1)
				
				for i in range(Y_validation.shape[1]-1):
		
					Y_raw, dec_state = net.decoder(dec_X, dec_state)
					Y = softmax_layer(Y_raw)
					dec_X = Y_validation[:, i].reshape(-1,1)
					index = dec_X.type(torch.int64)
					prob_i = torch.gather(Y, dim=2, index = index.unsqueeze(dim=2))
					prob_i = prob_i.squeeze(dim=2).squeeze(dim=0)
					prob = torch.mul(prob, prob_i)
				
				average_log_prob = torch.log(prob).mean()
				validation_score_average_log_prob.append(average_log_prob.cpu().numpy().item())
				print("epoch {}, average log prob: {}".format(epoch+1,average_log_prob))
				# save the print message to a a txt file in model name folder
				with open(model_name + '/print_message.txt', 'a') as f:
					f.write("epoch {}, average log prob: {}".format(epoch+1,average_log_prob) + '\n')


	print(f'loss {metric[0] / metric[1]:.3f}, {metric[1] / timer.stop():.1f} tokens/sec on {str(device)}')

	# save training loss plot
	plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
	x = list(range(1, num_epochs+1))
	file_time = time.strftime('%y%b%d_%I%M%p', time.gmtime())

	plt.plot(x,training_loss,label="weighted_average_loss")
	plt.xlabel('Epochs')
	plt.ylabel('Weighted average loss')
	plt.title("transformer <{}> training loss".format(model_name))
	plt.grid()
	plt.legend()
	plt.savefig(model_name + '/{}_lossplot_{}.png'.format(model_name, file_time))

	if X_validation is not None:
		# save validation hamming metic plot
		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = list(range(1, num_epochs+1))
		plt.plot(x,validation_score_hamming,label="hamming_distance")
		plt.xlabel('Epochs')
		plt.ylabel('Hamming distance sum')
		plt.title("transformer <{}> validation hamming scores".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_hamming_{}.png'.format(model_name, file_time))
		
		# save validation sequence accuracy plot
		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = list(range(1, num_epochs+1))
		plt.plot(x,validation_score_sequence_accuracy,label="sequence_accuracy")
		plt.xlabel('Epochs')
		plt.ylabel('Sequence accuracy')
		plt.title("transformer <{}> validation sequence accuracy".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_sequence_accuracy_{}.png'.format(model_name, file_time))

		# save validation average log prob plot
		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = list(range(1, num_epochs+1))
		plt.plot(x,validation_score_average_log_prob,label="average_log_prob")
		plt.xlabel('Epochs')
		plt.ylabel('Average log prob')
		plt.title("transformer <{}> validation average log prob".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_average_log_prob_{}.png'.format(model_name, file_time))

	# save model weights
	file_name = model_name + '_' + file_time
	torch.save(net.state_dict(), model_name + '/' + file_name)

	# save training_loss, validation_score_hamming, validation_score_sequence_accuracy
	np.savetxt(model_name + '/{}_training_loss_{}.csv'.format(model_name, file_time), training_loss, delimiter=',', fmt='%s')
	if X_validation is not None:
		np.savetxt(model_name + '/{}_validation_score_hamming_{}.csv'.format(model_name, file_time), validation_score_hamming, delimiter=',', fmt='%s')
		np.savetxt(model_name + '/{}_validation_score_sequence_accuracy_{}.csv'.format(model_name, file_time), validation_score_sequence_accuracy, delimiter=',', fmt='%s')
		np.savetxt(model_name + '/{}_validation_score_average_log_prob_{}.csv'.format(model_name, file_time), validation_score_average_log_prob, delimiter=',', fmt='%s')



def train_seq2seq_training_steps(net, X_train, X_valid_len, Y_train, Y_valid_len, sample_weights, lr, training_steps, batch_size, label_smoothing, amino_dict, device, optimizer=None, warmup=1, model_name = 'model_demo', X_validation=None, X_validation_valid_len=None, Y_validation=None, Y_validation_valid_len=None):
	"""Train a model for sequence to sequence."""

	def init_weights(module):
		if isinstance(module, (nn.Linear)):
			nn.init.normal_(module.weight, mean=0.0, std=0.01)
		if isinstance(module, (nn.Embedding)):
			nn.init.normal_(module.weight, mean=0.0, std=1e-5)

	net.apply(init_weights)
	net.to(device)
	if optimizer is None:
		optimizer = torch.optim.Adam(net.parameters(), lr=lr, betas=(0.9, 0.98), eps = 1.0e-9)
	else:
		optimizer = optimizer
	if isinstance(warmup, int):
		scheduler = torch.optim.lr_scheduler.LinearLR(optimizer, total_iters=warmup, start_factor=1/warmup)
	else:
		scheduler = warmup
	loss = MaskedSoftmaxCELoss()
	net.train()

	animator = d2l.Animator(xlabel='training steps', ylabel='loss')
	training_loss = []
	validation_score_hamming = []
	validation_score_sequence_accuracy = []
	validation_score_average_log_prob = []

	# here we modify the training steps to be the number of epochs
	current_step = 0
	training_steps_count = []
	
	while current_step < training_steps:

		net.train()

		timer = d2l.Timer()
		metric = d2l.Accumulator(2)  

		index = torch.randperm(X_train.shape[0])
		X_train_shuffled = X_train[index]
		X_valid_len_shuffled = X_valid_len[index]
		Y_train_shuffled = Y_train[index]
		Y_valid_len_shuffled = Y_valid_len[index]
		sample_weights_shuffled = sample_weights[index]

		X_train_batch = torch.split(X_train_shuffled, batch_size)
		X_valid_len_batch = torch.split(X_valid_len_shuffled, batch_size)
		Y_train_batch = torch.split(Y_train_shuffled, batch_size)
		Y_valid_len_batch = torch.split(Y_valid_len_shuffled, batch_size)
		sample_weights_batch = torch.split(sample_weights_shuffled, batch_size)

		for batch in zip(X_train_batch, X_valid_len_batch, Y_train_batch, Y_valid_len_batch, sample_weights_batch):
			
			optimizer.zero_grad()
			X_train_minibatch, X_valid_len_minibatch, Y_train_minibatch, Y_valid_len_minibatch, sample_weights_minibatch = [x.to(device) for x in batch]
			bos = torch.tensor([amino_dict['<bos>']] * Y_train_minibatch.shape[0],
								device=device).reshape(-1, 1)
			dec_input = torch.cat([bos, Y_train_minibatch[:, :-1]], 1)  
			
			Y_hat, _= net(X_train_minibatch, dec_input, X_valid_len_minibatch)
			l = loss(Y_hat, Y_train_minibatch, Y_valid_len_minibatch, sample_weights_minibatch, label_smoothing)
			l.backward() 
			d2l.grad_clipping(net, 1)
			num_tokens = Y_valid_len_minibatch.sum()
			optimizer.step()
			scheduler.step()
			with torch.no_grad():
				metric.add(l*num_tokens, num_tokens) 

			if current_step % 500 == 0:

				animator.add(current_step, (metric[0] / metric[1],))
				training_loss.append(metric[0] / metric[1])
				training_steps_count.append(current_step)
				print("training step {}, loss: {}".format(current_step, metric[0] / metric[1]))
				# save the print message to a a txt file in model name folder
				with open(model_name + '/print_message.txt', 'a') as f:
					f.write("training step {}, loss: {}".format(current_step, metric[0] / metric[1]) + '\n')

				metric = d2l.Accumulator(2)  

				if X_validation is not None:
					net.eval()

					with torch.no_grad():

						softmax_layer = nn.Softmax(dim=2)
						X_validation = X_validation.to(device)
						X_validation_valid_len = X_validation_valid_len.to(device)
						Y_validation = Y_validation.to(device)
						Y_validation_valid_len = Y_validation_valid_len.to(device)

						bos = torch.tensor([amino_dict['<bos>']] * Y_validation.shape[0],
											device=device).reshape(-1, 1)
						Y_true = torch.cat([bos, Y_validation[:, :-1]], 1)
						Y_true = Y_true.type(torch.int32)

						enc_outputs = net.encoder(X_validation, X_validation_valid_len)
						dec_state = net.decoder.init_state(enc_outputs, X_validation_valid_len)

						dec_X = torch.tensor([amino_dict['<bos>']] * X_validation.shape[0],
											device=device).reshape(-1, 1)

						for i in range(Y_validation.shape[1]-1):

								Y_raw, dec_state = net.decoder(dec_X, dec_state)
								Y = softmax_layer(Y_raw)

								dec_X = Y.argmax(dim=2)
								Y_pred = dec_X.type(torch.int32)

						Y_pred = torch.cat((net.decoder.seqX, dec_X), dim=1).type(torch.int32)
						
						hamming_scores = hamming_distance_training(Y_pred, Y_true)
						validation_score_hamming.append(hamming_scores.cpu().numpy().item())
						print("training step {}, hamming distance: {}".format(current_step,hamming_scores))
						# save the print message to a a txt file in model name folder
						with open(model_name + '/print_message.txt', 'a') as f:
							f.write("training step {}, hamming distance: {}".format(current_step,hamming_scores) + '\n')

						sequence_accuracy_list = [torch.equal(a, b) for a, b in zip(Y_pred, Y_true)]
						sequence_accuracy = sum(sequence_accuracy_list) / Y_true.shape[0]
						validation_score_sequence_accuracy.append(sequence_accuracy)
						print("training step {}, sequence accuracy: {}".format(current_step,sequence_accuracy))
						# save the print message to a a txt file in model name folder
						with open(model_name + '/print_message.txt', 'a') as f:
							f.write("training step {}, sequence accuracy: {}".format(current_step,sequence_accuracy) + '\n')

						enc_outputs = net.encoder(X_validation, X_validation_valid_len)
						dec_state = net.decoder.init_state(enc_outputs, X_validation_valid_len)
						bos = torch.tensor([amino_dict['<bos>']] * Y_validation.shape[0], device=device).reshape(-1, 1)
						dec_X = torch.tensor([amino_dict['<bos>']] * X_validation.shape[0],
											device=device).reshape(-1, 1)
						prob = torch.tensor([1] * Y_validation.shape[0], device=device).reshape(-1, 1)

						for i in range(Y_validation.shape[1]-1):
				
							Y_raw, dec_state = net.decoder(dec_X, dec_state)
							Y = softmax_layer(Y_raw)
							dec_X = Y_validation[:, i].reshape(-1,1)
							index = dec_X.type(torch.int64)
							prob_i = torch.gather(Y, dim=2, index = index.unsqueeze(dim=2))
							prob_i = prob_i.squeeze(dim=2).squeeze(dim=0)
							prob = torch.mul(prob, prob_i)
						
						average_log_prob = torch.log(prob).mean()
						validation_score_average_log_prob.append(average_log_prob.cpu().numpy().item())
						print("training step {}, average log prob: {}".format(current_step,average_log_prob))
						# save the print message to a a txt file in model name folder
						with open(model_name + '/print_message.txt', 'a') as f:
							f.write("training step {}, average log prob: {}".format(current_step,average_log_prob) + '\n')
						
			
					net.train()
			current_step += 1
			if current_step >= training_steps:
				break


	print(f'loss {metric[0] / metric[1]:.3f}, {metric[1] / timer.stop():.1f} tokens/sec on {str(device)}')

	plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
	x = np.array(training_steps_count)
	file_time = time.strftime('%y%b%d_%I%M%p', time.gmtime())
	plt.plot(x,training_loss,label="weighted_average_loss")
	plt.xlabel('Training steps')
	plt.ylabel('Weighted average loss')
	plt.title("transformer <{}> training loss".format(model_name))
	plt.grid()
	plt.legend()
	plt.savefig(model_name + '/{}_lossplot_{}.png'.format(model_name, file_time))

	if X_validation is not None:

		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = np.array(training_steps_count)
		plt.plot(x,validation_score_hamming,label="hamming_distance")
		plt.xlabel('Training steps')
		plt.ylabel('Hamming distance sum')
		plt.title("transformer <{}> validation hamming scores".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_hamming_{}.png'.format(model_name, file_time))

		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = np.array(training_steps_count)
		plt.plot(x,validation_score_sequence_accuracy,label="sequence_accuracy")
		plt.xlabel('Training steps')
		plt.ylabel('Sequence accuracy')
		plt.title("transformer <{}> validation sequence accuracy".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_sequence_accuracy_{}.png'.format(model_name, file_time))

		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = np.array(training_steps_count)
		plt.plot(x,validation_score_average_log_prob,label="average_log_prob")
		plt.xlabel('Training steps')
		plt.ylabel('Average log prob')
		plt.title("transformer <{}> validation average log prob".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_average_log_prob_{}.png'.format(model_name, file_time))
	

	file_name = model_name + '_' + file_time
	torch.save(net.state_dict(), model_name + '/' + file_name)
	
	np.savetxt(model_name + '/{}_training_loss_{}.csv'.format(model_name, file_time), training_loss, delimiter=',', fmt='%s')
	if X_validation is not None:
		np.savetxt(model_name + '/{}_validation_score_hamming_{}.csv'.format(model_name, file_time), validation_score_hamming, delimiter=',', fmt='%s')
		np.savetxt(model_name + '/{}_validation_score_sequence_accuracy_{}.csv'.format(model_name, file_time), validation_score_sequence_accuracy, delimiter=',', fmt='%s')
		np.savetxt(model_name + '/{}_validation_score_average_log_prob_{}.csv'.format(model_name, file_time), validation_score_average_log_prob, delimiter=',', fmt='%s')



def train_seq2seq_training_steps_DP(net, X_train, X_valid_len, Y_train, Y_valid_len, sample_weights, lr, training_steps, batch_size, label_smoothing, amino_dict, device, optimizer=None, warmup=1, model_name = 'model_demo', X_validation=None, X_validation_valid_len=None, Y_validation=None, Y_validation_valid_len=None):
	"""Train a model for sequence to sequence."""

	def init_weights(module):
		if isinstance(module, (nn.Linear)):
			nn.init.normal_(module.weight, mean=0.0, std=0.01)
		if isinstance(module, (nn.Embedding)):
			nn.init.normal_(module.weight, mean=0.0, std=1e-5)

	net.apply(init_weights)
	net.to(device)
	if optimizer is None:
		optimizer = torch.optim.Adam(net.parameters(), lr=lr, betas=(0.9, 0.98), eps = 1.0e-9)
	else:
		optimizer = optimizer
	if isinstance(warmup, int):
		scheduler = torch.optim.lr_scheduler.LinearLR(optimizer, total_iters=warmup, start_factor=1/warmup)
	else:
		scheduler = warmup
	loss = MaskedSoftmaxCELoss()
	net.train()

	animator = d2l.Animator(xlabel='training steps', ylabel='loss')
	training_loss = []
	validation_score_hamming = []
	validation_score_sequence_accuracy = []
	validation_score_average_log_prob = []

	# here we modify the training steps to be the number of epochs
	current_step = 0
	training_steps_count = []
	
	while current_step < training_steps:

		net.train()

		timer = d2l.Timer()
		metric = d2l.Accumulator(2)  

		index = torch.randperm(X_train.shape[0])
		X_train_shuffled = X_train[index]
		X_valid_len_shuffled = X_valid_len[index]
		Y_train_shuffled = Y_train[index]
		Y_valid_len_shuffled = Y_valid_len[index]
		sample_weights_shuffled = sample_weights[index]

		X_train_batch = torch.split(X_train_shuffled, batch_size)
		X_valid_len_batch = torch.split(X_valid_len_shuffled, batch_size)
		Y_train_batch = torch.split(Y_train_shuffled, batch_size)
		Y_valid_len_batch = torch.split(Y_valid_len_shuffled, batch_size)
		sample_weights_batch = torch.split(sample_weights_shuffled, batch_size)

		for batch in zip(X_train_batch, X_valid_len_batch, Y_train_batch, Y_valid_len_batch, sample_weights_batch):
			
			optimizer.zero_grad()
			X_train_minibatch, X_valid_len_minibatch, Y_train_minibatch, Y_valid_len_minibatch, sample_weights_minibatch = [x.to(device) for x in batch]
			bos = torch.tensor([amino_dict['<bos>']] * Y_train_minibatch.shape[0],
								device=device).reshape(-1, 1)
			dec_input = torch.cat([bos, Y_train_minibatch[:, :-1]], 1)  
			
			Y_hat, _= net(X_train_minibatch, dec_input, X_valid_len_minibatch)
			l = loss(Y_hat, Y_train_minibatch, Y_valid_len_minibatch, sample_weights_minibatch, label_smoothing)
			l.backward() 
			d2l.grad_clipping(net, 1)
			num_tokens = Y_valid_len_minibatch.sum()
			optimizer.step()
			scheduler.step()
			with torch.no_grad():
				metric.add(l*num_tokens, num_tokens) 

			if current_step % 500 == 0:

				animator.add(current_step, (metric[0] / metric[1],))
				training_loss.append(metric[0] / metric[1])
				training_steps_count.append(current_step)
				print("training step {}, loss: {}".format(current_step, metric[0] / metric[1]))
				# save the print message to a a txt file in model name folder
				with open(model_name + '/print_message.txt', 'a') as f:
					f.write("training step {}, loss: {}".format(current_step, metric[0] / metric[1]) + '\n')

				metric = d2l.Accumulator(2)  

				if X_validation is not None:
					net.eval()

					with torch.no_grad():

						softmax_layer = nn.Softmax(dim=2)
						X_validation = X_validation.to(device)
						X_validation_valid_len = X_validation_valid_len.to(device)
						Y_validation = Y_validation.to(device)
						Y_validation_valid_len = Y_validation_valid_len.to(device)

						bos = torch.tensor([amino_dict['<bos>']] * Y_validation.shape[0],
											device=device).reshape(-1, 1)
						Y_true = torch.cat([bos, Y_validation[:, :-1]], 1)
						Y_true = Y_true.type(torch.int32)

						enc_outputs = net.module.encoder(X_validation, X_validation_valid_len)
						dec_state = net.module.decoder.init_state(enc_outputs, X_validation_valid_len)

						dec_X = torch.tensor([amino_dict['<bos>']] * X_validation.shape[0],
											device=device).reshape(-1, 1)

						for i in range(Y_validation.shape[1]-1):

								Y_raw, dec_state = net.module.decoder(dec_X, dec_state)
								Y = softmax_layer(Y_raw)

								dec_X = Y.argmax(dim=2)
								Y_pred = dec_X.type(torch.int32)

						Y_pred = torch.cat((net.module.decoder.seqX, dec_X), dim=1).type(torch.int32)
						
						hamming_scores = hamming_distance_training(Y_pred, Y_true)
						validation_score_hamming.append(hamming_scores.cpu().numpy().item())
						print("training step {}, hamming distance: {}".format(current_step,hamming_scores))
						# save the print message to a a txt file in model name folder
						with open(model_name + '/print_message.txt', 'a') as f:
							f.write("training step {}, hamming distance: {}".format(current_step,hamming_scores) + '\n')

						sequence_accuracy_list = [torch.equal(a, b) for a, b in zip(Y_pred, Y_true)]
						sequence_accuracy = sum(sequence_accuracy_list) / Y_true.shape[0]
						validation_score_sequence_accuracy.append(sequence_accuracy)
						print("training step {}, sequence accuracy: {}".format(current_step,sequence_accuracy))
						# save the print message to a a txt file in model name folder
						with open(model_name + '/print_message.txt', 'a') as f:
							f.write("training step {}, sequence accuracy: {}".format(current_step,sequence_accuracy) + '\n')

						enc_outputs = net.module.encoder(X_validation, X_validation_valid_len)
						dec_state = net.module.decoder.init_state(enc_outputs, X_validation_valid_len)
						bos = torch.tensor([amino_dict['<bos>']] * Y_validation.shape[0], device=device).reshape(-1, 1)
						dec_X = torch.tensor([amino_dict['<bos>']] * X_validation.shape[0],
											device=device).reshape(-1, 1)
						prob = torch.tensor([1] * Y_validation.shape[0], device=device).reshape(-1, 1)

						for i in range(Y_validation.shape[1]-1):
				
							Y_raw, dec_state = net.module.decoder(dec_X, dec_state)
							Y = softmax_layer(Y_raw)
							dec_X = Y_validation[:, i].reshape(-1,1)
							index = dec_X.type(torch.int64)
							prob_i = torch.gather(Y, dim=2, index = index.unsqueeze(dim=2))
							prob_i = prob_i.squeeze(dim=2).squeeze(dim=0)
							prob = torch.mul(prob, prob_i)
						
						average_log_prob = torch.log(prob).mean()
						validation_score_average_log_prob.append(average_log_prob.cpu().numpy().item())
						print("training step {}, average log prob: {}".format(current_step,average_log_prob))
						# save the print message to a a txt file in model name folder
						with open(model_name + '/print_message.txt', 'a') as f:
							f.write("training step {}, average log prob: {}".format(current_step,average_log_prob) + '\n')
						
					
					net.train()
			current_step += 1
			if current_step >= training_steps:
				break


	print(f'loss {metric[0] / metric[1]:.3f}, {metric[1] / timer.stop():.1f} tokens/sec on {str(device)}')

	plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
	x = np.array(training_steps_count)
	file_time = time.strftime('%y%b%d_%I%M%p', time.gmtime())
	plt.plot(x,training_loss,label="weighted_average_loss")
	plt.xlabel('Training steps')
	plt.ylabel('Weighted average loss')
	plt.title("transformer <{}> training loss".format(model_name))
	plt.grid()
	plt.legend()
	plt.savefig(model_name + '/{}_lossplot_{}.png'.format(model_name, file_time))

	if X_validation is not None:

		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = np.array(training_steps_count)
		plt.plot(x,validation_score_hamming,label="hamming_distance")
		plt.xlabel('Training steps')
		plt.ylabel('Hamming distance sum')
		plt.title("transformer <{}> validation hamming scores".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_hamming_{}.png'.format(model_name, file_time))

		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = np.array(training_steps_count)
		plt.plot(x,validation_score_sequence_accuracy,label="sequence_accuracy")
		plt.xlabel('Training steps')
		plt.ylabel('Sequence accuracy')
		plt.title("transformer <{}> validation sequence accuracy".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_sequence_accuracy_{}.png'.format(model_name, file_time))

		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = np.array(training_steps_count)
		plt.plot(x,validation_score_average_log_prob,label="average_log_prob")
		plt.xlabel('Training steps')
		plt.ylabel('Average log prob')
		plt.title("transformer <{}> validation average log prob".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_average_log_prob_{}.png'.format(model_name, file_time))
	

	file_name = model_name + '_' + file_time
	torch.save(net.state_dict(), model_name + '/' + file_name)
	
	np.savetxt(model_name + '/{}_training_loss_{}.csv'.format(model_name, file_time), training_loss, delimiter=',', fmt='%s')
	if X_validation is not None:
		np.savetxt(model_name + '/{}_validation_score_hamming_{}.csv'.format(model_name, file_time), validation_score_hamming, delimiter=',', fmt='%s')
		np.savetxt(model_name + '/{}_validation_score_sequence_accuracy_{}.csv'.format(model_name, file_time), validation_score_sequence_accuracy, delimiter=',', fmt='%s')
		np.savetxt(model_name + '/{}_validation_score_average_log_prob_{}.csv'.format(model_name, file_time), validation_score_average_log_prob, delimiter=',', fmt='%s')



def train_prott5_DP(net, X_train, X_valid_len, Y_train, Y_valid_len, sample_weights, lr, training_steps, batch_size, label_smoothing, amino_dict, device, optimizer=None, warmup=1, model_name = 'model_demo', X_validation=None, X_validation_valid_len=None, Y_validation=None, Y_validation_valid_len=None):
	"""Train a model for sequence to sequence."""

	def init_weights(module):
		if isinstance(module, (nn.Embedding)):
			nn.init.normal_(module.weight, mean=0.0, std=1e-5)
	net.module.apply(init_weights)

	nn.init.normal_(net.module.dense.weight, 0, 0.01)

	net.to(device)
	if optimizer is None:
		optimizer = torch.optim.Adam(net.parameters(), lr=lr, betas=(0.9, 0.98), eps = 1.0e-9)
	else:
		optimizer = optimizer
	if isinstance(warmup, int):
		scheduler = torch.optim.lr_scheduler.LinearLR(optimizer, total_iters=warmup, start_factor=1/warmup)
	else:
		scheduler = warmup
	loss = MaskedSoftmaxCELoss()
	net.train()

	animator = d2l.Animator(xlabel='training steps', ylabel='loss')
	training_loss = []
	validation_score_hamming = []
	validation_score_sequence_accuracy = []
	validation_score_average_log_prob = []

	# here we modify the training steps to be the number of epochs
	current_step = 0
	training_steps_count = []

	while current_step < training_steps:

		net.train()

		timer = d2l.Timer()
		metric = d2l.Accumulator(2)  

		index = torch.randperm(X_train.shape[0])
		X_train_shuffled = X_train[index]
		X_valid_len_shuffled = X_valid_len[index]
		Y_train_shuffled = Y_train[index]
		Y_valid_len_shuffled = Y_valid_len[index]
		sample_weights_shuffled = sample_weights[index]

		X_train_batch = torch.split(X_train_shuffled, batch_size)
		X_valid_len_batch = torch.split(X_valid_len_shuffled, batch_size)
		Y_train_batch = torch.split(Y_train_shuffled, batch_size)
		Y_valid_len_batch = torch.split(Y_valid_len_shuffled, batch_size)
		sample_weights_batch = torch.split(sample_weights_shuffled, batch_size)

		for batch in zip(X_train_batch, X_valid_len_batch, Y_train_batch, Y_valid_len_batch, sample_weights_batch):
			
			optimizer.zero_grad()
			X_train_minibatch, X_valid_len_minibatch, Y_train_minibatch, Y_valid_len_minibatch, sample_weights_minibatch = [x.to(device) for x in batch]
			bos = torch.tensor([amino_dict['<bos>']] * Y_train_minibatch.shape[0],
								device=device).reshape(-1, 1)
			dec_input = torch.cat([bos, Y_train_minibatch[:, :-1]], 1)  
			
			Y_hat = net(X_train_minibatch, X_valid_len_minibatch, dec_input)
			l = loss(Y_hat, Y_train_minibatch, Y_valid_len_minibatch, sample_weights_minibatch, label_smoothing)
			l.backward() 
			d2l.grad_clipping(net, 1)
			num_tokens = Y_valid_len_minibatch.sum()
			optimizer.step()
			scheduler.step()
			with torch.no_grad():
				metric.add(l*num_tokens, num_tokens) 

			if current_step % 500 == 0:

				animator.add(current_step, (metric[0] / metric[1],))
				training_loss.append(metric[0] / metric[1])
				training_steps_count.append(current_step)
				print("training step {}, loss: {}".format(current_step, metric[0] / metric[1]))
				# save the print message to a a txt file in model name folder
				with open(model_name + '/print_message.txt', 'a') as f:
					f.write("training step {}, loss: {}".format(current_step, metric[0] / metric[1]) + '\n')

				metric = d2l.Accumulator(2)  

				if X_validation is not None:
					net.eval()

					with torch.no_grad():

						softmax_layer = nn.Softmax(dim=2)
						X_validation = X_validation.to(device)
						X_validation_valid_len = X_validation_valid_len.to(device)
						Y_validation = Y_validation.to(device)
						Y_validation_valid_len = Y_validation_valid_len.to(device)

						bos = torch.tensor([amino_dict['<bos>']] * Y_validation.shape[0],
											device=device).reshape(-1, 1)
						Y_true = torch.cat([bos, Y_validation[:, :-1]], 1)
						Y_true = Y_true.type(torch.int32)

						dec_X = torch.tensor([amino_dict['<bos>']] * X_validation.shape[0],
											device=device).reshape(-1, 1)
						# copy dec_X to a new tensor
						Y_pred = dec_X.clone()

						for i in range(Y_validation.shape[1]):

								Y_raw = net(X_validation, X_validation_valid_len, dec_X)
								Y = softmax_layer(Y_raw)

								dec_X = Y.argmax(dim=2)
								Y_pred = dec_X.type(torch.int32)

								# append the new prediction to the end of the tensor
								Y_pred = torch.cat((Y_pred, dec_X), dim=1)

						Y_pred = torch.cat((net.module.decoder.seqX, dec_X), dim=1).type(torch.int32)
						
						hamming_scores = hamming_distance_training(Y_pred, Y_true)
						validation_score_hamming.append(hamming_scores.cpu().numpy().item())
						print("training step {}, hamming distance: {}".format(current_step,hamming_scores))
						# save the print message to a a txt file in model name folder
						with open(model_name + '/print_message.txt', 'a') as f:
							f.write("training step {}, hamming distance: {}".format(current_step,hamming_scores) + '\n')

						sequence_accuracy_list = [torch.equal(a, b) for a, b in zip(Y_pred, Y_true)]
						sequence_accuracy = sum(sequence_accuracy_list) / Y_true.shape[0]
						validation_score_sequence_accuracy.append(sequence_accuracy)
						print("training step {}, sequence accuracy: {}".format(current_step,sequence_accuracy))
						# save the print message to a a txt file in model name folder
						with open(model_name + '/print_message.txt', 'a') as f:
							f.write("training step {}, sequence accuracy: {}".format(current_step,sequence_accuracy) + '\n')

						bos = torch.tensor([amino_dict['<bos>']] * Y_validation.shape[0], device=device).reshape(-1, 1)
						dec_X = torch.tensor([amino_dict['<bos>']] * X_validation.shape[0],
											device=device).reshape(-1, 1)
						prob = torch.tensor([1] * Y_validation.shape[0], device=device).reshape(-1, 1)

						for i in range(Y_validation.shape[1]-1):
				
							Y_raw = net(X_validation, X_validation_valid_len, dec_X)
							Y = softmax_layer(Y_raw)
							dec_X = Y_validation[:, i].reshape(-1,1)
							index = dec_X.type(torch.int64)
							prob_i = torch.gather(Y, dim=2, index = index.unsqueeze(dim=2))
							prob_i = prob_i.squeeze(dim=2).squeeze(dim=0)
							prob = torch.mul(prob, prob_i)
						
						average_log_prob = torch.log(prob).mean()
						validation_score_average_log_prob.append(average_log_prob.cpu().numpy().item())
						print("training step {}, average log prob: {}".format(current_step,average_log_prob))
						# save the print message to a a txt file in model name folder
						with open(model_name + '/print_message.txt', 'a') as f:
							f.write("training step {}, average log prob: {}".format(current_step,average_log_prob) + '\n')
						
					
					net.train()
			current_step += 1
			if current_step >= training_steps:
				break


	print(f'loss {metric[0] / metric[1]:.3f}, {metric[1] / timer.stop():.1f} tokens/sec on {str(device)}')

	plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
	x = np.array(training_steps_count)
	file_time = time.strftime('%y%b%d_%I%M%p', time.gmtime())
	plt.plot(x,training_loss,label="weighted_average_loss")
	plt.xlabel('Training steps')
	plt.ylabel('Weighted average loss')
	plt.title("transformer <{}> training loss".format(model_name))
	plt.grid()
	plt.legend()
	plt.savefig(model_name + '/{}_lossplot_{}.png'.format(model_name, file_time))

	if X_validation is not None:

		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = np.array(training_steps_count)
		plt.plot(x,validation_score_hamming,label="hamming_distance")
		plt.xlabel('Training steps')
		plt.ylabel('Hamming distance sum')
		plt.title("transformer <{}> validation hamming scores".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_hamming_{}.png'.format(model_name, file_time))

		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = np.array(training_steps_count)
		plt.plot(x,validation_score_sequence_accuracy,label="sequence_accuracy")
		plt.xlabel('Training steps')
		plt.ylabel('Sequence accuracy')
		plt.title("transformer <{}> validation sequence accuracy".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_sequence_accuracy_{}.png'.format(model_name, file_time))

		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = np.array(training_steps_count)
		plt.plot(x,validation_score_average_log_prob,label="average_log_prob")
		plt.xlabel('Training steps')
		plt.ylabel('Average log prob')
		plt.title("transformer <{}> validation average log prob".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_average_log_prob_{}.png'.format(model_name, file_time))


	file_name = model_name + '_' + file_time
	torch.save(net.state_dict(), model_name + '/' + file_name)

	np.savetxt(model_name + '/{}_training_loss_{}.csv'.format(model_name, file_time), training_loss, delimiter=',', fmt='%s')
	if X_validation is not None:
		np.savetxt(model_name + '/{}_validation_score_hamming_{}.csv'.format(model_name, file_time), validation_score_hamming, delimiter=',', fmt='%s')
		np.savetxt(model_name + '/{}_validation_score_sequence_accuracy_{}.csv'.format(model_name, file_time), validation_score_sequence_accuracy, delimiter=',', fmt='%s')
		np.savetxt(model_name + '/{}_validation_score_average_log_prob_{}.csv'.format(model_name, file_time), validation_score_average_log_prob, delimiter=',', fmt='%s')



def train_seq2seq_training_steps_DP_checkpoint(net, X_train, X_valid_len, Y_train, Y_valid_len, sample_weights, lr, training_steps, batch_size, label_smoothing, amino_dict, device, optimizer=None, warmup=1, model_name = 'model_demo', X_validation=None, X_validation_valid_len=None, Y_validation=None, Y_validation_valid_len=None, 
												resume = False, current_step = 0, target_step = 0):
	"""Train a model for sequence to sequence."""

	if resume == False:

		def init_weights(module):
			if isinstance(module, (nn.Linear)):
				nn.init.normal_(module.weight, mean=0.0, std=0.01)
			if isinstance(module, (nn.Embedding)):
				nn.init.normal_(module.weight, mean=0.0, std=1e-5)

		net.apply(init_weights)
		net.to(device)
		if optimizer is None:
			optimizer = torch.optim.Adam(net.parameters(), lr=lr, betas=(0.9, 0.98), eps = 1.0e-9)
		else:
			optimizer = optimizer
		if isinstance(warmup, int):
			scheduler = torch.optim.lr_scheduler.LinearLR(optimizer, total_iters=warmup, start_factor=1/warmup)
		else:
			scheduler = warmup
		loss = MaskedSoftmaxCELoss()

		net.to(device)
		net.train()

		animator = d2l.Animator(xlabel='training steps', ylabel='loss')
		training_loss = []
		validation_score_hamming = []
		validation_score_sequence_accuracy = []
		validation_score_average_log_prob = []
		validation_score_loss =[]

		# here we modify the training steps to be the number of epochs
		current_step = 0
		training_steps_count = []
	
	else:
		optimizer = optimizer
		scheduler = warmup
		loss = MaskedSoftmaxCELoss()
		net.to(device)
		net.train()

		animator = d2l.Animator(xlabel='training steps', ylabel='loss')
		training_loss = []
		validation_score_hamming = []
		validation_score_sequence_accuracy = []
		validation_score_average_log_prob = []

		current_step = current_step
		training_steps_count = []

	while current_step < training_steps:

		net.train()

		timer = d2l.Timer()
		metric = d2l.Accumulator(2)  

		index = torch.randperm(X_train.shape[0])
		X_train_shuffled = X_train[index]
		X_valid_len_shuffled = X_valid_len[index]
		Y_train_shuffled = Y_train[index]
		Y_valid_len_shuffled = Y_valid_len[index]
		sample_weights_shuffled = sample_weights[index]

		X_train_batch = torch.split(X_train_shuffled, batch_size)
		X_valid_len_batch = torch.split(X_valid_len_shuffled, batch_size)
		Y_train_batch = torch.split(Y_train_shuffled, batch_size)
		Y_valid_len_batch = torch.split(Y_valid_len_shuffled, batch_size)
		sample_weights_batch = torch.split(sample_weights_shuffled, batch_size)

		for batch in zip(X_train_batch, X_valid_len_batch, Y_train_batch, Y_valid_len_batch, sample_weights_batch):
			
			optimizer.zero_grad()
			X_train_minibatch, X_valid_len_minibatch, Y_train_minibatch, Y_valid_len_minibatch, sample_weights_minibatch = [x.to(device) for x in batch]
			bos = torch.tensor([amino_dict['<bos>']] * Y_train_minibatch.shape[0],
								device=device).reshape(-1, 1)
			dec_input = torch.cat([bos, Y_train_minibatch[:, :-1]], 1)  
			
			Y_hat, _= net(X_train_minibatch, dec_input, X_valid_len_minibatch)
			l = loss(Y_hat, Y_train_minibatch, Y_valid_len_minibatch, sample_weights_minibatch, label_smoothing)
			l.backward() 
			d2l.grad_clipping(net, 1)
			num_tokens = Y_valid_len_minibatch.sum()
			optimizer.step()
			scheduler.step()
			with torch.no_grad():
				metric.add(l*num_tokens, num_tokens) 

			if current_step % 500 == 0:

				animator.add(current_step, (metric[0] / metric[1],))
				training_loss.append(metric[0] / metric[1])
				training_steps_count.append(current_step)
				print("training step {}, loss: {}".format(current_step, metric[0] / metric[1]))
				# save the print message to a a txt file in model name folder
				with open(model_name + '/print_message.txt', 'a') as f:
					f.write("training step {}, loss: {}".format(current_step, metric[0] / metric[1]) + '\n')

				metric = d2l.Accumulator(2)

				# print learning rate
				print("learning rate: {}".format(scheduler.get_last_lr()))
				# save the print message to a a txt file in model name folder
				with open(model_name + '/print_message.txt', 'a') as f:
					f.write("learning rate: {}".format(scheduler.get_last_lr()) + '\n')

				if X_validation is not None:
					net.eval()

					with torch.no_grad():

						softmax_layer = nn.Softmax(dim=2)
						X_validation = X_validation.to(device)
						X_validation_valid_len = X_validation_valid_len.to(device)
						Y_validation = Y_validation.to(device)
						Y_validation_valid_len = Y_validation_valid_len.to(device)

						bos = torch.tensor([amino_dict['<bos>']] * Y_validation.shape[0],
											device=device).reshape(-1, 1)
						Y_true = torch.cat([bos, Y_validation[:, :-1]], 1)
						Y_true = Y_true.type(torch.int32)

						enc_outputs = net.module.encoder(X_validation, X_validation_valid_len)
						dec_state = net.module.decoder.init_state(enc_outputs, X_validation_valid_len)

						dec_X = torch.tensor([amino_dict['<bos>']] * X_validation.shape[0],
											device=device).reshape(-1, 1)

						for i in range(Y_validation.shape[1]-1):

								Y_raw, dec_state = net.module.decoder(dec_X, dec_state)
								Y = softmax_layer(Y_raw)

								dec_X = Y.argmax(dim=2)
								Y_pred = dec_X.type(torch.int32)

						Y_pred = torch.cat((net.module.decoder.seqX, dec_X), dim=1).type(torch.int32)
						
						hamming_scores = hamming_distance_training(Y_pred, Y_true)
						validation_score_hamming.append(hamming_scores.cpu().numpy().item())
						print("training step {}, hamming distance: {}".format(current_step,hamming_scores))
						# save the print message to a a txt file in model name folder
						with open(model_name + '/print_message.txt', 'a') as f:
							f.write("training step {}, hamming distance: {}".format(current_step,hamming_scores) + '\n')

						sequence_accuracy_list = [torch.equal(a, b) for a, b in zip(Y_pred, Y_true)]
						sequence_accuracy = sum(sequence_accuracy_list) / Y_true.shape[0]
						validation_score_sequence_accuracy.append(sequence_accuracy)
						print("training step {}, sequence accuracy: {}".format(current_step,sequence_accuracy))
						# save the print message to a a txt file in model name folder
						with open(model_name + '/print_message.txt', 'a') as f:
							f.write("training step {}, sequence accuracy: {}".format(current_step,sequence_accuracy) + '\n')

						# calculate the validation loss on the whole validation set
						net.train()
						validation_loss = 0
						for v_batch in zip(torch.split(X_validation, batch_size),
										torch.split(X_validation_valid_len, batch_size),
										torch.split(Y_validation, batch_size),
										torch.split(Y_validation_valid_len, batch_size)):
							X_validation_minibatch, X_validation_valid_len_minibatch, Y_validation_minibatch, Y_validation_valid_len_minibatch = [x.to(device) for x in v_batch]
							bos = torch.tensor([amino_dict['<bos>']] * Y_validation_minibatch.shape[0],
												device=device).reshape(-1, 1)
							dec_input = torch.cat([bos, Y_validation_minibatch[:, :-1]], 1)
							Y_hat_validation, _ = net(X_validation_minibatch, dec_input, X_validation_valid_len_minibatch)
							l_validation = loss(Y_hat_validation, Y_validation_minibatch, Y_validation_valid_len_minibatch, None, label_smoothing)
							validation_loss += l_validation.item() * Y_validation_valid_len_minibatch.sum().item()

						validation_loss /= Y_validation_valid_len.sum().item()
						validation_score_loss.append(validation_loss)
						print("training step {}, validation loss: {}".format(current_step,validation_loss))
						# save the print message to a a txt file in model name folder
						with open(model_name + '/print_message.txt', 'a') as f:
							f.write("training step {}, validation loss: {}".format(current_step,validation_loss) + '\n')
						net.eval()

						enc_outputs = net.module.encoder(X_validation, X_validation_valid_len)
						dec_state = net.module.decoder.init_state(enc_outputs, X_validation_valid_len)
						bos = torch.tensor([amino_dict['<bos>']] * Y_validation.shape[0], device=device).reshape(-1, 1)
						dec_X = torch.tensor([amino_dict['<bos>']] * X_validation.shape[0],
											device=device).reshape(-1, 1)
						prob = torch.tensor([1] * Y_validation.shape[0], device=device).reshape(-1, 1)

						for i in range(Y_validation.shape[1]-1):
				
							Y_raw, dec_state = net.module.decoder(dec_X, dec_state)
							Y = softmax_layer(Y_raw)
							dec_X = Y_validation[:, i].reshape(-1,1)
							index = dec_X.type(torch.int64)
							prob_i = torch.gather(Y, dim=2, index = index.unsqueeze(dim=2))
							prob_i = prob_i.squeeze(dim=2).squeeze(dim=0)
							prob = torch.mul(prob, prob_i)
						
						average_log_prob = torch.log(prob).mean()
						validation_score_average_log_prob.append(average_log_prob.cpu().numpy().item())
						print("training step {}, average log prob: {}".format(current_step,average_log_prob))
						# save the print message to a a txt file in model name folder
						with open(model_name + '/print_message.txt', 'a') as f:
							f.write("training step {}, average log prob: {}".format(current_step,average_log_prob) + '\n')
						
					
					net.train()
			current_step += 1
			if current_step >= training_steps:
				break


	print(f'loss {metric[0] / metric[1]:.3f}, {metric[1] / timer.stop():.1f} tokens/sec on {str(device)}')

	plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
	x = np.array(training_steps_count)
	file_time = time.strftime('%y%b%d_%I%M%p', time.gmtime())
	plt.plot(x,training_loss,label="weighted_average_loss")
	plt.xlabel('Training steps')
	plt.ylabel('Weighted average loss')
	plt.title("transformer <{}> training loss".format(model_name))
	plt.grid()
	plt.legend()
	plt.savefig(model_name + '/{}_lossplot.png'.format(model_name))

	if X_validation is not None:

		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = np.array(training_steps_count)
		plt.plot(x,validation_score_hamming,label="hamming_distance")
		plt.xlabel('Training steps')
		plt.ylabel('Hamming distance sum')
		plt.title("transformer <{}> validation hamming scores".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_hamming.png'.format(model_name))

		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = np.array(training_steps_count)
		plt.plot(x,validation_score_sequence_accuracy,label="sequence_accuracy")
		plt.xlabel('Training steps')
		plt.ylabel('Sequence accuracy')
		plt.title("transformer <{}> validation sequence accuracy".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_sequence_accuracy.png'.format(model_name))

		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = np.array(training_steps_count)
		plt.plot(x,validation_score_loss,label="validation_loss")
		plt.xlabel('Training steps')
		plt.ylabel('Validation loss')
		plt.title("transformer <{}> validation loss".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_loss.png'.format(model_name))

		plt.figure(figsize=(12, 8), facecolor=(1, 1, 1))
		x = np.array(training_steps_count)
		plt.plot(x,validation_score_average_log_prob,label="average_log_prob")
		plt.xlabel('Training steps')
		plt.ylabel('Average log prob')
		plt.title("transformer <{}> validation average log prob".format(model_name))
		plt.grid()
		plt.legend()
		plt.savefig(model_name + '/{}_validationplot_average_log_prob.png'.format(model_name))
	
	# save a checkpoint
	model_state_dict = net.state_dict()
	# remove the module prefix
	from collections import OrderedDict
	new_state_dict = OrderedDict()
	for k, v in model_state_dict.items():
		name = k[7:] #remove 'module'
		new_state_dict[name] = v

	checkpoint = {
		'model': new_state_dict,
		'optimizer': optimizer.state_dict(),
		'scheduler': scheduler.state_dict(),
		'current_step': current_step,
	}

	torch.save(checkpoint, model_name + '/{}_checkpoint_{}.pt'.format(model_name, file_time))

	file_name = model_name + '_' + file_time
	torch.save(net.state_dict(), model_name + '/' + file_name)

	
	np.savetxt(model_name + '/{}_training_loss.csv'.format(model_name), training_loss, delimiter=',', fmt='%s')
	if X_validation is not None:
		np.savetxt(model_name + '/{}_validation_score_hamming.csv'.format(model_name), validation_score_hamming, delimiter=',', fmt='%s')
		np.savetxt(model_name + '/{}_validation_score_sequence_accuracy.csv'.format(model_name), validation_score_sequence_accuracy, delimiter=',', fmt='%s')
		np.savetxt(model_name + '/{}_validation_score_loss.csv'.format(model_name), validation_score_loss, delimiter=',', fmt='%s')
		np.savetxt(model_name + '/{}_validation_score_average_log_prob.csv'.format(model_name), validation_score_average_log_prob, delimiter=',', fmt='%s')
