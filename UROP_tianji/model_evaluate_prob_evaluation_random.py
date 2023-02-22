from TransformerBeta import *
import torch.nn as nn
import numpy as np
import os 
import matplotlib.pyplot as plt

# 0.0 load the validation data
validation_dict = np.load('validation_l8_anti/validation_dict_fold0.npy', allow_pickle=True) # <------------------------------------------------------------------------change
validation_dict = validation_dict.tolist()
validation_list = []
for target, value_dict in validation_dict.items():
    for comp, count in value_dict.items():
        validation_list.append([target, comp, count])
validation_array = np.array(validation_list)

# HPC setting
num_resources = 5
counter = 0
array_index = int(os.environ['PBS_ARRAY_INDEX'])
array_index_list = list(range(1, num_resources+1))

fold_list = ["fold0", "fold1", "fold2", "fold3", "fold4"] # <------------------------------------------------------------------------change
for fold in fold_list:
	if array_index_list[counter%num_resources] == array_index:
		# only predict 1/5 of the validation data, thus divide the validation data into 5 folds
		# fold 0
		if fold == "fold0":
			validation_array = validation_array[0::5]
		# fold 1
		elif fold == "fold1":
			validation_array = validation_array[1::5]
		# fold 2
		elif fold == "fold2":
			validation_array = validation_array[2::5]
		# fold 3
		elif fold == "fold3":
			validation_array = validation_array[3::5]
		# fold 4
		elif fold == "fold4":
			validation_array = validation_array[4::5]
		else:
			print("Error: fold not found")
			break

		X_validation = validation_array[:, 0]
		Y_validation = validation_array[:, 1]
		count_array = validation_array[:, 2]



		# 0.1 load the model
		model_path = "AF2_transformer_l12h768_validation_anti_l7_100ksteps_22Dec11_1202AM" # <------------------------------------------------------------------------change

		query_size, key_size, value_size, num_hiddens = 512, 512, 512, 512
		num_layers, dropout = 6, 0.1
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

		state_dict = torch.load(model_path,map_location=('cpu'))
		from collections import OrderedDict
		new_state_dict = OrderedDict()
		for k, v in state_dict.items():
			name = k[7:] #remove 'module'
			new_state_dict[name] = v

		model_standard.load_state_dict(new_state_dict)

		model_use = model_standard # <------------------------------------------------------------------------change
		prediction_length = 8 # <------------------------------------------------------------------------change

		# 0.2 create a txt file to record the results
		if not os.path.exists('model_evaluation'):
			os.mkdir('model_evaluation')



		# 0.3 basic statistics of the validation data and model
		print('Standard model: total number of parameters: {}'.format(model_standard_total_params))
		print('Standard model: total number of trainable parameters: {}'.format(model_standard_total_trainable_params))
		print('Standard model: total number of validation data: {}'.format(len(X_validation)))
		with open('model_evaluation/prob_evaluation_random_{}.txt'.format(fold), 'w') as f:
			f.write('Standard model: total number of parameters: {}\n'.format(model_standard_total_params))
			f.write('Standard model: total number of trainable parameters: {}\n'.format(model_standard_total_trainable_params))
			f.write('Standard model: total number of validation data: {}\n'.format(len(X_validation)))



		# 1. sequence random validation average log probability
		np.random.seed(0)
		num_random_sequences = X_validation.shape[0]
		amino_list = list(amino_dict.keys())[4:]
		random_labels = ["".join(list(np.random.choice(amino_list, prediction_length))) for i in range(num_random_sequences)]

		num_minibatch = 100 # <------------------------------------------------------------------------change
		peptide_eval_random = evaluate_minibatch(model_use, validation_array[:, 0], random_labels, amino_dict, prediction_length + 2, device, num_minibatch=num_minibatch, print_info=True, output_file='model_evaluation/prob_evaluation_random_{}.txt'.format(fold))
		print("minibatch evaluatiion number: \n", len(peptide_eval_random))
		print("minibatch evaluatiion example data: \n", peptide_eval_random[:10])
		y = peptide_eval_random[:, 2].astype(np.float32)
		average_log_prob = np.mean(np.log(y))
		median_prob = np.median(y)
		print("Minibatch evaluatiion average log probability: \n", average_log_prob)
		print("Minibatch evaluatiion median probability: \n", median_prob)
		with open('model_evaluation/prob_evaluation_random_{}.txt'.format(fold), 'a') as f:
			f.write("minibatch evaluatiion number: \n {}\n".format(len(peptide_eval_random)))
			f.write("minibatch evaluatiion example data: \n {}\n".format(peptide_eval_random[:10]))
			f.write("Minibatch evaluatiion average log probability: \n {}\n".format(average_log_prob))
			f.write("Minibatch evaluatiion median probability: \n {}\n".format(median_prob))
		# 1.1 save the results
		np.save('model_evaluation/peptide_eval_random_{}.npy'.format(fold), peptide_eval_random)

	counter += 1