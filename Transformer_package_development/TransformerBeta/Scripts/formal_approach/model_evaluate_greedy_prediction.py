from TransformerBeta import *
import torch.nn as nn
import numpy as np
import os 
import matplotlib.pyplot as plt

# 0.0 load the validation data
validation_dict = np.load('validation_l7_anti/validation_dict_fold0.npy', allow_pickle=True) # <------------------------------------------------------------------------change
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
		model_path = "AF2_transformer_l12h1024_validation_anti_l7_22Nov07_0904AM" # <------------------------------------------------------------------------change

		query_size, key_size, value_size, num_hiddens = 1024, 1024, 1024, 1024
		num_layers, dropout = 12, 0.1
		lr, training_steps, batch_size, label_smoothing = 0.0002, 100000, 4096, 0.1
		ffn_num_input, ffn_num_hiddens, num_heads = 1024, 4096, 8

		norm_shape = [1024] # 1024 corresponds to the dim of such number to normalize
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
		prediction_length = 7 # <------------------------------------------------------------------------change

		# 0.2 create a txt file to record the results
		if not os.path.exists('model_evaluation'):
			os.mkdir('model_evaluation')



		# 0.3 basic statistics of the validation data and model
		print('Standard model: total number of parameters: {}'.format(model_standard_total_params))
		print('Standard model: total number of trainable parameters: {}'.format(model_standard_total_trainable_params))
		print('Standard model: total number of validation data: {}'.format(len(X_validation)))
		with open('model_evaluation/greedy_prediction_{}.txt'.format(fold), 'w') as f:
			f.write('Standard model: total number of parameters: {}\n'.format(model_standard_total_params))
			f.write('Standard model: total number of trainable parameters: {}\n'.format(model_standard_total_trainable_params))
			f.write('Standard model: total number of validation data: {}\n'.format(len(X_validation)))



		# 1. sequence validation accuracy
		num_minibatch = 100 # <------------------------------------------------------------------------change
		peptide_pred = predict_greedy_minibatch(model_use, validation_array[:, 0], amino_dict, prediction_length + 2, device, num_minibatch=num_minibatch, print_info=True, output_file='model_evaluation/greedy_prediction_{}.txt'.format(fold))
		print("Greedy minibatch prediction number: \n", len(peptide_pred))
		print("Greedy minibatch prediction example data: \n", peptide_pred[:10])

		correct = 0
		for pred, truth in zip(peptide_pred[:, 1], validation_array[:, 1]):
			if pred == truth:
				correct += 1
		print("Number of correct prediction: ", correct)
		print("Number of total prediction: ", len(peptide_pred))
		print("Greedy minibatch prediction accuracy: \n", correct / len(peptide_pred))
		with open('model_evaluation/greedy_prediction_{}.txt'.format(fold), 'a') as f:
			f.write("Greedy minibatch prediction number: \n {}\n".format(len(peptide_pred)))
			f.write("Greedy minibatch prediction example data: \n {}\n".format(peptide_pred[:10]))
			f.write("Number of correct prediction: {}\n".format(correct))
			f.write("Number of total prediction: {}\n".format(len(peptide_pred)))
			f.write("Greedy minibatch prediction accuracy: \n {}\n".format(correct / len(peptide_pred)))
		# 1.1 save the prediction result
		np.save('model_evaluation/peptide_pred_{}.npy'.format(fold), peptide_pred)

	counter += 1
