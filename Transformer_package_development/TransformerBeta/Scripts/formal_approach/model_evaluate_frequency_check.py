import numpy as np
import matplotlib.pyplot as plt

# load validation data
validation_dict = np.load('validation_l7_anti/validation_dict_fold0.npy', allow_pickle=True) # <------------------------------------------------------------------------change
validation_dict = validation_dict.tolist()
validation_list = []
for target, value_dict in validation_dict.items():
    for comp, count in value_dict.items():
        validation_list.append([target, comp, count])
validation_array = np.array(validation_list)

# load training data
train_dict = np.load('train_l7_anti/train_dict_fold0.npy', allow_pickle=True)
train_dict = train_dict.tolist()
train_list = []
for target, value_dict in train_dict.items():
    for comp, count in value_dict.items():
        train_list.append([target, comp, count])
train_array = np.array(train_list)

# load peptide_pred
peptide_pred = np.load('model_evaluation/peptide_pred.npy', allow_pickle=True)

amino_dict_freq = {'A': 0,
'C': 1,
'D': 2,
'E': 3,
'F': 4,
'G': 5,
'H': 6,
'I': 7,
'K': 8,
'L': 9,
'M': 10,
'N': 11,
'P': 12,
'Q': 13,
'R': 14,
'S': 15,
'T': 16,
'V': 17,
'W': 18,
'Y': 19}

# 1. calculate amino acid frequency in each position for validation data

# create a 2 by 2 array to store number of amino acids in each position
amino_count = np.zeros((len(validation_array[0][1]), len(amino_dict_freq)))

# count amino acids in each position
for i in range(len(validation_array)):
	for j in range(len(validation_array[i][1])):
		amino_count[j][amino_dict_freq[validation_array[i][1][j]]] += 1
        
# calculate amino acid frequency in each position
amino_freq = amino_count / np.sum(amino_count, axis=1)[:, None]

# 2. calculate amino acid frequency in each position for training data

# create a 2 by 2 array to store number of amino acids in each position
amino_count_train = np.zeros((len(train_array[0][1]), len(amino_dict_freq)))

# count amino acids in each position
for i in range(len(train_array)):
	for j in range(len(train_array[i][1])):
		amino_count_train[j][amino_dict_freq[train_array[i][1][j]]] += 1

# calculate amino acid frequency in each position
amino_freq_train = amino_count_train / np.sum(amino_count_train, axis=1)[:, None]

# 3. calculate amino acid frequency in each position for prediction

# create a 2 by 2 array to store number of amino acids in each position
amino_count_pred = np.zeros((len(validation_array[0][1]), len(amino_dict_freq)))

# count amino acids in each position
for i in range(len(peptide_pred)):
	for j in range(len(peptide_pred[i][1])):
		amino_count_pred[j][amino_dict_freq[peptide_pred[i][1][j]]] += 1

# calculate amino acid frequency in each position
amino_freq_pred = amino_count_pred / np.sum(amino_count_pred, axis=1)[:, None]

# 4. scatter plot the validation data frequency vs model generation frequency on a single plot
plt.figure(figsize=(24, 16), facecolor=(1, 1, 1))
plt.rcParams.update({'font.size': 25})
for i in range(len(amino_freq)):
	plt.scatter(amino_freq[i], amino_freq_pred[i], c='r', s=100, marker = "x")
plt.xlabel('Validation data amino acid frequency')
plt.ylabel('Model generation amino acid frequency')
plt.title('1-site frequencies of amino acids in validation data vs model generation')
plt.savefig('model_evaluation/1_site_freq_validation_vs_model_generation.png')
plt.show()

# 5. scatter plot the training data frequency vs model generation frequency on a single plot
plt.figure(figsize=(24, 16), facecolor=(1, 1, 1))
plt.rcParams.update({'font.size': 25})
for i in range(len(amino_freq_train)):
	plt.scatter(amino_freq_train[i], amino_freq_pred[i], c='r', s=100, marker = "x")
plt.xlabel('Training data amino acid frequency')
plt.ylabel('Model generation amino acid frequency')
plt.title('1-site frequencies of amino acids in training data vs model generation')
plt.savefig('model_evaluation/1_site_freq_training_vs_model_generation.png')
plt.show()

# 6. scatter plot the training data frequency vs validation data frequency on a single plot
plt.figure(figsize=(24, 16), facecolor=(1, 1, 1))
plt.rcParams.update({'font.size': 25})
for i in range(len(amino_freq_train)):
	plt.scatter(amino_freq_train[i], amino_freq[i], c='r', s=100, marker = "x")
plt.xlabel('Training data amino acid frequency')
plt.ylabel('Validation data amino acid frequency')
plt.title('1-site frequencies of amino acids in training data vs validation data')
plt.savefig('model_evaluation/1_site_freq_training_vs_validation.png')
plt.show()

# 1. calculate amino acid to amino acid frequency in two positions for validation data

# create a 4D array to store number of amino acids in each position
# 1st dimension: position
# 2nd dimension: position
# 3rd dimension: amino acid
# 4th dimension: amino acid
amino_count_2_sites = np.zeros((len(validation_array[0][1]), len(validation_array[0][1]), len(amino_dict_freq), len(amino_dict_freq)))

# count amino acid to amino acid in each position pair
for i in range(len(validation_array)):
	for j in range(len(validation_array[i][1])):
		for k in range(len(validation_array[i][1])):
			amino_count_2_sites[j][k][amino_dict_freq[validation_array[i][1][j]]][amino_dict_freq[validation_array[i][1][k]]] += 1
            
# calculate amino acid to amino acid frequency in each position pair
amino_freq_2_sites = amino_count_2_sites / np.sum(amino_count_2_sites, axis=(2, 3))[:, :, None, None]

# calculate correlation given by amino_freq_2_sites - product of amino_freq at each position
amino_freq_2_sites_corr = amino_freq_2_sites - np.einsum('ij,kl->ikjl', amino_freq, amino_freq)

# 2. calculate amino acid to amino acid frequency in two positions for training data

# create a 4D array to store number of amino acids in each position
# 1st dimension: position
# 2nd dimension: position
# 3rd dimension: amino acid
# 4th dimension: amino acid
amino_count_2_sites_train = np.zeros((len(train_array[0][1]), len(train_array[0][1]), len(amino_dict_freq), len(amino_dict_freq)))

# count amino acid to amino acid in each position pair
for i in range(len(train_array)):
	for j in range(len(train_array[i][1])):
		for k in range(len(train_array[i][1])):
			amino_count_2_sites_train[j][k][amino_dict_freq[train_array[i][1][j]]][amino_dict_freq[train_array[i][1][k]]] += 1

# calculate amino acid to amino acid frequency in each position pair
amino_freq_2_sites_train = amino_count_2_sites_train / np.sum(amino_count_2_sites_train, axis=(2, 3))[:, :, None, None]

# calculate correlation given by amino_freq_2_sites_train - product of amino_freq at each position
amino_freq_2_sites_train_corr = amino_freq_2_sites_train - np.einsum('ij,kl->ikjl', amino_freq_train, amino_freq_train)

# 3. calculate amino acid to amino acid frequency in two positions for prediction

# create a 4D array to store number of amino acids in each position
# 1st dimension: position
# 2nd dimension: position
# 3rd dimension: amino acid
# 4th dimension: amino acid
amino_count_2_sites_pred = np.zeros((len(validation_array[0][1]), len(validation_array[0][1]), len(amino_dict_freq), len(amino_dict_freq)))

# count amino acid to amino acid in each position pair
for i in range(len(peptide_pred)):
	for j in range(len(peptide_pred[i][1])):
		for k in range(len(peptide_pred[i][1])):
			amino_count_2_sites_pred[j][k][amino_dict_freq[peptide_pred[i][1][j]]][amino_dict_freq[peptide_pred[i][1][k]]] += 1

# calculate amino acid to amino acid frequency in each position pair
amino_freq_2_sites_pred = amino_count_2_sites_pred / np.sum(amino_count_2_sites_pred, axis=(2, 3))[:, :, None, None]

# calculate correlation given by amino_freq_2_sites - product of amino_freq at each position
amino_freq_2_sites_pred_corr = amino_freq_2_sites_pred - np.einsum('ij,kl->ikjl', amino_freq_pred, amino_freq_pred)

# 4. scatter plot the validation data frequency vs model generation frequency on a single plot
plt.figure(figsize=(24, 16), facecolor=(1, 1, 1))
plt.rcParams.update({'font.size': 25})
for i in range(len(amino_freq_2_sites_corr)):
	for j in range(len(amino_freq_2_sites_corr[i])):
		plt.scatter(amino_freq_2_sites_corr[i][j], amino_freq_2_sites_pred_corr[i][j], c='r', s=100, marker = "x")
plt.xlabel('Validation data amino acid frequency')
plt.ylabel('Model generation amino acid to amino acid frequency correlation')
plt.title('2-site connected correlation of amino acids in validation data vs model generation')
plt.savefig('model_evaluation/2_site_corr_validation_vs_model_generation.png')
plt.show()

# 5. scatter plot the training data frequency vs model generation frequency on a single plot
plt.figure(figsize=(24, 16), facecolor=(1, 1, 1))
plt.rcParams.update({'font.size': 25})
for i in range(len(amino_freq_2_sites_train_corr)):
	for j in range(len(amino_freq_2_sites_train_corr[i])):
		plt.scatter(amino_freq_2_sites_train_corr[i][j], amino_freq_2_sites_pred_corr[i][j], c='r', s=100, marker = "x")
plt.xlabel('Training data amino acid frequency')
plt.ylabel('Model generation amino acid to amino acid frequency correlation')
plt.title('2-site connected correlation of amino acids in training data vs model generation')
plt.savefig('model_evaluation/2_site_corr_training_vs_model_generation.png')
plt.show()

# 6. scatter plot the training data frequency vs validation data frequency on a single plot
plt.figure(figsize=(24, 16), facecolor=(1, 1, 1))
plt.rcParams.update({'font.size': 25})
for i in range(len(amino_freq_2_sites_train_corr)):
	for j in range(len(amino_freq_2_sites_train_corr[i])):
		plt.scatter(amino_freq_2_sites_train_corr[i][j], amino_freq_2_sites_corr[i][j], c='r', s=100, marker = "x")
plt.xlabel('Training data amino acid frequency')
plt.ylabel('Validation data amino acid to amino acid frequency correlation')
plt.title('2-site connected correlation of amino acids in training data vs validation data')
plt.savefig('model_evaluation/2_site_corr_training_vs_validation.png')
plt.show()