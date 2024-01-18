import math
import numpy as np
import json

def beta_strand(cif_file, pLDDT_file=None, PAE_file=None, pLLDDT_min=0, PAE_max=100, least_length=2):
	"""
	Version 4
	Convert one sturcture to beta strands that face each other
	"""

	# open pLDDT file
	if pLDDT_file is not None:
		with open(pLDDT_file, 'r') as f:
			pLDDT = json.load(f)
		pLDDT_residueNumber = pLDDT['residueNumber']
		pLDDT_confidenceScore = pLDDT['confidenceScore']
		pLDDT_confidenceCategory = pLDDT['confidenceCategory']

	# open PAE file
	if PAE_file is not None:
		with open(PAE_file, 'r') as f:
			PAE = json.load(f)
		PAE_predicted_aligned_error = PAE[0]['predicted_aligned_error']
		PAE_max_predicted_aligned_error= PAE[0]['max_predicted_aligned_error']

	# extract all amino residue
	# extract beta strand (amino num + pieces)
	# extract coordinates of carbon alpha for each amino residue in beta_strand_amino_num
	amino_residue = []
	beta_strand_amino_num = []
	beta_strand_pieces = []
	beta_strand_amino_coordinates = []
	with open(cif_file, 'r') as f:
		for line in f:
			if line.startswith('ATOM') and line.split()[3] == 'CA':
				amino_residue.append(line.split()[-1])
			if 'STRN' in line and line.startswith('A'):
				start_index = int(line.split()[2])
				end_index = int(line.split()[9])
				# appends values in bewteen start and end index
				beta_piece = list(range(start_index, end_index + 1))
				if len(beta_piece) >= least_length:
					beta_strand_amino_num += beta_piece
					beta_strand_pieces.append(beta_piece) 
			if line.startswith('ATOM') and line.split()[3] == 'CA' and int(line.split()[-2]) in beta_strand_amino_num:
				beta_strand_amino_coordinates.append([line.split()[10], line.split()[11], line.split()[12]])

	# return empty list if less than 2 beta strands
	if len(beta_strand_pieces) < 2:
		beta_strand = []
		beta_strand = np.array(beta_strand, dtype=object)
		return beta_strand

	# calculate distance between each pair of beta strand amino num
	beta_strand_amino_distance = np.zeros((len(beta_strand_amino_coordinates), len(beta_strand_amino_coordinates)))
	for i in range(len(beta_strand_amino_coordinates)):
		for j in range(len(beta_strand_amino_coordinates)):
			beta_strand_amino_distance[i][j] = math.sqrt((float(beta_strand_amino_coordinates[i][0]) - float(beta_strand_amino_coordinates[j][0]))**2 
														+ (float(beta_strand_amino_coordinates[i][1]) - float(beta_strand_amino_coordinates[j][1]))**2 
														+ (float(beta_strand_amino_coordinates[i][2]) - float(beta_strand_amino_coordinates[j][2]))**2)
	
	# save pLDDT and PAE only for beta strand
	if pLDDT_file is not None and PAE_file is not None:
		beta_strand_pLDDT = []
		beta_strand_PAE = np.zeros((len(beta_strand_amino_num), len(beta_strand_amino_num)))
		for i in range(len(beta_strand_amino_num)):
			beta_strand_pLDDT.append(pLDDT_confidenceScore[beta_strand_amino_num[i] - 1])
			for j in range(len(beta_strand_amino_num)):
				beta_strand_PAE[i][j] = PAE_predicted_aligned_error[beta_strand_amino_num[i] - 1][beta_strand_amino_num[j] - 1]
	elif pLDDT_file is not None and PAE_file is None:
		beta_strand_pLDDT = []
		for i in range(len(beta_strand_amino_num)):
			beta_strand_pLDDT.append(pLDDT_confidenceScore[beta_strand_amino_num[i] - 1])
	elif pLDDT_file is None and PAE_file is not None:
		beta_strand_PAE = np.zeros((len(beta_strand_amino_num), len(beta_strand_amino_num)))
		for i in range(len(beta_strand_amino_num)):
			for j in range(len(beta_strand_amino_num)):
				beta_strand_PAE[i][j] = PAE_predicted_aligned_error[beta_strand_amino_num[i] - 1][beta_strand_amino_num[j] - 1]

	# loop thourgh pieces of beta strands
	beta_strand = []
	type = 0 # type 0 is parallel and type 1 is anti-parallel
	for piece in beta_strand_pieces:
		beta_strand_target = piece.copy()

		# 1. locate amino residues that are closest to beta_strand_target but not in the same piece
		# 1. locate amino residues that are 2nd closest to beta_strand_target but not in the same piece
		beta_strand_target_closest_residues = []
		beta_strand_target_2nd_closest_residues = []
		for amino_num in beta_strand_target:
			# exlude amino residues in the same piece
			index = beta_strand_amino_num.index(amino_num)
			beta_strand_amino_distance_exclude_current_piece = beta_strand_amino_distance.copy()
			for i in beta_strand_target:
				index_i = beta_strand_amino_num.index(i)
				beta_strand_amino_distance_exclude_current_piece[index][index_i] = 1000
			# find the closest amino num in beta_strand_amino_num and it is not in the same piece
			closest_amino_num = beta_strand_amino_num[np.argmin(beta_strand_amino_distance_exclude_current_piece[index])]
			beta_strand_target_closest_residues.append(closest_amino_num)
			# find the 2nd closest amino num in beta_strand_amino_num and it is not in the same piece
			second_closest_amino_num = beta_strand_amino_num[np.argsort(beta_strand_amino_distance_exclude_current_piece[index])[1]]
			beta_strand_target_2nd_closest_residues.append(second_closest_amino_num)
		
		# 2. find the closest beta strand by a majority vote of the closest residues
		beta_strand_closest = []
		beta_strand_closest_num = 0
		for piece_comp in beta_strand_pieces:
			closest_num = 0
			for amino_num in beta_strand_target_closest_residues:
				if amino_num in piece_comp:
					closest_num += 1
			if closest_num > beta_strand_closest_num:
				beta_strand_closest = piece_comp.copy()
				beta_strand_closest_num = closest_num

		# 3. build a potential contact list
		beta_strand_target_potential_contact = []
		for i, j in zip(beta_strand_target_closest_residues, beta_strand_target_2nd_closest_residues):
			if i in beta_strand_closest:
				beta_strand_target_potential_contact.append(i)
				continue
			elif j in beta_strand_closest:
				beta_strand_target_potential_contact.append(j)
				continue
		contact_num_min = min(beta_strand_target_potential_contact)
		contact_num_max = max(beta_strand_target_potential_contact)
		contact_length = contact_num_max - contact_num_min + 1
			
		if beta_strand_target_potential_contact[-1] > beta_strand_target_potential_contact[0]:
			type = 0
		else:
			type = 1
			# reverse the beta_strand_closest
			beta_strand_closest = beta_strand_closest[::-1]

		# 4. build beta_strand_comp
		if contact_length > len(beta_strand_target) or contact_length > len(beta_strand_closest):
			contact_length = min(len(beta_strand_target), len(beta_strand_closest))
		if contact_length != len(beta_strand_target) or contact_length != len(beta_strand_closest):
			# break beta_strand_target into pieces of contact_length
			# break beta_strand_closest into pieces of contact_length
			beta_strand_target_list = [beta_strand_target[i:i+contact_length] for i in range(0, len(beta_strand_target)-contact_length+1)]
			beta_strand_closest_list = [beta_strand_closest[i:i+contact_length] for i in range(0, len(beta_strand_closest)-contact_length+1)]
			# find the beta_strand_comp that has shortest distance between beta_strand_target
			beta_strand_target = []
			beta_strand_comp = []
			shortest_distance = 1000
			distance = 0
			for i in range(len(beta_strand_target_list)):
				for j in range(len(beta_strand_closest_list)):
					for k in range(contact_length):
						index = beta_strand_amino_num.index(beta_strand_target_list[i][k])
						index_comp = beta_strand_amino_num.index(beta_strand_closest_list[j][k])
						distance += beta_strand_amino_distance[index][index_comp]
					if distance < shortest_distance:
						shortest_distance = distance
						beta_strand_target = beta_strand_target_list[i].copy()
						beta_strand_comp = beta_strand_closest_list[j].copy()
					distance = 0
		elif contact_length == len(beta_strand_target) and contact_length == len(beta_strand_closest):
			beta_strand_target = beta_strand_target.copy()
			beta_strand_comp = beta_strand_closest.copy()

		# 5. check all residues satisfy the following conditions:
		# pLDDt >= pLLDDT_min
		# PAE <= PAE_max
		# length of beta_strand > least_length
		if pLDDT_file is not None and PAE_file is not None:
			for i in range(len(beta_strand_target)):
				# find the index of beta_strand_target[i] in beta_strand_amino_num
				index = beta_strand_amino_num.index(beta_strand_target[i])
				# find the index of beta_strand_comp[i] in beta_strand_amino_num
				index_comp = beta_strand_amino_num.index(beta_strand_comp[i])
				# if pLDDt < pLLDDT_min or PAE > PAE_max, remove the residue
				if beta_strand_pLDDT[index] < pLLDDT_min or beta_strand_pLDDT[index_comp] < pLLDDT_min or beta_strand_PAE[index, index_comp] > PAE_max:
					# remove the whole piece
					beta_strand_target = []
					beta_strand_comp = []
					break
		elif pLDDT_file is not None and PAE_file is None:
			for i in range(len(beta_strand_target)):
				# find the index of beta_strand_target[i] in beta_strand_amino_num
				index = beta_strand_amino_num.index(beta_strand_target[i])
				# find the index of beta_strand_comp[i] in beta_strand_amino_num
				index_comp = beta_strand_amino_num.index(beta_strand_comp[i])
				# if pLDDt < pLLDDT_min, remove the residue
				if beta_strand_pLDDT[index] < pLLDDT_min or beta_strand_pLDDT[index_comp] < pLLDDT_min:
					# remove the whole piece
					beta_strand_target = []
					beta_strand_comp = []
					break
		elif pLDDT_file is None and PAE_file is not None:
			for i in range(len(beta_strand_target)):
				# find the index of beta_strand_target[i] in beta_strand_amino_num
				index = beta_strand_amino_num.index(beta_strand_target[i])
				# find the index of beta_strand_comp[i] in beta_strand_amino_num
				index_comp = beta_strand_amino_num.index(beta_strand_comp[i])
				# if PAE > PAE_max, remove the residue
				if beta_strand_PAE[index, index_comp] > PAE_max:
					# remove the whole piece
					beta_strand_target = []
					beta_strand_comp = []
					break
		if len(beta_strand_target) < least_length:
			continue

		# 6. replace amino num with actual amino residue if the length is longer than 2
		beta_strand_target_residue = ''
		beta_strand_comp_residue = ''
		for amino_num in beta_strand_target:
			beta_strand_target_residue += amino_residue[amino_num - 1]
		for amino_num in beta_strand_comp:
			beta_strand_comp_residue += amino_residue[amino_num - 1]
		beta_strand.append([[beta_strand_target[0], beta_strand_target[-1]], [beta_strand_comp[0], beta_strand_comp[-1]], beta_strand_target_residue, beta_strand_comp_residue, type])

	# convert the list to array
	beta_strand = np.array(beta_strand, dtype=object)

	return beta_strand

if __name__ == '__main__':
	# please speficify:
	cif_file = 'AF-A0A004-F1-model_v3.cif' # must specify !!!
	pLDDT_file = 'AF-A0A004-F1-confidence_v3.json' # default is None
	PAE_file = 'AF-A0A004-F1-predicted_aligned_error_v3.json' # default is None
	pLLDDT_min = 70 # default is 0 if pLDDT_file is not None
	PAE_max = 10 # default is 100 if PAE_file is not None
	least_length = 3 # default is 2

	# run the function
	beta_strand_output = beta_strand(cif_file, pLDDT_file, PAE_file, pLLDDT_min, PAE_max, least_length)

	# np save the output
	file_name = cif_file.split('.')[0]
	np.save(file_name + '_beta_strand.npy', beta_strand_output)
