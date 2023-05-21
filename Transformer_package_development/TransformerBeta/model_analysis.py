import numpy as np
from scipy.stats import percentileofscore
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def calculate_net_charge(peptide_list, reference_list=None):
    amino_positive = ['R', 'K', 'H']
    amino_negative = ['D', 'E']
    charge_list = []

    for peptide in peptide_list:
        charge = 0
        for amino_acid in peptide:
            if amino_acid in amino_positive:
                charge += 1
            elif amino_acid in amino_negative:
                charge -= 1
        charge_list.append(charge)

    if reference_list is not None:
        reference_charge_list = [calculate_net_charge([ref_seq])[0] for ref_seq in reference_list]
        percentile_list = [percentileofscore(reference_charge_list, charge) for charge in charge_list]
        return charge_list, percentile_list
    else:
        return charge_list
    
kyte_doolittle_scale = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

def calculate_hydrophobicity(peptide_list, scale, reference_list=None):
    hydrophobicity_list = []

    for peptide in peptide_list:
        hydrophobicity = np.mean([scale[residue] for residue in peptide if residue in scale])
        hydrophobicity_list.append(hydrophobicity)

    if reference_list is not None:
        reference_hydrophobicity_list = calculate_hydrophobicity(reference_list, scale)
        percentile_list = [percentileofscore(reference_hydrophobicity_list, hydrophobicity) for hydrophobicity in hydrophobicity_list]
        return hydrophobicity_list, percentile_list
    else:
        return hydrophobicity_list
    
molecular_weights = {
    'A': 89.094,   # Alanine
    'R': 174.203,  # Arginine
    'N': 132.119,  # Asparagine
    'D': 133.104,  # Aspartic acid
    'C': 121.154,  # Cysteine
    'E': 147.131,  # Glutamic acid
    'Q': 146.146,  # Glutamine
    'G': 75.067,   # Glycine
    'H': 155.156,  # Histidine
    'I': 131.175,  # Isoleucine
    'L': 131.175,  # Leucine
    'K': 146.189,  # Lysine
    'M': 149.208,  # Methionine
    'F': 165.192,  # Phenylalanine
    'P': 115.132,  # Proline
    'S': 105.093,  # Serine
    'T': 119.120,  # Threonine
    'W': 204.228,  # Tryptophan
    'Y': 181.191,  # Tyrosine
    'V': 117.148,  # Valine
}

def calculate_molecular_weights(peptide_list, scale, reference_list=None):
    molecular_weight_list = []

    for peptide in peptide_list:
        molecular_weight = np.mean([scale[residue] for residue in peptide if residue in scale])
        molecular_weight_list.append(molecular_weight)

    if reference_list is not None:
        reference_molecular_weight_list = calculate_molecular_weights(reference_list, scale)
        percentile_list = [percentileofscore(reference_molecular_weight_list, molecular_weight) for molecular_weight in molecular_weight_list]
        return molecular_weight_list, percentile_list
    else:
        return molecular_weight_list


def calculate_isoelectric_points(peptide_list, reference_list=None):
    isoelectric_point_list = []

    for peptide in peptide_list:
        analysis = ProteinAnalysis(peptide)
        isoelectric_point = analysis.isoelectric_point()
        isoelectric_point_list.append(isoelectric_point)

    if reference_list is not None:
        reference_isoelectric_point_list = calculate_isoelectric_points(reference_list)
        percentile_list = [percentileofscore(reference_isoelectric_point_list, isoelectric_point) for isoelectric_point in isoelectric_point_list]
        return isoelectric_point_list, percentile_list
    else:
        return isoelectric_point_list
    
def calculate_aromaticity(peptide_list, reference_list=None):
    aromaticity_list = []

    for peptide in peptide_list:
        analysis = ProteinAnalysis(peptide)
        aromaticity = analysis.aromaticity()
        aromaticity_list.append(aromaticity)

    if reference_list is not None:
        reference_aromaticity_list = calculate_aromaticity(reference_list)
        percentile_list = [percentileofscore(reference_aromaticity_list, aromaticity) for aromaticity in aromaticity_list]
        return aromaticity_list, percentile_list
    else:
        return aromaticity_list
    
def calculate_novelty_scores(peptide_list, reference_list):
    novelty_scores = []

    for peptide_search in peptide_list:
        min_dist = 100
        for train_peptide in reference_list:
            for j in range(len(peptide_search)):
                dist = 0
                for k in range(len(peptide_search)):
                    if peptide_search[k] != train_peptide[k]:
                        dist += 1
                if dist < min_dist:
                    min_dist = dist
        novelty_scores.append(min_dist)

    return np.array(novelty_scores)

def calculate_novelty_scores_and_reference_targets(peptide_list, reference_list, target, target_reference_list):
    novelty_scores = []
    target_novelty_scores = []

    for peptide_search in peptide_list:
        min_dist = 100
        min_dist_indices = []
        for idx, train_peptide in enumerate(reference_list):
            dist = 0
            for k in range(len(peptide_search)):
                if peptide_search[k] != train_peptide[k]:
                    dist += 1
            if dist < min_dist:
                min_dist = dist
                min_dist_indices = [idx]
            elif dist == min_dist:
                min_dist_indices.append(idx)

        novelty_scores.append(min_dist)
        
        min_target_dist = 100
        for idx in min_dist_indices:
            reference_target = target_reference_list[idx]
            target_dist = sum(a != b for a, b in zip(target, reference_target))
            if target_dist < min_target_dist:
                min_target_dist = target_dist

        target_novelty_scores.append(min_target_dist)

    return np.array(novelty_scores), np.array(target_novelty_scores)

def generate_output_table(peptide_candidates, peptide_candidates_prob, reference_list, output_file='output_table.xlsx', cluster_labels='', target = None, target_reference_list = None, target_novelty_score = ''):
    print('Clustering analysis Done')
    ranks = np.arange(1, len(peptide_candidates) + 1) 
    print('rank analysis Done')
    charge, charge_percentile = calculate_net_charge(peptide_candidates, reference_list)
    print('Net Charge Analysis Done')
    hydrophobicity, hydrophobicity_percentile = calculate_hydrophobicity(peptide_candidates, kyte_doolittle_scale, reference_list)
    print('Hydrophobicity Analysis Done')
    molecular_weight, molecular_weight_percentile = calculate_molecular_weights(peptide_candidates, molecular_weights, reference_list)
    print('Molecular Weight Analysis Done')
    isoelectric_point, isoelectric_point_percentile = calculate_isoelectric_points(peptide_candidates, reference_list)
    print('Isoelectric Point Analysis Done')
    aromaticity, aromaticity_percentile = calculate_aromaticity(peptide_candidates, reference_list)
    print('Aromaticity Analysis Done')
    if target == None and target_reference_list == None and target_novelty_score == '':
        novelty_score = calculate_novelty_scores(peptide_candidates, reference_list)
        print('Novelty Score Analysis Done')
    else:
        novelty_score, target_novelty_score = calculate_novelty_scores_and_reference_targets(peptide_candidates, reference_list, target, target_reference_list)
        print('Novelty Score and Target Novelty Score Analysis Done')
    # reverse each of the sequence in peptide_candidates
    peptide_candidates_synthesis = [peptide[::-1] for peptide in peptide_candidates]

    data = {
        'Rank': ranks,
        'Cluster Label': cluster_labels,
        'Target Sequence': '',
        'Designed Complementary Peptide': peptide_candidates,
        'Synthesis Complementary Peptide': peptide_candidates_synthesis,
        'Probability': peptide_candidates_prob,
        'Novelty Score': novelty_score,
        'Target Novelty Score': target_novelty_score,
        'CamSol Solubility Score': '',
        'Net Charge': charge,
        'Net Charge Percentile': charge_percentile,
        'Hydrophobicity': hydrophobicity,
        'Hydrophobicity Percentile': hydrophobicity_percentile,
        'Molecular Weight': molecular_weight,
        'Molecular Weight Percentile': molecular_weight_percentile,
        'Isoelectric Point': isoelectric_point,
        'Isoelectric Point Percentile': isoelectric_point_percentile,
        'Aromaticity': aromaticity,
        'Aromaticity Percentile': aromaticity_percentile
    }

    df = pd.DataFrame(data)
    df.to_excel(output_file, index=False)