import numpy as np
import pandas as pd
import random


def simulate_data_congas(ncells = 1500, n_segments = 10, subclones = 3,mixture_proportion_rna = np.array([0.3,0.6,0.1]),
                         mixture_proportion_atac = np.array([0.4,0.4,0.2]),
                         simulate_rna = True, simulate_atac = True, multiome = False, seed = 3):
    random.seed(seed)
    np.random.seed(seed)
    
    
    Majors = [np.random.randint(1, 5, n_segments) for i in range(subclones)]
    minors = [np.random.randint(0, Majors[i], n_segments) for i in range(subclones)] 
    totals = [Majors[i] + minors[i] for i in range(subclones)] 


    # sample normalization factors 
    shape, scale = 2., 2.
    # MAF parameters 
    Ntrials = 50

    if simulate_rna or multiome:
        norm_factors_rna = np.random.gamma(shape, scale, ncells)
        shape_segment_rna, scale_segment_rna = 20.,2.
    if simulate_atac or multiome:
        norm_factors_atac = np.random.gamma(shape, scale, ncells)
        shape_segment_atac, scale_segment_atac = 10.,2.
    # Generate cells 
    if simulate_rna or multiome:
        subclones_matrices_exp_rna = []
        subclones_matrices_MAF_rna = []
        for i in range(subclones):
            tmp_list_expr = []
            tmp_list_MAF = []
            for j in range(n_segments):
                tmp_list_expr.append(np.random.gamma(shape = shape_segment_rna,scale =  scale_segment_rna,size = ncells) * norm_factors_rna * totals[i][j])
                tmp_list_MAF.append(generate_MAF_given_CN(Majors[i][j], minors[i][j], Ntrials, ncells))
            subclones_matrices_exp_rna.append(np.concatenate(tmp_list_expr).reshape([ncells, n_segments]).round())
            subclones_matrices_MAF_rna.append(np.concatenate(tmp_list_MAF).reshape([ncells, n_segments]))
        expression_rna = generate_proportion_matrix(subclones_matrices_exp_rna, mixture_proportion_rna)
        MAF_rna = generate_proportion_matrix(subclones_matrices_MAF_rna, mixture_proportion_rna)
        labels_rna = generate_labels(ncells, mixture_proportion_rna)
    else:
        expression_rna = None
        MAF_rna = None
        labels_rna = none

    if simulate_atac or multiome:
        subclones_matrices_exp_atac = []
        subclones_matrices_MAF_atac = []
        for i in range(subclones):
            tmp_list_expr = []
            tmp_list_MAF = []
            for j in range(n_segments):
                tmp_list_expr.append(np.random.gamma(shape = shape_segment_atac,scale =  scale_segment_atac,size = ncells) * norm_factors_atac * totals[i][j])
                tmp_list_MAF.append(generate_MAF_given_CN(Majors[i][j], minors[i][j], Ntrials, ncells))
            subclones_matrices_exp_atac.append(np.concatenate(tmp_list_expr).reshape([ncells, n_segments]).round())
            subclones_matrices_MAF_atac.append(np.concatenate(tmp_list_MAF).reshape([ncells, n_segments])) 
        if multiome:
            mixture_proportion_atac = mixture_proportion_rna
        expression_atac = generate_proportion_matrix(subclones_matrices_exp_atac, mixture_proportion_atac)
        MAF_atac = generate_proportion_matrix(subclones_matrices_MAF_atac, mixture_proportion_atac)
        labels_atac = generate_labels(ncells, mixture_proportion_atac)
    else:
        expression_atac = None
        MAF_atac = None
        labels_atac = None

    mixture_proportion_bulk = (mixture_proportion_atac + mixture_proportion_rna) / 2
    Major_bulk = (np.array(Majors) * np.expand_dims(mixture_proportion_bulk, axis = 1)).sum(axis = 0).round()
    minor_bulk = (np.array(minors) * np.expand_dims(mixture_proportion_bulk, axis = 1)).sum(axis = 0).round()
    total_bulk = Major_bulk + minor_bulk

    return_df = {"atac" : {"exp": expression_atac, "MAF": MAF_atac, "labels" : labels_atac, 
                           "norm_factors" : norm_factors_atac, "theta_shape_atac" : np.ones(n_segments) * shape_segment_atac,
                          "theta_rate_atac" : np.ones(n_segments) * 1/scale_segment_atac },
                "rna" : {"exp": expression_rna, "MAF": MAF_rna, "labels" : labels_rna, "norm_factors" : norm_factors_rna, "theta_shape_rna" : np.ones(n_segments) * shape_segment_rna,
                          "theta_rate_rna" : np.ones(n_segments) * 1/scale_segment_rna },
                 "CNA_subclones" : {"Major" : pd.DataFrame(Majors), "minor" : pd.DataFrame(minors), "total" : pd.DataFrame(totals)},
                 "CNA_bulk" : pd.DataFrame({"Major" : Major_bulk, "minor" : minor_bulk, "total" : total_bulk} ), 
                "hyperparams" : {"ncells" : ncells, "is_multiome" : multiome}}
    return return_df



def generate_MAF_given_CN(Major,minor, ntrials = 50, ncells = 1000):
    mean = (minor / (Major + minor)) + 1e-16
    return np.random.beta(mean * ntrials, (1 - mean) * ntrials, ncells)      

def generate_proportion_matrix(matrices, proportions):
    samples = []
    for matrix, proportion in zip(matrices, proportions):
            # Determine the number of samples to draw from this matrix
        # Determine the number of rows to slice from this matrix
        num_rows_from_matrix = int(matrix.shape[0] * proportion)
        
        # Slice the first num_rows_from_matrix rows from the matrix
        sliced_matrix = matrix[:num_rows_from_matrix, :]
        
        # Add the sliced matrix to our samples list
        samples.append(sliced_matrix)

        # Concatenate all samples
        output_matrix = np.concatenate(samples, axis=0)

    return output_matrix

def generate_labels(num_samples, proportions):    
    # Create an empty list to store the labels
    labels = []

    # For each proportion
    for i, proportion in enumerate(proportions):
        # Determine the number of samples for this label
        num_samples_for_label = int(num_samples * proportion)
        
        # Assign a label to these samples
        labels += [i] * num_samples_for_label

    return np.array(labels)
    