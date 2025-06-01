import pickle
import gzip

def decode_from_latent_space(latent_points, grammar_model):

    decode_attempts = 500
    decoded_RNAs = []
    for i in range(decode_attempts):
        current_decoded_RNAs = grammar_model.decode(latent_points)
        current_decoded_RNAs = [ x if x != '' else 'Sequence too long' for x in current_decoded_RNAs ]
        decoded_RNAs.append(current_decoded_RNAs)

    # We see which ones are decoded by rdkit
    
    RNA_series = []
    for i in range(decode_attempts):
        RNA_series.append([])
        for j in range(latent_points.shape[ 0 ]):
            rna_seq = np.array([ decoded_RNAs[ i ][ j ] ]).astype('str')[ 0 ]
            #print("rna_seq")
            #print(rna_seq)
            if rna_seq == 'Sequence too long': #is None:
                RNA_series[ i ].append(None)
                #print("RNA_series_None")
                #print(RNA_series[i])
            else:
                RNA_series[ i ].append(rna_seq)
                #print("RNA_series")
                #print(RNA_series[i])
    
    import collections

    print("RNA_series")
    print(RNA_series[i])

    decoded_RNAs = np.array(decoded_RNAs)
    RNA_series = np.array(RNA_series)
    #rdkit_molecules = np.array(rdkit_molecules)

    final_rnas = []
    g_c_scores = []
    MFE_scores = []
    Whole_rnas = []
    Whole_MFE_scores = []
    Whole_gcs = []
    Whole_gcs_scores = []
    
    #final_rnas.append(rna_seq)

    #target_structure = "(((((((....(((...........)))((((((((..(((((((((((((((((((...(((((......))))).)))))).)))))))))))))..))))))))..)))))))" #"(((((......)))))"
    #len_t = len(target_structure)

    for i in range(latent_points.shape[ 0 ]):
        #rna_seq = "A"
        #structure, energy = RNA.fold(rna_seq)
        #min_energy = energy
        min_H_dist = 100000   # Minimum Hamming distance
        min_energy = 0
        
        for j in range(decode_attempts):
            if ( ~np.equal(RNA_series[ j, i ], None)):
                structure, energy = RNA.fold(RNA_series[ j, i ])
                Whole_rnas.append(RNA_series[ j, i ])
                Whole_MFE_scores.append(energy)
                tokens = RNA_vae.tokenization(RNA_series[ j, i ])
                Whole_gcs_scores.append(abs(RNA_vae.g_c(tokens)-0.5))
                Whole_gcs.append(RNA_vae.g_c(tokens))
                '''
                # Computing the hamming distance between the target sequence and the structure:
                len_s = len(structure)
                min_len = min(len_t,len_s) # Min length
                max_len = max(len_t,len_s) # Max length

                H_dist = 0 # Hamming distance
                for t in range(min_len):
                    if (target_structure[t] != structure[t]):
                        H_dist = H_dist + 1
    
                H_dist = H_dist + (max_len - min_len)
                # print("H_dist", H_dist)
                #target_structure_scores.append(H_dist)
                #print(i)
                
                
                if (H_dist<min_H_dist):
                    min_H_dist = H_dist
                    rna_seq = RNA_series[ j, i ]
                    min_energy = energy
                '''

                if (energy<min_energy):
                    rna_seq = RNA_series[ j, i ]
                    min_energy = energy
                


        '''    
        structure, energy = RNA.fold(valid_RNA_final[ i ])
        aux = RNA_series[ ~np.equal(RNA_series[ :, i ], None) , i ]
        print("aux")
        print(aux)
        if len(aux) > 0:
            #rna_seq = list(aux.items())[ np.argmax(list(aux.values())) ][ 0 ]
            rna_seq = list(aux.items())[0][0]
            print("rna_seq")
            print(rna_seq)

        else:
            rna_seq = None
        '''
        
        tokens = RNA_vae.tokenization(rna_seq)
        g_c_scores.append(abs(RNA_vae.g_c(tokens)-0.5))
        #print(i)
        
        print("rna_seq")
        print(rna_seq)
        print("min_H_dist")
        print(min_H_dist)
        final_rnas.append(rna_seq)
        MFE_scores.append(min_energy)

    return final_rnas, MFE_scores, g_c_scores, Whole_rnas, Whole_MFE_scores, Whole_gcs, Whole_gcs_scores

# We define the functions used to load and save objects

def save_object(obj, filename):

    """
    Function that saves an object to a file using pickle
    """

    result = pickle.dumps(obj)
    with gzip.GzipFile(filename, 'wb') as dest: dest.write(result)
    dest.close()


def load_object(filename):

    """
    Function that loads an object from a file using pickle
    """

    with gzip.GzipFile(filename, 'rb') as source: result = source.read()
    ret = pickle.loads(result)
    source.close()

    return ret

from sparse_gp import SparseGP

import scipy.stats    as sps

import numpy as np

random_seed = int(np.loadtxt('../random_seed.txt'))
np.random.seed(random_seed)

alpha = 1
print("alpha",alpha)

# We load the data

X = np.loadtxt('../../latent_features_and_targets_grammar/latent_faetures_t.txt')

MFE_scores = np.loadtxt('../../latent_features_and_targets_grammar/MFE_scores_t.txt')
g_c_scores = np.loadtxt('../../latent_features_and_targets_grammar/g_c_scores_t.txt')

MFE_scores_normalized = (np.array(MFE_scores) - np.mean(MFE_scores)) / np.std(MFE_scores)
g_c_scores_normalized = (np.array(g_c_scores) - np.mean(g_c_scores)) / np.std(g_c_scores)

targets = MFE_scores_normalized + alpha*g_c_scores_normalized




#y = -np.loadtxt('../../latent_features_and_targets_grammar/targets_t.txt')
y = -targets
y = y.reshape((-1, 1))

n = X.shape[ 0 ]
permutation = np.random.choice(n, n, replace = False)

X_train = X[ permutation, : ][ 0 : np.int(np.round(0.9 * n)), : ]
X_test = X[ permutation, : ][ np.int(np.round(0.9 * n)) :, : ]

y_train = y[ permutation ][ 0 : np.int(np.round(0.9 * n)) ]
y_test = y[ permutation ][ np.int(np.round(0.9 * n)) : ]

import os.path

np.random.seed(random_seed)

iteration = 0

while iteration < 10:

    # We fit the GP

    np.random.seed(iteration * random_seed)
    M = 500
    sgp = SparseGP(X_train, 0 * X_train, y_train, M)
    sgp.train_via_ADAM(X_train, 0 * X_train, y_train, X_test, X_test * 0,  \
        y_test, minibatch_size = 10 * M, max_iterations = 50, learning_rate = 0.0005)

    pred, uncert = sgp.predict(X_test, 0 * X_test)
    error = np.sqrt(np.mean((pred - y_test)**2))
    testll = np.mean(sps.norm.logpdf(pred - y_test, scale = np.sqrt(uncert)))
    print('Test RMSE: ') 
    print(error)
    print('Test ll: ') 
    print(testll)

    pred, uncert = sgp.predict(X_train, 0 * X_train)
    error = np.sqrt(np.mean((pred - y_train)**2))
    trainll = np.mean(sps.norm.logpdf(pred - y_train, scale = np.sqrt(uncert)))
    print('Train RMSE: ') 
    print(error)
    print('Train ll: ') 
    print(trainll)

    # We load the decoder to obtain the molecules

    #from rdkit.Chem import MolFromSmiles, MolToSmiles
    #from rdkit.Chem import Draw
    #import image
    #import copy
    #import time

    import sys
    sys.path.insert(0, '../../../')
    import RNA_vae
    import RNA
    grammar_weights = '../../../pretrained/RNA_vae_grammar_SLF_h_100k100_c234_L10_E100_batchB.hdf5' #RNA_vae_grammar_h_conditional_switch100_c234_L10_E100_batchB
    grammar_model = RNA_vae.RNAGrammarModel(grammar_weights)

    # We pick the next 50 inputs

    next_inputs = sgp.batched_greedy_ei(50, np.min(X_train, 0), np.max(X_train, 0))

    valid_RNA_final, MFE_final, g_c_final, Whole_rnas, Whole_MFE_scores, Whole_gcs, Whole_gcs_scores = decode_from_latent_space(next_inputs, grammar_model)

    #from rdkit.Chem import Descriptors
    #from rdkit.Chem import MolFromSmiles, MolToSmiles

    new_features = next_inputs

    save_object(valid_RNA_final, "results/valid_smiles{}.dat".format(iteration))
    '''
    logP_values = np.loadtxt('../../latent_features_and_targets_grammar/logP_values.txt')
    SA_scores = np.loadtxt('../../latent_features_and_targets_grammar/SA_scores.txt')
    cycle_scores = np.loadtxt('../../latent_features_and_targets_grammar/cycle_scores.txt')
    SA_scores_normalized = (np.array(SA_scores) - np.mean(SA_scores)) / np.std(SA_scores)
    logP_values_normalized = (np.array(logP_values) - np.mean(logP_values)) / np.std(logP_values)
    cycle_scores_normalized = (np.array(cycle_scores) - np.mean(cycle_scores)) / np.std(cycle_scores)
    '''
    #target_structure_scores = np.loadtxt('../../latent_features_and_targets_grammar/target_structure_scores.txt')
    #MFE_scores = np.loadtxt('../../latent_features_and_targets_grammar/MFE_scores_t.txt')
    #g_c_scores = np.loadtxt('../../latent_features_and_targets_grammar/g_c_scores_t.txt')
    #structure_scores = np.loadtxt('../../latent_features_and_targets_grammar/structure_scores.txt')
    #target_structure_scores_normalized = (np.array(target_structure_scores) - np.mean(target_structure_scores)) / np.std(target_structure_scores)
    #MFE_scores_normalized = (np.array(MFE_scores) - np.mean(MFE_scores)) / np.std(MFE_scores)
    #g_c_scores_normalized = (np.array(g_c_scores) - np.mean(g_c_scores)) / np.std(g_c_scores)
             
    
    #targets = MFE_scores_normalized + alpha*g_c_scores_normalized # target_structure_scores_normalized # + logP_values_normalized + cycle_scores_normalized

    #import sascorer
    #import networkx as nx
    #from rdkit.Chem import rdmolops
    #target_structure = "(((((((....(((...........)))((((((((..(((((((((((((((((((...(((((......))))).)))))).)))))))))))))..))))))))..)))))))" #"(((((......)))))"
    #len_t = len(target_structure)

    scores = []
    #current_scores = []
    for i in range(len(valid_RNA_final)):
        if valid_RNA_final[ i ] is not None:
            structure, energy = RNA.fold(valid_RNA_final[ i ])
            '''
            # Computing the hamming distance between the target sequence and the structure:
            len_s = len(structure)
            min_len = min(len_t,len_s) # Min length
            max_len = max(len_t,len_s) # Max length

            H_dist = 0 # Hamming distance
            for t in range(min_len):
                if (target_structure[t] != structure[t]):
                    H_dist = H_dist + 1
    
            H_dist = H_dist + (max_len - min_len)
            current_target_structure_scores = H_dist
            '''
            current_MFE_scores = MFE_final[i]
            current_g_c_scores = g_c_final[i] 
            #current_MFE_score = energy
            #current_log_P_value = Descriptors.MolLogP(MolFromSmiles(valid_smiles_final[ i ]))
            #current_SA_score = -sascorer.calculateScore(MolFromSmiles(valid_smiles_final[ i ]))
            #cycle_list = nx.cycle_basis(nx.Graph(rdmolops.GetAdjacencyMatrix(MolFromSmiles(valid_smiles_final[ i ]))))
            '''
            if len(cycle_list) == 0:
                cycle_length = 0
            else:
                cycle_length = max([ len(j) for j in cycle_list ])
            if cycle_length <= 6:
                cycle_length = 0
            else:
                cycle_length = cycle_length - 6

            current_cycle_score = -cycle_length
            '''
         
            #current_target_structure_scores_normalized = (np.array(current_target_structure_scores) - np.mean(target_structure_scores)) / np.std(target_structure_scores)
            current_MFE_scores_normalized = (np.array(current_MFE_scores) - np.mean(MFE_scores)) / np.std(MFE_scores)
            current_g_c_scores_normalized = (np.array(current_g_c_scores) - np.mean(g_c_scores)) / np.std(g_c_scores)
            #current_SA_score_normalized = (current_SA_score - np.mean(SA_scores)) / np.std(SA_scores)
            #current_log_P_value_normalized = (current_log_P_value - np.mean(logP_values)) / np.std(logP_values)
            #current_cycle_score_normalized = (current_cycle_score - np.mean(cycle_scores)) / np.std(cycle_scores)

            #score = (current_SA_score_normalized + current_log_P_value_normalized + current_cycle_score_normalized)
            score = current_MFE_scores_normalized+alpha*current_g_c_scores_normalized # + current_target_structure_scores_normalized
        else:
            score = -max(y)[ 0 ]

        scores.append(-score)
        print(i)

    print(valid_RNA_final)
    print(scores)

    save_object(scores, "results/scores{}.dat".format(iteration))

    with open("results/Whole_rnas{}.txt".format(iteration), 'w') as txt_file:  # Specify encoding here
        for i in range(len(Whole_rnas)):
            txt_file.write(Whole_rnas[i])
            txt_file.write('\n')

    with open("results/Whole_MFE_scores{}.txt".format(iteration), 'w') as txt_file:  # Specify encoding here
        for i in range(len(Whole_MFE_scores)):
            txt_file.write(str(Whole_MFE_scores[i]))
            txt_file.write('\n')
            

    with open("results/Whole_gcs{}.txt".format(iteration), 'w') as txt_file:  # Specify encoding here
        for i in range(len(Whole_gcs)):
            txt_file.write(str(Whole_gcs[i]))
            txt_file.write('\n')

    with open("results/Whole_gcs_scores{}.txt".format(iteration), 'w') as txt_file:  # Specify encoding here
        for i in range(len(Whole_gcs_scores)):
            txt_file.write(str(Whole_gcs_scores[i]))
            txt_file.write('\n')
            
    if len(new_features) > 0:
        X_train = np.concatenate([ X_train, new_features ], 0)
        y_train = np.concatenate([ y_train, np.array(scores)[ :, None ] ], 0)

    iteration += 1
    
    print(iteration)
