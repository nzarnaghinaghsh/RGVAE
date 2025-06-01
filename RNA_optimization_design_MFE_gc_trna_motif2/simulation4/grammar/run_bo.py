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
    motif_scores1 = []
    motif_scores2 = []
    motif_scores3 = []
    
    
    Whole_rnas = []
    Whole_MFE_scores = []
    Whole_gcs = []
    Whole_gcs_scores = []
    Whole_motif_scores1 = []
    Whole_motif_scores2 = []
    Whole_motif_scores3 = []
    Whole_h_dist = []
    
    '''
    motif_scores4 = []
    motif_scores5 = []
    motif_scores6 = []
    motif_scores7 = []
    motif_scores8 = []
    motif_scores9 = []
    '''
    #final_rnas.append(rna_seq)

    target_structure = "......(((....)))" #".....((((((....))))))..........." #"(((((......)))))"
    len_t = len(target_structure)

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
                
                temp = RNA_series[ j, i ]
                if (len(temp)>=1 and temp[0]=="U"):
                    Whole_motif_scores1.append(0)
                else:
                    Whole_motif_scores1.append(1)
                if (len(temp)>=2 and temp[1]=="C"):
                    Whole_motif_scores2.append(0)
                else:
                    Whole_motif_scores2.append(1)
                if (len(temp)>=3 and temp[2]=="G"):
                    Whole_motif_scores3.append(0)
                else:
                    Whole_motif_scores3.append(1)
                
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
                Whole_h_dist.append(H_dist)


                if (H_dist<min_H_dist):
                    min_H_dist = H_dist
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
        temp = rna_seq
        if (temp[0]=="U"):
            motif_scores1.append(0)
        else:
            motif_scores1.append(1)
        if (temp[1]=="C"):
            motif_scores2.append(0)
        else:
            motif_scores2.append(1)
        if (temp[2]=="G"):
            motif_scores3.append(0)
        else:
            motif_scores3.append(1)
        '''
        if (temp[3]=="U"):
            motif_scores4.append(0)
        else:
            motif_scores4.append(1)
        if (temp[4]=="C"):
            motif_scores5.append(0)
        else:
            motif_scores5.append(1)
        if (temp[5]=="G"):
            motif_scores6.append(0)
        else:
            motif_scores6.append(1)
        if (temp[27:28]=="AA"):
            motif_scores7.append(0)
        else:
            motif_scores7.append(1)
        if (temp[29:30]=="UU"):
            motif_scores8.append(0)
        else:
            motif_scores8.append(1)
        if (temp[31]=="C"):
            motif_scores9.append(0)
        else:
            motif_scores9.append(1)
        '''
        #if (rna_seq[0:4]=="CUCGA" and rna_seq[21:31]=="UACAGAAAUUC"):
        #    motif_scores.append(0)
        #else:
        #    motif_scores.append(1)
        #print(i)
        
        print("rna_seq")
        print(rna_seq)
        print("min_H_dist")
        print(min_H_dist)
        final_rnas.append(rna_seq)
        MFE_scores.append(min_energy)

    return final_rnas, MFE_scores, g_c_scores, motif_scores1, motif_scores2, motif_scores3, Whole_rnas, Whole_MFE_scores, Whole_gcs, Whole_gcs_scores, Whole_motif_scores1, Whole_motif_scores2, Whole_motif_scores3, Whole_h_dist

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

# We load the data

X = np.loadtxt('../../latent_features_and_targets_grammar/latent_faetures_t.txt')
#y = -np.loadtxt('../../latent_features_and_targets_grammar/targets.txt')

target_structure_scores = np.loadtxt('../../latent_features_and_targets_grammar/target_structure_scores.txt')
MFE_scores = np.loadtxt('../../latent_features_and_targets_grammar/MFE_scores_t.txt')
g_c_scores = np.loadtxt('../../latent_features_and_targets_grammar/g_c_scores_t.txt')
motif_scores1 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores1.txt')
motif_scores2 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores2.txt')
motif_scores3 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores3.txt')
'''
motif_scores4 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores4.txt')
motif_scores5 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores5.txt')
motif_scores6 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores6.txt')
motif_scores7 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores7.txt')
motif_scores8 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores8.txt')
motif_scores9 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores9.txt')
'''

MFE_scores_normalized = (np.array(MFE_scores) - np.mean(MFE_scores)) / np.std(MFE_scores)

g_c_scores_normalized = (np.array(g_c_scores) - np.mean(g_c_scores)) / np.std(g_c_scores)

motif_scores_normalized1 = (np.array(motif_scores1) - np.mean(motif_scores1)) / np.std(motif_scores1)
motif_scores_normalized2 = (np.array(motif_scores2) - np.mean(motif_scores2)) / np.std(motif_scores2)
motif_scores_normalized3 = (np.array(motif_scores3) - np.mean(motif_scores3)) / np.std(motif_scores3)
'''
motif_scores_normalized4 = (np.array(motif_scores4) - np.mean(motif_scores4)) / np.std(motif_scores4)
motif_scores_normalized5 = (np.array(motif_scores5) - np.mean(motif_scores5)) / np.std(motif_scores5)
motif_scores_normalized6 = (np.array(motif_scores6) - np.mean(motif_scores6)) / np.std(motif_scores6)
motif_scores_normalized7 = (np.array(motif_scores7) - np.mean(motif_scores7)) / np.std(motif_scores7)
motif_scores_normalized8 = (np.array(motif_scores8) - np.mean(motif_scores8)) / np.std(motif_scores8)
motif_scores_normalized9 = (np.array(motif_scores9) - np.mean(motif_scores9)) / np.std(motif_scores9)
'''

target_structure_scores_normalized = (np.array(target_structure_scores) - np.mean(target_structure_scores)) / np.std(target_structure_scores)

targets = target_structure_scores_normalized + g_c_scores_normalized + motif_scores_normalized1 + motif_scores_normalized2 + motif_scores_normalized3 

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

alpha = 6
iteration = 0
while iteration < 15:

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
    #grammar_weights = '../../../pretrained/RNA_vae_grammar_h_conditional_switch100_c234_L10_E100_batchB.hdf5'
    grammar_weights = '../../../pretrained/RNA_vae_grammar_SLF_h_100k100_c234_L10_E100_batchB.hdf5'
    grammar_model = RNA_vae.RNAGrammarModel(grammar_weights)

    # We pick the next 50 inputs

    next_inputs = sgp.batched_greedy_ei(50, np.min(X_train, 0), np.max(X_train, 0))

    valid_RNA_final, MFE_final, g_c_final, motif_final1, motif_final2, motif_final3, Whole_rnas, Whole_MFE_scores, Whole_gcs, Whole_gcs_scores, Whole_motif_scores1, Whole_motif_scores2, Whole_motif_scores3, Whole_h_dist  = decode_from_latent_space(next_inputs, grammar_model)
    
    
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

    
    #structure_scores = np.loadtxt('../../latent_features_and_targets_grammar/structure_scores.txt')
    target_structure_scores_normalized = (np.array(target_structure_scores) - np.mean(target_structure_scores)) / np.std(target_structure_scores)
    MFE_scores_normalized = (np.array(MFE_scores) - np.mean(MFE_scores)) / np.std(MFE_scores)
    g_c_scores_normalized = (np.array(g_c_scores) - np.mean(g_c_scores)) / np.std(g_c_scores)
    motif_scores_normalized1 = (np.array(motif_scores1) - np.mean(motif_scores1)) / np.std(motif_scores1)
    motif_scores_normalized2 = (np.array(motif_scores2) - np.mean(motif_scores2)) / np.std(motif_scores2)
    motif_scores_normalized3 = (np.array(motif_scores3) - np.mean(motif_scores3)) / np.std(motif_scores3)
    '''
    motif_scores_normalized4 = (np.array(motif_scores4) - np.mean(motif_scores4)) / np.std(motif_scores4)
    motif_scores_normalized5 = (np.array(motif_scores5) - np.mean(motif_scores5)) / np.std(motif_scores5)
    motif_scores_normalized6 = (np.array(motif_scores6) - np.mean(motif_scores6)) / np.std(motif_scores6)
    motif_scores_normalized7 = (np.array(motif_scores7) - np.mean(motif_scores7)) / np.std(motif_scores7)
    motif_scores_normalized8 = (np.array(motif_scores8) - np.mean(motif_scores8)) / np.std(motif_scores8)
    motif_scores_normalized9 = (np.array(motif_scores9) - np.mean(motif_scores9)) / np.std(motif_scores9)
    '''
             
    
    #targets = target_structure_scores_normalized + MFE_scores_normalized + g_c_scores_normalized + alpha*( motif_scores_normalized1 + motif_scores_normalized2 + motif_scores_normalized3 + motif_scores_normalized4 + motif_scores_normalized5 + motif_scores_normalized6) #+ motif_scores_normalized7 + motif_scores_normalized8 + motif_scores_normalized9 
    targets = target_structure_scores_normalized + g_c_scores_normalized + motif_scores_normalized1 + motif_scores_normalized2 + motif_scores_normalized3 #+ motif_scores_normalized4 + motif_scores_normalized5 + motif_scores_normalized6) #+ motif_scores_normalized7 + motif_scores_normalized8 + motif_scores_normalized9 

    #import sascorer
    #import networkx as nx
    #from rdkit.Chem import rdmolops
    target_structure = "......(((....)))" #".....((((((....))))))..........." #"(((((......)))))"
    len_t = len(target_structure)

    scores = []
    #current_scores = []
    for i in range(len(valid_RNA_final)):
        if valid_RNA_final[ i ] is not None:
            structure, energy = RNA.fold(valid_RNA_final[ i ])
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
            current_MFE_scores = MFE_final[i]
            current_g_c_scores = g_c_final[i] 
            current_motif_scores1 = motif_final1[i]
            current_motif_scores2 = motif_final2[i]
            current_motif_scores3 = motif_final3[i]
            '''
            current_motif_scores4 = motif_final4[i]
            current_motif_scores5 = motif_final5[i]
            current_motif_scores6 = motif_final6[i]
            current_motif_scores7 = motif_final7[i]
            current_motif_scores8 = motif_final8[i]
            current_motif_scores9 = motif_final9[i]
            '''
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
         
            current_target_structure_scores_normalized = (np.array(current_target_structure_scores) - np.mean(target_structure_scores)) / np.std(target_structure_scores)
            current_MFE_scores_normalized = (np.array(current_MFE_scores) - np.mean(MFE_scores)) / np.std(MFE_scores)
            current_g_c_scores_normalized = (np.array(current_g_c_scores) - np.mean(g_c_scores)) / np.std(g_c_scores)
            current_motif_scores_normalized1 = (np.array(current_motif_scores1) - np.mean(motif_scores1)) / np.std(motif_scores1)
            current_motif_scores_normalized2 = (np.array(current_motif_scores2) - np.mean(motif_scores2)) / np.std(motif_scores2)
            current_motif_scores_normalized3 = (np.array(current_motif_scores3) - np.mean(motif_scores3)) / np.std(motif_scores3)
            '''
            current_motif_scores_normalized4 = (np.array(current_motif_scores4) - np.mean(motif_scores4)) / np.std(motif_scores4)
            current_motif_scores_normalized5 = (np.array(current_motif_scores5) - np.mean(motif_scores5)) / np.std(motif_scores5)
            current_motif_scores_normalized6 = (np.array(current_motif_scores6) - np.mean(motif_scores6)) / np.std(motif_scores6)
            current_motif_scores_normalized7 = (np.array(current_motif_scores7) - np.mean(motif_scores7)) / np.std(motif_scores7)
            current_motif_scores_normalized8 = (np.array(current_motif_scores8) - np.mean(motif_scores8)) / np.std(motif_scores8)
            current_motif_scores_normalized9 = (np.array(current_motif_scores9) - np.mean(motif_scores9)) / np.std(motif_scores9)
            '''
            
            #current_SA_score_normalized = (current_SA_score - np.mean(SA_scores)) / np.std(SA_scores)
            #current_log_P_value_normalized = (current_log_P_value - np.mean(logP_values)) / np.std(logP_values)
            #current_cycle_score_normalized = (current_cycle_score - np.mean(cycle_scores)) / np.std(cycle_scores)

            #score = (current_SA_score_normalized + current_log_P_value_normalized + current_cycle_score_normalized)
            score = current_target_structure_scores_normalized+current_g_c_scores_normalized+ alpha*(current_motif_scores_normalized1+current_motif_scores_normalized2+current_motif_scores_normalized3)#+current_motif_scores_normalized7+current_motif_scores_normalized8+current_motif_scores_normalized9
        else:
            score = -max(y)[ 0 ]

        scores.append(-score)
        print(i)

    print(valid_RNA_final)
    print(scores)

    save_object(scores, "results/scores{}.dat".format(iteration))

    with open("results/Whole_rnas{}.txt".format(iteration), 'w') as txt_file:  # Specify encoding here
        for i in range(len(Whole_rnas)):
            txt_file.write(str(Whole_rnas[i]))
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

    with open("results/Whole_motif_scores1{}.txt".format(iteration), 'w') as txt_file:  # Specify encoding here
        for i in range(len(Whole_motif_scores1)):
            txt_file.write(str(Whole_motif_scores1[i]))
            txt_file.write('\n')

    with open("results/Whole_motif_scores2{}.txt".format(iteration), 'w') as txt_file:  # Specify encoding here
        for i in range(len(Whole_motif_scores2)):
            txt_file.write(str(Whole_motif_scores2[i]))
            txt_file.write('\n')

    with open("results/Whole_motif_scores3{}.txt".format(iteration), 'w') as txt_file:  # Specify encoding here
        for i in range(len(Whole_motif_scores3)):
            txt_file.write(str(Whole_motif_scores3[i]))
            txt_file.write('\n')

    with open("results/Whole_h_dist{}.txt".format(iteration), 'w') as txt_file:  # Specify encoding here
        for i in range(len(Whole_h_dist)):
            txt_file.write(str(Whole_h_dist[i]))
            txt_file.write('\n')
  

    if len(new_features) > 0:
        X_train = np.concatenate([ X_train, new_features ], 0)
        y_train = np.concatenate([ y_train, np.array(scores)[ :, None ] ], 0)

    iteration += 1
    
    print(iteration)
