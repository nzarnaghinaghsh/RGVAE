#from rdkit.Chem import Descriptors
#from rdkit.Chem import rdmolops
#Target length = [100,150]

#import sascorer

import numpy as np  
import RNA

# We load the RNA data

#fname = '../../data/trna_down_sample_U.stk'
fname = '../../data/trna200maxlen_102k_mod.txt'

with open(fname) as f:
    RNA_data = f.readlines()

for i in range(len(RNA_data)):
    RNA_data[ i ] = RNA_data[ i ].strip()

# We load the auto-encoder

import sys
sys.path.insert(0, '../../')
import RNA_vae
#grammar_weights = '../../pretrained/RNA_vae_grammar_h_conditional_switch100_c234_L10_E100_batchB.hdf5'
grammar_weights = '../../pretrained/RNA_vae_grammar_SLF_h_100k100_c234_L10_E100_batchB.hdf5'
grammar_model = RNA_vae.RNAGrammarModel(grammar_weights)

#from rdkit.Chem import MolFromSmiles, MolToSmiles
#from rdkit.Chem import Draw
#import image
#import copy
#import time

'''
smiles_rdkit = []
for i in range(len(smiles)):
    smiles_rdkit.append(MolToSmiles(MolFromSmiles(smiles[ i ])))
    print(i)

logP_values = []
for i in range(len(smiles)):
    logP_values.append(Descriptors.MolLogP(MolFromSmiles(smiles_rdkit[ i ])))
    print(i)
'''

#target_structure = "(((((((....(((...........)))((((((((..(((((((((((((((((((...(((((......))))).)))))).)))))))))))))..))))))))..)))))))"
#len_t = len(target_structure)

MFE_scores = []
structure_scores = []
length_scores = []
length = []
U_length = 150 #Upper length limit
D_length = 100 #Down length limit

target_structure_scores = []
for i in range(len(RNA_data)):
    """Compute the native structure (in dot bracket notation) of a one hot encoded sequence."""
    #structure, energy = RNA.fold(RNA_data[i])
    #MFE_scores.append(energy)
    #structure_scores.append(structure)
    len_R = len(RNA_data[i]) #length of the RNA sequence
    length.append(len_R)
    if ((len_R>(D_length-1)) and (len_R<(U_length+1))):
        length_scores.append(0)
    elif (len_R>U_length):
        length_scores.append((len_R - U_length))
    elif (len_R<D_length):
        length_scores.append((D_length - len_R))
    '''
    # Computing the hamming distance between the target sequence and the structure:
    len_s = len(structure)
    min_len = min(len_t,len_s) # Min length
    max_len = max(len_t,len_s) # Max length
    
    H_dist = 0 # Hamming distance
    for j in range(min_len):
        if (target_structure[j] != structure[j]):
            H_dist = H_dist + 1
    
    H_dist = H_dist + (max_len - min_len)
    # print("H_dist", H_dist)
    target_structure_scores.append(H_dist)
    print(i)
    '''

#import networkx as nx

g_c_scores = []
'''
for i in range(len(RNA_data)):
    """Compute the native structure (in dot bracket notation) of a one hot encoded sequence."""
    tokens = RNA_vae.tokenization(RNA_data[i])
    g_c_scores.append(abs(RNA_vae.g_c(tokens)-0.5))
    print(i)
'''

'''
cycle_scores = []
for i in range(len(smiles)):
    cycle_list = nx.cycle_basis(nx.Graph(rdmolops.GetAdjacencyMatrix(MolFromSmiles(smiles_rdkit[ i ]))))
    if len(cycle_list) == 0:
        cycle_length = 0
    else:
        cycle_length = max([ len(j) for j in cycle_list ])
    if cycle_length <= 6:
        cycle_length = 0
    else:
        cycle_length = cycle_length - 6
    cycle_scores.append(-cycle_length)
    print(i)

SA_scores_normalized = (np.array(SA_scores) - np.mean(SA_scores)) / np.std(SA_scores)
logP_values_normalized = (np.array(logP_values) - np.mean(logP_values)) / np.std(logP_values)
cycle_scores_normalized = (np.array(cycle_scores) - np.mean(cycle_scores)) / np.std(cycle_scores)
'''
MFE_scores_normalized = (np.array(MFE_scores) - np.mean(MFE_scores)) / np.std(MFE_scores)

g_c_scores_normalized = (np.array(g_c_scores) - np.mean(g_c_scores)) / np.std(g_c_scores)

length_scores_normalized = (np.array(length_scores) - np.mean(length_scores)) / np.std(length_scores)

#target_structure_scores_normalized = (np.array(target_structure_scores) - np.mean(target_structure_scores)) / np.std(target_structure_scores)

#latent_points = grammar_model.encode(RNA_data)

# We store the results

#latent_points = np.array(latent_points)
#np.savetxt('latent_faetures.txt', latent_points)

#targets = MFE_scores_normalized + g_c_scores_normalized + length_scores_normalized# + target_structure_scores_normalized 
#np.savetxt('targets.txt', targets)
#np.savetxt('logP_values.txt', np.array(logP_values))
#np.savetxt('MFE_scores.txt', np.array(MFE_scores))
#np.savetxt('g_c_scores.txt', np.array(g_c_scores))
np.savetxt('length_scores_t.txt', np.array(length_scores))
np.savetxt('length_t.txt', np.array(length))
#np.savetxt('target_structure_scores.txt', np.array(target_structure_scores))
#np.savetxt('structure_scores.txt', "%s\n" %structure_scores)



