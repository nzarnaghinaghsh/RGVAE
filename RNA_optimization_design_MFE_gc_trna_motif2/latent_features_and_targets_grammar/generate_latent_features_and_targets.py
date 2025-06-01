#from rdkit.Chem import Descriptors
#from rdkit.Chem import rdmolops

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

target_structure = "......(((....)))"
len_t = len(target_structure)

MFE_scores = []
structure_scores = []
target_structure_scores = []
motif_scores1 = []
motif_scores2 = []
motif_scores3 = []
'''
motif_scores4 = []
motif_scores5 = []
motif_scores6 = []
motif_scores7 = []
motif_scores8 = []
motif_scores9 = []
'''

for i in range(len(RNA_data)):
    """Compute the native structure (in dot bracket notation) of a one hot encoded sequence."""
    structure, energy = RNA.fold(RNA_data[i])
    #MFE_scores.append(energy)
    structure_scores.append(structure)
    print("energy")
    print(type(energy))
    
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
    print("H_dist")
    print(type(H_dist))
    print(i)

#import networkx as nx

g_c_scores = []
for i in range(len(RNA_data)):
    """Compute the native structure (in dot bracket notation) of a one hot encoded sequence."""
    tokens = RNA_vae.tokenization(RNA_data[i])
    #g_c_scores.append(abs(RNA_vae.g_c(tokens)-0.5))
    print(i)
    temp = RNA_data[i]
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
    if (temp[6]=="AA"):
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
#MFE_scores_normalized = (np.array(MFE_scores) - np.mean(MFE_scores)) / np.std(MFE_scores)

#g_c_scores_normalized = (np.array(g_c_scores) - np.mean(g_c_scores)) / np.std(g_c_scores)

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

#latent_points = grammar_model.encode(RNA_data)

# We store the results

#latent_points = np.array(latent_points)
#np.savetxt('latent_faetures.txt', latent_points)

#targets = target_structure_scores_normalized + MFE_scores_normalized + g_c_scores_normalized + motif_scores_normalized1 + motif_scores_normalized2 + motif_scores_normalized3 + motif_scores_normalized4 + motif_scores_normalized5 + motif_scores_normalized6 
#np.savetxt('targets.txt', targets)
#np.savetxt('logP_values.txt', np.array(logP_values))
#np.savetxt('MFE_scores.txt', np.array(MFE_scores))
#np.savetxt('g_c_scores.txt', np.array(g_c_scores))
np.savetxt('motif_scores1.txt', np.array(motif_scores1))
np.savetxt('motif_scores2.txt', np.array(motif_scores2))
np.savetxt('motif_scores3.txt', np.array(motif_scores3))
'''
np.savetxt('motif_scores4.txt', np.array(motif_scores4))
np.savetxt('motif_scores5.txt', np.array(motif_scores5))
np.savetxt('motif_scores6.txt', np.array(motif_scores6))
np.savetxt('motif_scores7.txt', np.array(motif_scores7))
np.savetxt('motif_scores8.txt', np.array(motif_scores8))
np.savetxt('motif_scores9.txt', np.array(motif_scores9))
'''
np.savetxt('target_structure_scores.txt', np.array(target_structure_scores))
#np.savetxt('structure_scores.txt', "%s\n" %structure_scores)



