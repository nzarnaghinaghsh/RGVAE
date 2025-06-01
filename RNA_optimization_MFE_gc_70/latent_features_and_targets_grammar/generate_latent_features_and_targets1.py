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
#RNA_data = RNA_data2 [0:20000]
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


#target_structure = "(((((((....(((...........)))((((((((..(((((((((((((((((((...(((((......))))).)))))).)))))))))))))..))))))))..)))))))"
#len_t = len(target_structure)

MFE_scores = []
structure_scores = []

target_structure_scores = []
g_c_scores = []

for i in range(len(RNA_data)):
    #print("i",i)
    """Compute the native structure (in dot bracket notation) of a one hot encoded sequence."""
    # MFE:
    structure, energy = RNA.fold(RNA_data[i])
    MFE_scores.append(energy)
    structure_scores.append(structure)
    # GC_50:
    tokens = RNA_vae.tokenization(RNA_data[i])
    g_c_scores.append(abs(RNA_vae.g_c(tokens)-0.7))

print("MFE is calculated")

print("g_c_scores is calculated")

MFE_scores_normalized = (np.array(MFE_scores) - np.mean(MFE_scores)) / np.std(MFE_scores)

g_c_scores_normalized = (np.array(g_c_scores) - np.mean(g_c_scores)) / np.std(g_c_scores)

#target_structure_scores_normalized = (np.array(target_structure_scores) - np.mean(target_structure_scores)) / np.std(target_structure_scores)

#latent_points = grammar_model.encode(RNA_data)
print("latent points are calculated")
# We store the results

#latent_points = np.array(latent_points)
#np.savetxt('latent_faetures1.txt', latent_points)

targets = MFE_scores_normalized + g_c_scores_normalized # + target_structure_scores_normalized 
np.savetxt('targets1.txt', targets)
#np.savetxt('logP_values.txt', np.array(logP_values))
np.savetxt('MFE_scores1.txt', np.array(MFE_scores))
np.savetxt('g_c_scores1.txt', np.array(g_c_scores))
#np.savetxt('target_structure_scores.txt', np.array(target_structure_scores))
#np.savetxt('structure_scores.txt', "%s\n" %structure_scores)



