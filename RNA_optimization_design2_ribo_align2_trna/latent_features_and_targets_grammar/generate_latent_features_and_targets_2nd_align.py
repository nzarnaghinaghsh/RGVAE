# Computing the base pair probability matrix and the edit distance to the alignment of the target structure.
#from rdkit.Chem import Descriptors
#from rdkit.Chem import rdmolops

#import sascorer

import sys
#sys.path.insert(0, '..')
##import RNA_vae
import numpy as np
import h5py
import RNA
#import perl
#from perlfunc import perl5lib, perlfunc, perlreq
import math

import sys
sys.path.append("/usr/local/ViennaRNA")

#! /usr/bin/env python

#from psfile import EPSFile
#from psfile import PSFile




def compute_distance(mat1, mat2):
    """Compute Euclidean distance between upper triangles of two matrices."""
    """It computes the distance of the upper triangle."""
    mat1 = np.array(mat1)
    mat2 = np.array(mat2)
    n = min(mat1.shape[0], mat2.shape[0])
    triu_indices = np.triu_indices(n, k=1)
    vec1 = mat1[triu_indices]
    vec2 = mat2[triu_indices]
    #vec1 = mat1[:n, :n][triu_indices]
    #vec2 = mat2[:n, :n][triu_indices]
    return np.linalg.norm(vec1 - vec2)

def sliding_window_distance(binary_mat, prob_mat):
    """Align matrices of different sizes using a sliding window approach."""
    binary_mat2 = np.array(binary_mat)
    prob_mat2 = np.array(prob_mat)
    N = binary_mat2.shape[0]
    M = prob_mat2.shape[0]
    min_distance = np.inf
    window0 = np.zeros((N, N))
    #if (M>=N):
    #    window0 = binary_mat2
    #else:
    #    window0 = prob_mat2
    
    # Slide the smaller matrix over the larger one along the diagonal
    window = np.zeros((N, N),dtype=float)
    
    if (M>=N):
        for i in range(M - N + 1):
            window = np.copy(prob_mat2[i:i+N, i:i+N])
            
            
            distance = compute_distance(binary_mat2, window)
            if distance < min_distance:
                min_distance = distance
                window0 = window
    else:
        
        for i in range(N - M + 1):
            window = np.zeros((N, N),dtype=float)
            window[i:i+M, i:i+M] = prob_mat2[:, :]
            '''
            if (i<3 or i==(N-M)):
                print("i",i)
                print("window shape = ",np.shape(window))
                print("window type = ",type(window))
                print("window sum = ",np.sum(window))
                print("window = ",window)
                print("prob_mat2",prob_mat2)
                print("prob_mat2 shape = ",np.shape(prob_mat2))
                print("window prob_mat2 = ",type(prob_mat2))
                print("prob_mat2 sum = ",np.sum(prob_mat2))
            '''
            #window = binary_mat2[i:i+M, i:i+M]
            distance = compute_distance(window, binary_mat2)
            if distance < min_distance:
                min_distance = distance
                window0 = window
    
    temp = abs(M-N)
    return window0, min_distance + ((temp)*(M+N+1))/2 

def align_and_distance(mat_a, mat_b):
    """Universal alignment and distance function."""
    '''
    if np.array(mat_a).shape == np.array(mat_b).shape:
        return compute_distance(np.array(mat_a), np.array(mat_b))
    else:
    '''
    return sliding_window_distance(mat_a, mat_b)





def profile_edit_dist(T1,T2): #length is the length of each row

    length1 = T1[0]
    length2 = T2[0]
    i_point = [[0 for i in range(length2+1)] for j in range(length1+1)]
    j_point = [[0 for i in range(length2+1)] for j in range(length1+1)]
    distance = [[0 for i in range(length2+1)] for j in range(length1+1)]
    
    edit_backtrack = 1
    '''
    for i in range (0, length1+1):
        for i in range (0, length1+1):
            i_point[i][j] = 0
            j_point[i][j] = 0
    '''
    for i in range (1, length1+1):
        distance[i][0] = distance[i - 1][0] + PrfEditCost(i, 0, T1, T2)
        if (edit_backtrack):
            i_point[i][0] = i - 1
            j_point[i][0] = 0
    '''
    print("First")
    print("i_point")
    print(i_point)
    print("j_point")
    print(j_point)
    '''
    
    
    for j in range (1, length2+1):
        distance[0][j] = distance[0][j - 1] + PrfEditCost(0, j, T1, T2)
        if (edit_backtrack):
            i_point[0][j] = 0
            j_point[0][j] = j - 1
    '''
    print("second")
    print("i_point")
    print(i_point)
    print("j_point")
    print(j_point)
    '''
    
    
    for i in range (1, length1+1):
        for j in range (1, length2+1):
            minus   = distance[i - 1][j] + PrfEditCost(i, 0, T1, T2)
            plus    = distance[i][j - 1] + PrfEditCost(0, j, T1, T2)
            change  = distance[i - 1][j - 1] + PrfEditCost(i, j, T1, T2)
            distance[i][j] = min(minus, min(plus, change))
            if (edit_backtrack):
                if (distance[i][j] == change):
                    i_point[i][j] = i - 1
                    j_point[i][j] = j - 1
                elif (distance[i][j] == plus):
                    i_point[i][j] = i
                    j_point[i][j] = j - 1
                else:
                    i_point[i][j] = i - 1
                    j_point[i][j] = j
                    
    
    temp = distance[length1][length2]
    
    alignment = [[0 for i in range(length1+length2+1)] for j in range(2)]
    
    #print("i_point")
    #print(i_point)
    #print("j_point")
    #print(j_point)

    if (edit_backtrack):
        pos = length1 + length2
        i   = length1
        j   = length2
        #print("lengths")
        #print(length1)
        #print(length2)
        while ((i > 0) or (j > 0)):
            i1  = i_point[i][j]
            j1  = j_point[i][j]
            if (((i - i1) == 1) and ((j - j1) == 1)):
                #substituition:
                alignment[0][pos] = i
                alignment[1][pos] = j
            if (((i - i1) == 1) and (j == j1)):
                #Deletion in [1]:
                alignment[0][pos] = i
                alignment[1][pos] = 0
            if ((i == i1) and ((j - j1) == 1)):
                #Deletion in [0]:
                alignment[0][pos] = 0
                alignment[1][pos] = j
            pos-=1
            #print("pos")
            #print(pos)
            #print(length1)
            #print(length2)
            #print("i1")
            #print(i1)
            #print("j1")
            #print(j1)
            
            i = i1
            j = j1
            
        for i in range(pos+1, length1 + length2+1):
            alignment[0][i - pos] = alignment[0][i]
            alignment[1][i - pos] = alignment[1][i]
        alignment[0][0] = length1 + length2 - pos;   #length of alignment
        
    
    
    return temp, alignment


def PrfEditCost(i, j, T1, T2):
    kmax = T1[1]
    if (T2[1] != kmax):
        print("inconsistent Profiles in PrfEditCost")
    if (i == 0):
        dist = 0
        for k in range (kmax):
            dist += T2[j*kmax + k]
    if (j == 0):
        dist = 0
        for k in range (kmax):
            dist += T1[i*kmax + k]
    if ((i > 0) and (j > 0)):
        dist = 2
        for k in range (kmax):
            #print("error")
            #print(T1[i * kmax + k])
            #print(T2[j * kmax + k])
            dist -= 2*math.sqrt(T1[i * kmax + k]* T2[j * kmax + k])
    return dist
        

    
def Make_bp_profile_bppm(bppm, length):
    L = 3 
    #P is the profile
    #P[i*3+0] unpaired, P[i*3+1] upstream, P[i*3+2] downstream p
    #indices start at 1 use first entries to store length and dimension
    P = [0 for i in range((length+2)*L)]
    #print("length of P")
    #print(len(P))
    #print(length)
    P[0]  = length
    P[1]  = L
    for i in range (1,length):
        for j in range (i+1,length+1):
            #print("P[i * L + 1]")
            #print(P[i * L + 1])
            #print("bppm[i][j]")
            #print(bppm[i][j])
            #print(i)
            #print(j)
            
            if (j==length):
                P[j * L + 2]  += 0
                P[j * L + 2] += 0
            else:
                #print("i")
                #print(i)
                #print("j")
                #print(j)
                P[j * L + 2]  += bppm[i][j]
                P[i * L + 1]  += bppm[i][j]
    for i in range (1,length+1):
        P[i * 3 + 0] = 1 - P[i * 3 + 1] - P[i * 3 + 2]
    
    return P






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
seq = "AACUUAUCAAGAGAAGUGGAGGGACUGGCCCAAAGAAGCUUCGGCAACAUUGUAUCAUGUGCCAAUUCCAGUAACCGAGAAGGUUAGAAGAUAAGGU"
structure0, energy0 = RNA.fold(seq)
print("{} [ {:6.2f} ]".format(structure0, energy0))

'''
fc0 = RNA.fold_compound(seq)
(propensity0, ensemble_energy0) = fc0.pf()
basepair_probs0 = fc0.bpp()
dim0 = len(basepair_probs0[0])

base_profile0 = Make_bp_profile_bppm(basepair_probs0, dim0)
'''
basepair_probs0 = [[0 for i in range(97)] for j in range(97)]
basepair_probs02 = [[0 for i in range(97)] for j in range(97)]

basepair_probs0[3][95] = 1
basepair_probs0[4][94] = 1
basepair_probs0[5][93] = 1
basepair_probs0[6][92] = 1
basepair_probs0[7][91] = 1
basepair_probs0[8][90] = 1

basepair_probs0[21][69] = 1
basepair_probs0[22][68] = 1
basepair_probs0[23][67] = 1
basepair_probs0[24][66] = 1
basepair_probs0[26][64] = 1
basepair_probs0[27][63] = 1
basepair_probs0[28][62] = 1
basepair_probs0[29][61] = 1


basepair_probs0[38][45] = 1
basepair_probs0[39][44] = 1
basepair_probs0[48][59] = 1
basepair_probs0[49][58] = 1
basepair_probs0[50][57] = 1



basepair_probs0[72][86] = 1
basepair_probs0[73][85] = 1
basepair_probs0[74][84] = 1
basepair_probs0[75][83] = 1
basepair_probs0[76][82] = 1



basepair_probs02[3][95] = 1
basepair_probs02[4][94] = 1
basepair_probs02[5][93] = 1
basepair_probs02[6][92] = 1
basepair_probs02[7][91] = 1
basepair_probs02[8][90] = 1

basepair_probs02[13][42] = 1
basepair_probs02[14][41] = 1
basepair_probs02[15][40] = 1
basepair_probs02[16][39] = 1
basepair_probs02[17][38] = 1
basepair_probs02[18][37] = 1

basepair_probs02[21][31] = 1
basepair_probs02[22][30] = 1
basepair_probs02[23][29] = 1

basepair_probs02[43][63] = 1
basepair_probs02[44][62] = 1
basepair_probs02[45][61] = 1
basepair_probs02[47][60] = 1
basepair_probs02[48][59] = 1
basepair_probs02[49][58] = 1
basepair_probs02[50][57] = 1

basepair_probs02[66][89] = 1
basepair_probs02[67][88] = 1
basepair_probs02[68][87] = 1
basepair_probs02[72][86] = 1
basepair_probs02[73][85] = 1
basepair_probs02[74][84] = 1
basepair_probs02[75][83] = 1

for i in range(97):
    for j in range(i,97):
        basepair_probs0[j][i] = basepair_probs0[i][j]
        basepair_probs02[j][i] = basepair_probs02[i][j]

dim0 = len(basepair_probs0[0])

#base_profile0 = Make_bp_profile_bppm(basepair_probs0, dim0)
#base_profile02 = Make_bp_profile_bppm(basepair_probs02, dim0)

'''
                     1234567890 1234567890 1234567890 1234567890 1234567890 1234567890 1234567890 1234567890 1234567890 1234567
target_structure1 = "..((((((.. .......... ((((.((((. .......((. ...)).(((( ......)))) )))).)))). .(((((.... .)))))...) ))))).."
target_structure2 = "..((((((.. ..(((((... (((.....)) )......))) ))(((.(((( ......)))) )))..(((.. .(((((.... .))))))))) ))))).."
'''
target_structure1 = "..((((((............((((.((((........((....)).((((......)))))))).))))..(((((.....)))))...)))))).."
target_structure2 = "..((((((....(((((...(((.....)))......)))))(((.((((......)))))))..(((...(((((.....)))))))))))))).."

len_t = len(target_structure1)

MFE_scores = []
structure_scores = []
target_structure_scores0 = []
target_structure_scores02 = []
motif_scores1 = []
motif_scores2 = []
motif_scores3 = []
motif_scores4 = []
motif_scores5 = []
motif_scores6 = []
motif_scores7 = []
motif_scores8 = []
motif_scores9 = []

for i in range(len(RNA_data)):
    """Compute the native structure (in dot bracket notation) of a one hot encoded sequence."""
    #if ((len(RNA_data[i])>10) and (len(RNA_data[i])<150)):
    #print("RNA_data[i]=", RNA_data[i])
    structure, energy = RNA.fold(RNA_data[i])
    #MFE_scores.append(energy)
    structure_scores.append(structure)
    #print("energy")
    #print(type(energy))
    
    fc1 = RNA.fold_compound(RNA_data[i])
    (propensity1, ensemble_energy1) = fc1.pf()
    basepair_probs1 = fc1.bpp()
    dim1 = len(basepair_probs1[0])

    #base_profile1 = Make_bp_profile_bppm(basepair_probs1, dim1)
    
    #dist0, alignment0 = profile_edit_dist(base_profile0,base_profile1)
    #dist02, alignment02 = profile_edit_dist(base_profile02,base_profile1)
    '''
    print("edit distance")
    print(dist)
    print("alignment")
    print(alignment)
    # Computing the hamming distance between the target sequence and the structure:
    len_s = len(structure)
    min_len = min(len_t,len_s) # Min length
    max_len = max(len_t,len_s) # Max length
    
    H1_dist = 0 # Hamming distance
    for j in range(min_len):
        if (target_structure1[j] != structure[j]):
            H1_dist = H1_dist + 1
    
    H1_dist = H1_dist + (max_len - min_len)
    
    H2_dist = 0 # Hamming distance
    for j in range(min_len):
        if (target_structure2[j] != structure[j]):
            H2_dist = H2_dist + 1
    
    H2_dist = H2_dist + (max_len - min_len)
    H_dist = min(H1_dist,H2_dist)
    # print("H_dist", H_dist)
    target_structure_scores.append(H_dist)
    print("H_dist")
    print(type(H_dist))
    '''
    basepair_probs11 = [list(row) for row in basepair_probs1]
    #print("type(basepair_probs0)",type(basepair_probs0))
    #print("type(basepair_probs11)",type(basepair_probs11))
    #print("basepair_probs11",basepair_probs11)
    #total = sum(sum(row) for row in basepair_probs11)
    #print("sum(basepair_probs11)",total)
    #print("type(basepair_probs02)",type(basepair_probs02))
    #print("align results)=",  align_and_distance(basepair_probs0,basepair_probs11))
    window1, distance1 = align_and_distance(basepair_probs0,basepair_probs11)
    window2, distance2 = align_and_distance(basepair_probs02,basepair_probs11)
    target_structure_scores0.append(distance1)
    target_structure_scores02.append(distance2)
    print(i)

#import networkx as nx

'''
g_c_scores = []
for i in range(len(RNA_data)):
    """Compute the native structure (in dot bracket notation) of a one hot encoded sequence."""
    tokens = RNA_vae.tokenization(RNA_data[i])
    g_c_scores.append(abs(RNA_vae.g_c(tokens)-0.5))
    print(i)
    temp = RNA_data[i]
    if (temp[0:1]=="CU"):
        motif_scores1.append(0)
    else:
        motif_scores1.append(1)
    if (temp[2:3]=="CG"):
        motif_scores2.append(0)
    else:
        motif_scores2.append(1)
    if (temp[4]=="A"):
        motif_scores3.append(0)
    else:
        motif_scores3.append(1)
    if (temp[21:22]=="UA"):
        motif_scores4.append(0)
    else:
        motif_scores4.append(1)
    if (temp[23:24]=="CA"):
        motif_scores5.append(0)
    else:
        motif_scores5.append(1)
    if (temp[25:26]=="GA"):
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

'''
g_c_scores_normalized = (np.array(g_c_scores) - np.mean(g_c_scores)) / np.std(g_c_scores)

motif_scores_normalized1 = (np.array(motif_scores1) - np.mean(motif_scores1)) / np.std(motif_scores1)
motif_scores_normalized2 = (np.array(motif_scores2) - np.mean(motif_scores2)) / np.std(motif_scores2)
motif_scores_normalized3 = (np.array(motif_scores3) - np.mean(motif_scores3)) / np.std(motif_scores3)
motif_scores_normalized4 = (np.array(motif_scores4) - np.mean(motif_scores4)) / np.std(motif_scores4)
motif_scores_normalized5 = (np.array(motif_scores5) - np.mean(motif_scores5)) / np.std(motif_scores5)
motif_scores_normalized6 = (np.array(motif_scores6) - np.mean(motif_scores6)) / np.std(motif_scores6)
motif_scores_normalized7 = (np.array(motif_scores7) - np.mean(motif_scores7)) / np.std(motif_scores7)
motif_scores_normalized8 = (np.array(motif_scores8) - np.mean(motif_scores8)) / np.std(motif_scores8)
motif_scores_normalized9 = (np.array(motif_scores9) - np.mean(motif_scores9)) / np.std(motif_scores9)
'''

target_structure_scores_normalized0 = (np.array(target_structure_scores0) - np.mean(target_structure_scores0)) / np.std(target_structure_scores0)
target_structure_scores_normalized02 = (np.array(target_structure_scores02) - np.mean(target_structure_scores02)) / np.std(target_structure_scores02)

#latent_points = grammar_model.encode(RNA_data)

# We store the results

#latent_points = np.array(latent_points)
#np.savetxt('latent_faetures.txt', latent_points)

targets = target_structure_scores_normalized0 + target_structure_scores_normalized02 #+  MFE_scores_normalized 
#np.savetxt('targets.txt', targets)
#np.savetxt('logP_values.txt', np.array(logP_values))
#np.savetxt('MFE_scores.txt', np.array(MFE_scores))
'''
np.savetxt('g_c_scores.txt', np.array(g_c_scores))
np.savetxt('motif_scores1.txt', np.array(motif_scores1))
np.savetxt('motif_scores2.txt', np.array(motif_scores2))
np.savetxt('motif_scores3.txt', np.array(motif_scores3))
np.savetxt('motif_scores4.txt', np.array(motif_scores4))
np.savetxt('motif_scores5.txt', np.array(motif_scores5))
np.savetxt('motif_scores6.txt', np.array(motif_scores6))
np.savetxt('motif_scores7.txt', np.array(motif_scores7))
np.savetxt('motif_scores8.txt', np.array(motif_scores8))
np.savetxt('motif_scores9.txt', np.array(motif_scores9))
'''
np.savetxt('target_structure_scores0.txt', np.array(target_structure_scores0))
np.savetxt('target_structure_scores02.txt', np.array(target_structure_scores02))
#np.savetxt('structure_scores.txt', "%s\n" %structure_scores)



