import pickle
import gzip

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
    
    
    Whole_rnas = []
    Whole_target_struct_dist1 = []
    Whole_target_struct_dist2 = []
    Whole_target_struct_dist3 = []
    
    
    #final_rnas.append(rna_seq)

    target_structure1 = "..((((((............((((.((((........((....)).((((......)))))))).))))..(((((.....)))))...)))))).."
    target_structure2 = "..((((((....(((((...(((.....)))......)))))(((.((((......)))))))..(((...(((((.....)))))))))))))).."
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

    for i in range(latent_points.shape[ 0 ]):
        #rna_seq = "A"
        #structure, energy = RNA.fold(rna_seq)
        #min_energy = energy
        min_H_dist = 100000   # Minimum Hamming distance
        min_H_dist2 = 100000
        min_energy0 = 0
        min_energy2 = 0
        min_energy = 0
        
        for j in range(decode_attempts):
            if (( ~np.equal(RNA_series[ j, i ], None)) and (len(RNA_series[ j, i ])<200)):
                structure, energy = RNA.fold(RNA_series[ j, i ])
                Whole_rnas.append(RNA_series[ j, i ])
                # Computing the hamming distance between the target sequence and the structure:
                
                fc1 = RNA.fold_compound(RNA_series[ j, i ])
                (propensity1, ensemble_energy1) = fc1.pf()
                basepair_probs1 = fc1.bpp()
                dim1 = len(basepair_probs1[0])
                #base_profile1 = Make_bp_profile_bppm(basepair_probs1, dim1)
                
    
                #dist0, alignment0 = profile_edit_dist(base_profile0,base_profile1)
                #dist02, alignment02 = profile_edit_dist(base_profile02,base_profile1)
                
                window1, distance1 = align_and_distance(basepair_probs0,basepair_probs1)
                window2, distance2 = align_and_distance(basepair_probs02,basepair_probs1)
                Whole_target_struct_dist1.append(distance1)
                Whole_target_struct_dist2.append(distance2)
                
                '''
                len_s = len(structure)
                min_len = min(len_t,len_s) # Min length
                max_len = max(len_t,len_s) # Max length

                H1_dist = 0 # Hamming distance
                for t in range(min_len):
                    if (target_structure1[t] != structure[t]):
                        H1_dist = H1_dist + 1
    
                H1_dist = H1_dist + (max_len - min_len)
                # print("H_dist", H_dist)
                #target_structure_scores.append(H_dist)
                #print(i)
                
                H2_dist = 0 # Hamming distance
                for t in range(min_len):
                    if (target_structure2[t] != structure[t]):
                        H2_dist = H2_dist + 1
    
                H2_dist = H2_dist + (max_len - min_len)
                


                if (min(H1_dist,H2_dist) <min_H_dist):
                    min_H_dist = min(H1_dist,H2_dist)
                    rna_seq = RNA_series[ j, i ]
                    min_energy = energy
                    
                '''
                if ((distance1 <min_H_dist) or (distance2 <min_H_dist2)):
                    min_H_dist = distance1
                    min_H_dist2 = distance2
                    rna_seq = RNA_series[ j, i ]
                    min_energy = energy
                ''' 
                if (dist02 <min_H_dist2):
                    min_H_dist2 = dist02
                    rna_seq = RNA_series[ j, i ]
                    min_energy = energy
                '''

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
        '''
        g_c_scores.append(abs(RNA_vae.g_c(tokens)-0.5))
        temp = rna_seq
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
        #if (rna_seq[0:4]=="CUCGA" and rna_seq[21:31]=="UACAGAAAUUC"):
        #    motif_scores.append(0)
        #else:
        #    motif_scores.append(1)
        #print(i)
        '''
        
        print("rna_seq")
        print(rna_seq)
        print("min_H_dist")
        print(min_H_dist)
        final_rnas.append(rna_seq)
        MFE_scores.append(min_energy)
        target_structure_scores0.append(min_H_dist)
        target_structure_scores02.append(min_H_dist2)

    return final_rnas, MFE_scores, target_structure_scores0, target_structure_scores02, Whole_rnas, Whole_target_struct_dist1, Whole_target_struct_dist2

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
target_structure_scores0 = np.loadtxt('../../latent_features_and_targets_grammar/target_structure_scores0.txt')
target_structure_scores1 = np.loadtxt('../../latent_features_and_targets_grammar/target_structure_scores02.txt')

target_structure_scores_normalized0 = (np.array(target_structure_scores0) - np.mean(target_structure_scores0)) / np.std(target_structure_scores0)
target_structure_scores_normalized02 = (np.array(target_structure_scores1) - np.mean(target_structure_scores1)) / np.std(target_structure_scores1)

targets = target_structure_scores_normalized0 + target_structure_scores_normalized02

y = -targets
#y = -np.loadtxt('../../latent_features_and_targets_grammar/targets.txt')
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

    valid_RNA_final, MFE_final, target_structure_final0, target_structure_final02, Whole_rnas, Whole_target_struct_dist1, Whole_target_struct_dist2 = decode_from_latent_space(next_inputs, grammar_model)


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
    #target_structure_scores0 = np.loadtxt('../../latent_features_and_targets_grammar/target_structure_scores0.txt')
    #target_structure_scores02 = np.loadtxt('../../latent_features_and_targets_grammar/target_structure_scores02.txt')
    #MFE_scores = np.loadtxt('../../latent_features_and_targets_grammar/MFE_scores.txt')
    '''
    g_c_scores = np.loadtxt('../../latent_features_and_targets_grammar/g_c_scores.txt')
    motif_scores1 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores1.txt')
    motif_scores2 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores2.txt')
    motif_scores3 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores3.txt')
    motif_scores4 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores4.txt')
    motif_scores5 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores5.txt')
    motif_scores6 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores6.txt')
    motif_scores7 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores7.txt')
    motif_scores8 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores8.txt')
    motif_scores9 = np.loadtxt('../../latent_features_and_targets_grammar/motif_scores9.txt')
    '''
    #structure_scores = np.loadtxt('../../latent_features_and_targets_grammar/structure_scores.txt')
    #target_structure_scores_normalized0 = (np.array(target_structure_scores0) - np.mean(target_structure_scores0)) / np.std(target_structure_scores0)
    #target_structure_scores_normalized02 = (np.array(target_structure_scores1) - np.mean(target_structure_scores1)) / np.std(target_structure_scores1)
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
             
    
    #targets = target_structure_scores_normalized0 + target_structure_scores_normalized02 #+ MFE_scores_normalized   

    #import sascorer
    #import networkx as nx
    #from rdkit.Chem import rdmolops
    #target_structure = ".....((((((....))))))..........." #"(((((......)))))"
    #len_t = len(target_structure)

    scores = []
    #current_scores = []
    for i in range(len(valid_RNA_final)):
        if valid_RNA_final[ i ] is not None:
            '''
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
            '''
            current_target_structure_scores0 = target_structure_final0 [i]
            current_target_structure_scores02 = target_structure_final02 [i]
            current_MFE_scores = MFE_final[i]
            '''
            current_g_c_scores = g_c_final[i] 
            current_motif_scores1 = motif_final1[i]
            current_motif_scores2 = motif_final2[i]
            current_motif_scores3 = motif_final3[i]
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
         
            current_target_structure_scores_normalized0 = (np.array(current_target_structure_scores0) - np.mean(target_structure_scores0)) / np.std(target_structure_scores0)
            current_target_structure_scores_normalized02 = (np.array(current_target_structure_scores02) - np.mean(target_structure_scores1)) / np.std(target_structure_scores1)
            #current_MFE_scores_normalized = (np.array(current_MFE_scores) - np.mean(MFE_scores)) / np.std(MFE_scores)
            '''
            current_g_c_scores_normalized = (np.array(current_g_c_scores) - np.mean(g_c_scores)) / np.std(g_c_scores)
            current_motif_scores_normalized1 = (np.array(current_motif_scores1) - np.mean(motif_scores1)) / np.std(motif_scores1)
            current_motif_scores_normalized2 = (np.array(current_motif_scores2) - np.mean(motif_scores2)) / np.std(motif_scores2)
            current_motif_scores_normalized3 = (np.array(current_motif_scores3) - np.mean(motif_scores3)) / np.std(motif_scores3)
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
            score = current_target_structure_scores_normalized0 + current_target_structure_scores_normalized02#+current_MFE_scores_normalized
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

    with open("results/Whole_target_struct_dist1{}.txt".format(iteration), 'w') as txt_file:  # Specify encoding here
        for i in range(len(Whole_target_struct_dist1)):
            txt_file.write(str(Whole_target_struct_dist1[i]))
            txt_file.write('\n')
            
    with open("results/Whole_target_struct_dist2{}.txt".format(iteration), 'w') as txt_file:  # Specify encoding here
        for i in range(len(Whole_target_struct_dist2)):
            txt_file.write(str(Whole_target_struct_dist2[i]))
            txt_file.write('\n')

    if len(new_features) > 0:
        X_train = np.concatenate([ X_train, new_features ], 0)
        y_train = np.concatenate([ y_train, np.array(scores)[ :, None ] ], 0)

    iteration += 1
    
    print(iteration)
