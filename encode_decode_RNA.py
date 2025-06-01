import sys
#sys.path.insert(0, '..')
import RNA_vae
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
    for i in range (0,length):
        for j in range (i,length):
            #print("P[i * L + 1]")
            #print(P[i * L + 1])
            #print("bppm[i][j]")
            #print(bppm[i][j])
            #print(i)
            #print(j)
            P[(i+1) * L + 1]  += bppm[i][j]
            P[(j+1) * L + 2]  += bppm[i][j]
    for i in range (1,length+1):
        P[i * 3 + 0] = 1 - P[i * 3 + 1] - P[i * 3 + 2]
    
    return P
    






d = RNA.duplexfold('AGAGGGAGCTGACGGATGACGCCCCCGCGCCACGCCGC', 'TATGCACGTGGC')
print("RNA duplex")
print(d.structure)
print(d.energy)
print(d.i)
print(d.j)


seq = "GGGAGAAUAGCUCAGUUGGGAGAGCGUUUGACUGAAGUCUACUUCGAUCAGCAACUACCAUUAGGUAAUCAGAAGGUCACCGGUUCGAUCCCGGUUUGUCCCA"
seq2 = "GAGGAAAGUCCCGCCUCCAGAUCAAGGGAAGUCCCGCGA"
structure, energy = RNA.fold(seq)
structure2, energy2 = RNA.fold(seq2)

print("{} [ {:6.2f} ]".format(structure, energy))
print("{} [ {:6.2f} ]".format(structure2, energy2))
# 1. load grammar VAE
grammar_weights = "pretrained/RNA_vae_grammar_h_conditional_switch100_c234_L10_E100_batchB.hdf5"
grammar_model = RNA_vae.RNAGrammarModel(grammar_weights)

# 2. let's encode and decode some example SMILES strings


fc = RNA.fold_compound(seq)
(propensity, ensemble_energy) = fc.pf()
basepair_probs = fc.bpp()
dim = len(basepair_probs)
print("dim")
print(dim)
print(type(basepair_probs))
#basepair_probs2 = basepair_probs

dim1 = len(basepair_probs[0])


fc2 = RNA.fold_compound(seq2)
(propensity2, ensemble_energy2) = fc2.pf()
basepair_probs2 = fc2.bpp()
dim2 = len(basepair_probs2)
print("dim2")
print(dim2)
print(type(basepair_probs2))


dim2 = len(basepair_probs2[0])



#fd = PSFile("basepair_probs1.ps", paper="letter")
#fd.write("")
# ... write PostScript commands to `fd` ...


'''

print("base-pairs matrix")
for i in range(len(basepair_probs)):
    dim1 = len(basepair_probs[i])
    print("dim1")
    print(dim1)
    for j in range(len(basepair_probs[i])):
        #print ("pr(%d,%d) = %g" % (i, j, basepair_probs[i][j]))
        #print("i")
        #print(i)
        #print("j")
        #print(j)
        print (i, j, basepair_probs[i][j])   # print probabilities
        #fd.write("basepair_probs[i][j]")
    #fd.write("\n")

#fd.close()


#fd2 = PSFile("basepair_probs2.ps", paper="letter")
#fd.write("")
# ... write PostScript commands to `fd` ...




print("base-pairs matrix")
for i in range(len(basepair_probs)):
    dim1 = len(basepair_probs[i])
    print("dim1")
    print(dim1)
    for j in range(len(basepair_probs[i])):
        #print ("pr(%d,%d) = %g" % (i, j, basepair_probs[i][j]))
        #print("i")
        #print(i)
        #print("j")
        #print(j)
        print (i, j, basepair_probs[i][j])  # print probabilities
        #fd2.write("basepair_probs[i][j]")
    #fd2.write("\n")

#fd2.close()
'''


base_profile1 = Make_bp_profile_bppm(basepair_probs, dim1)
base_profile2 = Make_bp_profile_bppm(basepair_probs2, dim2)
dist, alignment = profile_edit_dist(base_profile1,base_profile1)
print("edit distance")
print(dist)
print("alignment")
print(alignment)



'''
@perlfunc
@perlreq('pmcomp.pl')
def usage(basepair_probs1, basepair_probs2):
    pass  # Empty body

print(myfunc(1, 3))  # Should print 4
'''

#profile = *Make_bp_profile_bppm ( basepair_probs, dim)

#dist =  profile_edit_distance ( profile, profile)

#print("dist")
#print(dist)


'''
f = open('data/benchmark_edited.stk','r')  # The processed data
L = []

count = -1
for line in f:
    line = line.strip()
    L.append(line)
f.close()

smiles = list(L[1])

print(smiles)
'''


smiles = ["GAGGAAAGUCCCGCCUCCAGAUCAAGGGAAGUCCCGCGA",
          "GGGACAAGGGUAGUACCCUUGGCAACUGCACAGAAAACUU",
          "ACCCCUAAAUAUUCAAUGAGGAUUUGAUUCGACUCUUACC"]


# z: encoded latent points
# NOTE: this operation returns the mean of the encoding distribution
# if you would like it to sample from that distribution instead
# replace line 83 in molecule_vae.py with: return self.vae.encoder.predict(one_hot)
z1 = grammar_model.encode(smiles)

print(z1)

# mol: decoded SMILES string
# NOTE: decoding is stochastic so calling this function many
# times for the same latent point will return different answers



for mol,real in zip(grammar_model.decode(z1),smiles):
    print(mol) # + '  ' + real)


'''
# 3. the character VAE (https://github.com/maxhodak/keras-molecules)
# works the same way, let's load it
char_weights = "pretrained/zRNA_vae_grammar_h_conditional100_c234_L2_E100_batchB.hdf5"
char_model = molecule_vae.RNACharacterModel(char_weights)

# 4. encode and decode
z2 = char_model.encode(smiles)
for mol in char_model.decode(z2):
    print mol

'''









