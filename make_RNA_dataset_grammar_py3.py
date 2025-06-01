from __future__ import print_function
import nltk
import pdb
# import zinc_grammar
import RNA_grammar
import numpy as np
import h5py
#import molecule_vae
import RNA_vae



#f = open('data/250k_rndm_zinc_drugs_clean.smi','r')
f = open('data/trna200maxlen_102k_mod.txt','r') #trna_down_sample_U.stk','r')  # The processed data
L = []

count = -1
for line in f:
    line = line.strip()
    L.append(line)
f.close()

'''
print("data")
print (L[0])
print (L[1])
print (L[2])
'''

MAX_LEN=199 #188 #64 #277
#NCHARS = len(zinc_grammar.GCFG.productions())
NCHARS = len(RNA_grammar.GCFG.productions())

def to_one_hot(smiles):
    """ Encode a list of smiles strings to one-hot vectors """
    assert type(smiles) == list
    # The grammar
    '''
    grammar = """
    S -> S S 
    S -> 'G' S 'C' 
    S -> 'C' S 'G' 
    S -> 'G' S 
    S -> 'A' S 'U' 
    S -> 'U' S 
    S -> 'A' 
    S -> 'U' S 'A' 
    S -> 'C' S 
    S -> 'A' S 
    S -> 'G' 
    S -> 'C' 
    S -> 'U' 
    Nothing -> None
    """
    '''
    
    #'''
    grammar = """
    S -> L S | L
    L -> 'A' F 'U' | 'A' | 'U' F 'A' | 'U' | 'C' F 'G' | 'C' | 'G' F 'C' | 'G'
    F -> 'A' F 'U' | 'U' F 'A' | 'C' F 'G' | 'G' F 'C' | L S
    Nothing -> None
    """
    #'''
    # Make a chartparser
    cfg = nltk.CFG.fromstring(grammar)
    
    prod_map = {}
    #for ix, prod in enumerate(zinc_grammar.GCFG.productions()):
    for ix, prod in enumerate(cfg.productions()):
        prod_map[prod] = ix
    
    #tokenize = molecule_vae.get_zinc_tokenizer(zinc_grammar.GCFG)
    #tokenize = RNA_vae.get_RNA_tokenizer(RNA_grammar.GCFG)
    #tokens = list(map(tokenize, smiles))
    
    tokens = [RNA_vae.tokenization(item) for item in smiles]
    print("tokens")
    print(tokens)
    
    # smiles = ['CCC']
    #print('smiles')
    #print(smiles)
    #print('tokenize')
    #print(tokenize)
    
    #print('tokens')
    #print(tokens)
    #print('smiles')
    print(smiles[0])
    #parser = nltk.ChartParser(zinc_grammar.GCFG)
    #parser = nltk.ChartParser(RNA_grammar.GCFG)
    parser = nltk.ViterbiParser(RNA_grammar.GCFG)
    print('parse')
    #print(parser)
    #print(np.shape(tokens))
    #l1, l2 = np.shape(tokens)
    #parse_trees = [parser.parse(t).next() for t in tokens]
    #print(type(tokens))
    #print("length")
    #print(len(tokens))
    
    #print(tokens)
    parse_trees = [next(parser.parse(t)) for t in tokens]
    #for itr in range(l1):
        
    productions_seq = [tree.productions() for tree in parse_trees]
    #print("product")
    #print("productions_seq")
    #print(productions_seq)
    #print(np.shape(productions_seq))
    
    #print(productions_seq[0])
    indices = [np.array([prod_map[prod] for prod in entry], dtype=int) for entry in productions_seq]
    print("indices")
    #print(indices[0])
    print(len(indices))
    one_hot = np.zeros((len(indices), MAX_LEN, NCHARS), dtype=np.float32)
    #print(indices[1])
    for i in range(len(indices)):
        num_productions = len(indices[i])
        print("num_productions")
        print(num_productions)
        #print("i")
        #print(i)
        #print("indices")
        #print(indices)
        #print(len(indices[i]))
        #print(np.shape(indices))
        one_hot[i][np.arange(num_productions),indices[i]] = 1.
        one_hot[i][np.arange(num_productions, MAX_LEN),-1] = 1.
    return one_hot, num_productions


#print("NCHARS")
#print(NCHARS)
# lim = len(L) # upper limit
DLim = 0 # Down limit
UpLim = len(L) # upper limit
lim = UpLim - DLim
gap = 1
OH = np.zeros((lim,MAX_LEN,NCHARS))
max_num_prod = 0
for i in range(DLim, UpLim, gap):
    print('Processing: i=[' + str(i) + ':' + str(i+gap) + ']')
    onehot, num_productions = to_one_hot(L[i:i+gap])
    if (num_productions>max_num_prod):
        max_num_prod = num_productions
    #print("one-hot")
    #print(np.shape(onehot))
    #print(np.shape(OH))
    OH[i:i+gap,:,:] = onehot

print("Max number of productions = ")
print(max_num_prod)

h5f = h5py.File('RNA_grammar_dataset_test_Py3.h5','w')
h5f.create_dataset('data', data=OH)
h5f.close()

