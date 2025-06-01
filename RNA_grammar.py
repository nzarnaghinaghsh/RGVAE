import nltk
import numpy as np
import six
import pdb

# the RNA grammar
'''
gram = """
    S -> S S [0.099123]
    S -> 'G' S 'C' [0.187488]
    S -> 'C' S 'G' [0.143455]
    S -> 'G' S [0.100518]
    S -> 'A' S 'U' [0.112871]
    S -> 'U' S [0.079398]
    S -> 'A' [0.035166]
    S -> 'U' S 'A' [0.104104]
    S -> 'C' S [0.027396]
    S -> 'A' S [0.026599]
    S -> 'G' [0.024308]
    S -> 'C' [0.01056]
    S -> 'U' [0.049014]
    Nothing -> None  [1.0]
"""
'''
'''
gram = """S -> L S [0.529331]
L -> F [0.298078]
F -> 'G' F 'C' [0.289603]
F -> 'C' F 'G' [0.175691]
F -> L S [0.265841]
L -> 'G' [0.217019]
F -> 'A' F 'U' [0.09519]
F -> 'U' F 'A' [0.173675]
S -> L [0.470669]
L -> 'A' [0.150331]
L -> 'U' [0.258356]
L -> 'C' [0.076215]
Nothing -> None  [1.0]
"""
'''

#'''
gram = """S -> L S [0.5] | L [0.5]
L -> 'A' F 'U' [0.125] | 'A' [0.125] | 'U' F 'A' [0.125] | 'U' [0.125] | 'C' F 'G' [0.125] | 'C' [0.125] | 'G' F 'C' [0.125] | 'G' [0.125]
F -> 'A' F 'U' [0.2] | 'U' F 'A' [0.2] | 'C' F 'G' [0.2] | 'G' F 'C' [0.2] | L S  [0.2]
Nothing -> None  [1.0]
"""
#'''


'''
gram = """S -> S S [0.5] | 'A' S [0.5] | 'U' S [0.5] | 'C' S [0.5] | 'G' S [0.5] | 'A' [0.5] | 'U' [0.5] | 'C' [0.5] | 'G' [0.5] | 'A' X1 [0.5] | 'U' X2 [0.5] | 'C' X3 [0.5] | 'G' X4 [0.5]
X1 -> S 'U' [1] 
X2 -> S 'A' [1] 
X3 -> S 'G' [1] 
X4 -> S 'C' [1] 
Nothing -> None  [1.0]
"""
'''

'''
gram = """S -> L S | L
L -> 'A' F 'U' | 'A' | 'U' F 'A' | 'U' | 'C' F 'G' | 'C' | 'G' F 'C' | 'G'
F -> 'A' F 'U' | 'U' F 'A' | 'C' F 'G' | 'G' F 'C' | L S
Nothing -> Nones
"""
'''

'''
gram = """S -> L F
L -> 'A' L 'U' | 'A' 
F -> 'C' L 'G' | 'C' 
Nothing -> Nones
"""
'''

'''
gram = """S -> L F
L -> 'A' L 'U' | 'A' | 'U' L 'A' | 'U'
F -> 'C' L 'G' | 'C' | 'G' L 'C' | 'G'
Nothing -> Nones
"""
'''
'''
gram = """S -> 'A' S | 'U' S | 'C' S | 'G' S | 'A' S 'U' S | 'U' S 'A' S | 'C' S 'G' S | 'G' S 'C' S | Nones
Nothing -> Nones
"""
'''


'''
gram = """S -> LS
S -> L
LS -> L
LS -> S
L -> AFU
L -> UFA
L -> GFC
L -> CFG
L -> 'A'
L -> 'U'
L -> 'C'
L -> 'G'
F -> AFU
F -> UFA
F -> GFC
F -> CFG
F -> LS
AFU -> 'A'
AFU -> F
AFU -> 'U'
UFA -> 'U'
UFA -> F
UFA -> 'A'
GFC -> 'G'
GFC -> F
GFC -> 'C'
CFG -> 'C'
CFG -> F
CFG -> 'G'
Nothing -> Nones
"""
'''

# form the CFG and get the start symbol
GCFG = nltk.PCFG.fromstring(gram)
start_index = GCFG.productions()[0].lhs()

# collect all lhs symbols (non-terminal), and the unique set of them  (lhs: left hand side)
all_lhs = [a.lhs().symbol() for a in GCFG.productions()]
lhs_list = []
for a in all_lhs:
    if a not in lhs_list:
        lhs_list.append(a)

D = len(GCFG.productions())

# this map tells us the rhs (terminal) symbol indices for each production rule (rhs: right hand side)
rhs_map = [None]*D
count = 0
for a in GCFG.productions():
    rhs_map[count] = []
    for b in a.rhs():
        if not isinstance(b,six.string_types):
            s = b.symbol()
            rhs_map[count].extend(list(np.where(np.array(lhs_list) == s)[0]))
    count = count + 1

masks = np.zeros((len(lhs_list),D))
count = 0

# this tells us for each lhs (non-terminal) symbol which productions rules should be masked (Decoding)
# masks is a matrix. Each row of it corresponds to a production rule and must be a logit vector.
# A one in this vector corresponds a specific rule and shows the index of that rule. 
# Each rule corresponds to a logit vector.
for sym in lhs_list:
    is_in = np.array([a == sym for a in all_lhs], dtype=int).reshape(1,-1)
    masks[count] = is_in
    count = count + 1

# this tells us the indices where the masks are equal to 1
index_array = []
for i in range(masks.shape[1]):
    index_array.append(np.where(masks[:,i]==1)[0][0])
ind_of_ind = np.array(index_array)

max_rhs = max([len(l) for l in rhs_map])

# rules 29 and 31 aren't used in the zinc data so we 
# 0 their masks so they can never be selected (in RNA, we do not need it)
# masks[:,29] = 0
# masks[:,31] = 0
