import nltk
import numpy as np
import six
import pdb

# the RNA CFG grammar

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
grammar = """
S -> L S | L
L -> 'A' F 'U' | 'A' | 'U' F 'A' | 'U' | 'C' F 'G' | 'C' | 'G' F 'C' | 'G'
F -> 'A' F 'U' | 'U' F 'A' | 'C' F 'G' | 'G' F 'C' | L S
Nothing -> None
"""
'''

# Make a chartparser
cfg = nltk.CFG.fromstring(grammar)
