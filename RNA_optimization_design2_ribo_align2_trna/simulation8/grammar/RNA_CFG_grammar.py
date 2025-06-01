import nltk
import numpy as np
import six
import pdb

# the RNA CFG grammar
grammar = """
S -> L S | L
L -> 'A' F 'U' | 'A' | 'U' F 'A' | 'U' | 'C' F 'G' | 'C' | 'G' F 'C' | 'G'
F -> 'A' F 'U' | 'U' F 'A' | 'C' F 'G' | 'G' F 'C' | L S
Nothing -> Nones
"""
# Make a chartparser
cfg = nltk.CFG.fromstring(grammar)
