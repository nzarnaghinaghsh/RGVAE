import numpy as np
import matplotlib.pyplot as plt


'''
fname = 'Whole_rnas4.txt'

with open(fname) as f:
    RNA_seq = f.readlines()
    
for i in range(len(RNA_seq)):
    RNA_seq[ i ] = RNA_seq[ i ].strip()


'''


fname = 'Whole_MFE_scores_t.txt'

with open(fname) as f:
    RNA_data = f.readlines()
    
for i in range(len(RNA_data)):
    RNA_data[ i ] = RNA_data[ i ].strip()


RNA_gc = []
for i in range(len(RNA_data)):
    if (float(RNA_data[i]) != 0):
        RNA_gc.append(float(RNA_data[i]))


plt.figure()  # Create a new figure
plt.hist(RNA_gc, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the MFE score \n of the generated sequences')
plt.savefig('MFE_SLF4.png', dpi=500)


'''


fname = 'Whole_rnas1.txt'

with open(fname) as f:
    RNA_seq = f.readlines()
    
for i in range(len(RNA_seq)):
    RNA_seq[ i ] = RNA_seq[ i ].strip()





fname = 'Whole_MFE_scores1.txt'

with open(fname) as f:
    RNA_data = f.readlines()
    
for i in range(len(RNA_data)):
    RNA_data[ i ] = RNA_data[ i ].strip()


RNA_gc = []
for i in range(len(RNA_data)):
    if (len(RNA_seq[ i ])>1):
        RNA_gc.append(float(RNA_data[i]))


plt.figure()  # Create a new figure
plt.hist(RNA_gc, bins=50, color='red', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the MFE score \n of the generated sequences')
plt.savefig('MFE_SLF1.png', dpi=500)


fname = 'Whole_target_struct_dist2_t.txt'

with open(fname) as f:
    RNA_data = f.readlines()
    
for i in range(len(RNA_data)):
    RNA_data[ i ] = RNA_data[ i ].strip()


RNA_gc = []
for i in range(len(RNA_data)):
    if (len(RNA_seq[ i ])>1):
        RNA_gc.append(float(RNA_data[i]))


plt.figure()  # Create a new figure
plt.hist(RNA_gc, bins=50, color='red', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the second alignment distance \n of the generated sequences')
plt.savefig('Whole_target_struct_dist2_t.png', dpi=500)



fname = 'Whole_target_struct_dist3_t.txt'

with open(fname) as f:
    RNA_data = f.readlines()
    
for i in range(len(RNA_data)):
    RNA_data[ i ] = RNA_data[ i ].strip()


RNA_gc = []
for i in range(len(RNA_data)):
    if (len(RNA_seq[ i ])>1):
        RNA_gc.append(float(RNA_data[i]))


plt.figure()  # Create a new figure
plt.hist(RNA_gc, bins=50, color='red', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the sum of the alignment distances \n of the generated sequences')
plt.savefig('Whole_target_struct_dist3_t.png', dpi=500)
'''

