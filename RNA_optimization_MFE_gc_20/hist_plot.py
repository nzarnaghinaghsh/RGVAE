import numpy as np
import matplotlib.pyplot as plt



fname = 'Whole_rnas_tt.txt'

with open(fname) as f:
    RNA_seq = f.readlines()
    
for i in range(len(RNA_seq)):
    RNA_seq[ i ] = RNA_seq[ i ].strip()





fname = 'Whole_MFE_scores_tt.txt'

with open(fname) as f:
    RNA_data = f.readlines()
    
for i in range(len(RNA_data)):
    RNA_data[ i ] = RNA_data[ i ].strip()


RNA_gc = []
for i in range(len(RNA_data)):
    if (len(RNA_seq[ i ])>1):
        RNA_gc.append(float(RNA_data[i]))


plt.figure()  # Create a new figure
plt.hist(RNA_gc, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the MFE score of the generated \n sequences for the target GC-content of 20%')
plt.savefig('MFE_SLF_g20.png', dpi=500)



fname = 'Whole_gcs_tt.txt'

with open(fname) as f:
    RNA_data = f.readlines()
    
for i in range(len(RNA_data)):
    RNA_data[ i ] = RNA_data[ i ].strip()


RNA_gc = []
for i in range(len(RNA_data)):
    if (len(RNA_seq[ i ])>1):
        RNA_gc.append(float(RNA_data[i]))


plt.figure()  # Create a new figure
plt.hist(RNA_gc, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the GC-content of the generated \n sequences for the target GC-content of 20%')
plt.savefig('GC_SLF_g20.png', dpi=500)

'''
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
plt.hist(RNA_gc, bins=50, color='skyblue', edgecolor='black')
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

