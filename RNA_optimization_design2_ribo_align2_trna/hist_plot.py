import numpy as np
import matplotlib.pyplot as plt

fname = 'Whole_target_struct_dist3_tt.txt'
frna = 'Whole_rnas_tt.txt'

with open(frna) as f:
    RNA_seq = f.readlines()
    
for i in range(len(RNA_seq)):
    RNA_seq[ i ] = RNA_seq[ i ].strip()

with open(fname) as f:
    RNA_data = f.readlines()
    
for i in range(len(RNA_data)):
    RNA_data[ i ] = RNA_data[ i ].strip()

RNA_align = []
for i in range(len(RNA_data)):
    if (i<10):
        print("RNA_seq[i]",RNA_seq[i])
        print("len(RNA_seq[i])",len(RNA_seq[i]))
    if ((len(RNA_seq[i]) >1)):
        RNA_align.append(float(RNA_data[i]))
    '''
    if ((float(RNA_data[i]) !=0) and (float(RNA_data[i]) !=1) ):
        #print("i",i)
        #print("float(RNA_data[i])",float(RNA_data[i]))
        #print(float(RNA_data[i]))
        #print(float(RNA_data[i])<1)
        RNA_gc.append(float(RNA_data[i]))
    '''

plt.figure()  # Create a new figure
plt.hist(RNA_align, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the summation of the alignment scores \n of the generated sequences using RGVAE model')
plt.savefig('Whole_2nd_align_mod_RGVAE_len_g2.png', dpi=500)

num = 1869062
RNA_align_sorted = sorted(RNA_align)
plt.figure()  # Create a new figure
plt.hist(RNA_align_sorted[0:186906], bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the summation of the alignment scores \n of the top %10 generated sequences using RGVAE model')
plt.savefig('Whole_2nd_align_mod2_RGVAE_top10p.png', dpi=500)



num = 1869062
RNA_align_sorted = sorted(RNA_align)
plt.figure()  # Create a new figure
plt.hist(RNA_align_sorted[0:18690], bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the summation of the alignment scores \n of the top %1 generated sequences using RGVAE model')
plt.savefig('Whole_2nd_align_mod2_RGVAE_top1p.png', dpi=500)



RNA_align_len = []
for i in range(len(RNA_align)):
    if (len(RNA_seq[i])>80 and len(RNA_seq[i])<110):
        RNA_align_len.append(float(RNA_data[i]))
    
plt.figure()  # Create a new figure
plt.hist(RNA_align_len, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the summation of the alignment \n scores of the generated sequences using RGVAE with \n the length between 80 and 110.')
plt.savefig('hist_2nd_align_mod2_RGVAE_len80and110.png', dpi=500) 
 



RNA_align_len = []
for i in range(len(RNA_align)):
    if (len(RNA_seq[i])>80 and len(RNA_seq[i])<100):
        RNA_align_len.append(float(RNA_data[i]))
    
plt.figure()  # Create a new figure
plt.hist(RNA_align_len, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the summation of the alignment \n scores of the generated sequences using RGVAE with \n the length between 80 and 100.')
plt.savefig('hist_2nd_align_mod2_RGVAE_len80and100.png', dpi=500) 




RNA_align_len = []
for i in range(len(RNA_align)):
    if (len(RNA_seq[i])>90 and len(RNA_seq[i])<100):
        RNA_align_len.append(float(RNA_data[i]))
    
plt.figure()  # Create a new figure
plt.hist(RNA_align_len, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the summation of the alignment \n scores of the generated sequences using RGVAE with \n the length between 90 and 100.')
plt.savefig('hist_2nd_align_mod2_RGVAE_len90and100.png', dpi=500) 