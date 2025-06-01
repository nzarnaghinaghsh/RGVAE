import numpy as np
import matplotlib.pyplot as plt

fname = 'target_structure_scores0.txt'

with open(fname) as f:
    RNA_data1 = f.readlines()
    
for i in range(len(RNA_data1)):
    RNA_data1[ i ] = RNA_data1[ i ].strip()

fseq = 'trna200maxlen_102k_mod.txt'

with open(fseq) as f:
    RNA_seq = f.readlines()
    
for i in range(len(RNA_seq)):
    RNA_seq[ i ] = RNA_seq[ i ].strip()


fname = 'target_structure_scores02.txt'

with open(fname) as f:
    RNA_data2 = f.readlines()
  
RNA_data = []
for i in range(len(RNA_data2)):
    RNA_data2[ i ] = RNA_data2[ i ].strip()
    RNA_data.append(float(RNA_data1[ i ]) + float(RNA_data2[ i ]))




RNA_align = []
for i in range(len(RNA_data)):
    RNA_align.append((RNA_data[i]))
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
plt.title('Histogram of the summation of the alignment \n scores of the training sequences')
plt.savefig('hist_2nd_align_mod2_train.png', dpi=500)



num = 101754
RNA_align_sorted = sorted(RNA_align)
plt.figure()  # Create a new figure
plt.hist(RNA_align_sorted[0:10175], bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the summation of the alignment scores \n of the top %10 of the training sequences')
plt.savefig('hist_2nd_align_mod2_train_top10p.png', dpi=500)


num = 101754
RNA_align_sorted = sorted(RNA_align)
plt.figure()  # Create a new figure
plt.hist(RNA_align_sorted[0:1017], bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the summation of the alignment scores \n of the top %1 of the training sequences')
plt.savefig('hist_2nd_align_mod_train_top1p.png', dpi=500)


RNA_align_len = []
for i in range(len(RNA_align)):
    if (len(RNA_seq[i])>80 and len(RNA_seq[i])<110):
        RNA_align_len.append((RNA_align[i]))
    
plt.figure()  # Create a new figure
plt.hist(RNA_align_len, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the summation of the alignment \n scores of the training sequences with the length \n between 80 and 110.')
plt.savefig('hist_2nd_align_mod2_train_len80and110.png', dpi=500) 
    




RNA_align_len = []
for i in range(len(RNA_align)):
    if (len(RNA_seq[i])>80 and len(RNA_seq[i])<100):
        RNA_align_len.append((RNA_align[i]))
    
plt.figure()  # Create a new figure
plt.hist(RNA_align_len, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the summation of the alignment \n scores of the training sequences with the length \n between 80 and 100.')
plt.savefig('hist_2nd_align_mod2_train_len80and100.png', dpi=500)



RNA_align_len = []
for i in range(len(RNA_align)):
    if (len(RNA_seq[i])>90 and len(RNA_seq[i])<100):
        RNA_align_len.append((RNA_align[i]))
    
plt.figure()  # Create a new figure
plt.hist(RNA_align_len, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the summation of the alignment \n scores of the training sequences with the length \n between 90 and 100.')
plt.savefig('hist_2nd_align_mod2_train_len90and100.png', dpi=500)
