import numpy as np
import matplotlib.pyplot as plt

'''
fname = 'g_c_scores_t.txt'

with open(fname) as f:
    RNA_data = f.readlines()
    
for i in range(len(RNA_data)):
    RNA_data[ i ] = RNA_data[ i ].strip()


RNA_gc = []
for i in range(len(RNA_data)):
    if ((float(RNA_data[i]) !=0) and (float(RNA_data[i]) !=1) ):
        #print("i",i)
        #print("float(RNA_data[i])",float(RNA_data[i]))
        #print(float(RNA_data[i]))
        #print(float(RNA_data[i])<1)
        RNA_gc.append(float(RNA_data[i]))


plt.figure()  # Create a new figure
plt.hist(RNA_gc, bins=50, color='red', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the g_c_scores_t of the training sequences')
plt.savefig('g_c_scores_t.png', dpi=500)

'''



fname = 'selected_results_just_MFE.txt'

with open(fname) as f:
    RNA_data = f.readlines()
    
for i in range(len(RNA_data)):
    RNA_data[ i ] = RNA_data[ i ].strip()


RNA_gc = []
for i in range(len(RNA_data)):
    if ((float(RNA_data[i]) !=0) and (float(RNA_data[i]) !=1) ):
        #print("i",i)
        #print("float(RNA_data[i])",float(RNA_data[i]))
        #print(float(RNA_data[i]))
        #print(float(RNA_data[i])<1)
        RNA_gc.append(float(RNA_data[i]))

plt.figure()  # Create a new figure
plt.hist(RNA_gc, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the MFE of the selected generated sequences with constraints')
plt.savefig('just_MFE_scores_t.png', dpi=500)






'''
fname = 'length_t.txt'

with open(fname) as f:
    RNA_data = f.readlines()
    
for i in range(len(RNA_data)):
    RNA_data[ i ] = RNA_data[ i ].strip()


RNA_gc = []
for i in range(len(RNA_data)):
    if ((float(RNA_data[i]) !=0) and (float(RNA_data[i]) !=1) ):
        #print("i",i)
        #print("float(RNA_data[i])",float(RNA_data[i]))
        #print(float(RNA_data[i]))
        #print(float(RNA_data[i])<1)
        RNA_gc.append(float(RNA_data[i]))

plt.figure()  # Create a new figure
plt.hist(RNA_gc, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the length_scores_t of the training sequences')
plt.savefig('length_scores_t.png', dpi=500)
'''


