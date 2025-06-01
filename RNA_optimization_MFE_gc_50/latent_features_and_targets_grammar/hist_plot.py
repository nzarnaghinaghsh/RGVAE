import numpy as np
import matplotlib.pyplot as plt

fname = 'g_c.txt'

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
plt.title('Histogram of the G-C content of the training sequences')
plt.savefig('GC_train.png', dpi=500)







fname = 'MFE_scores_t.txt'

with open(fname) as f:
    RNA_data = f.readlines()
    
for i in range(len(RNA_data)):
    RNA_data[ i ] = RNA_data[ i ].strip()


RNA_gc = []
MFE_min = 100
for i in range(len(RNA_data)):
    if ((float(RNA_data[i]) !=0) and (float(RNA_data[i]) !=1) ):
        #print("i",i)
        #print("float(RNA_data[i])",float(RNA_data[i]))
        #print(float(RNA_data[i]))
        #print(float(RNA_data[i])<1)
        RNA_gc.append(float(RNA_data[i]))
        if (float(RNA_data[i]) < MFE_min):
            MFE_min = float(RNA_data[i])

plt.figure()  # Create a new figure
plt.hist(RNA_gc, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Histogram of the MFE of the training sequences')
plt.savefig('MFE_train.png', dpi=500)


print("min MFE of the training data =", MFE_min)