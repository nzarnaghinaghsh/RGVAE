import numpy as np
import matplotlib.pyplot as plt

fname = 'Whole_gcs0.txt'

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
plt.title('Histogram of the Whole_gcs_scores0 of the generated sequences')
plt.savefig('Whole_gcs_scores0.png', dpi=500)





fname = 'Whole_gcs9.txt'

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
plt.title('Histogram of the Whole_gcs_scores9 of the generated sequences')
plt.savefig('Whole_gcs_scores9.png', dpi=500)







fname = 'Whole_length9.txt'

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
plt.title('Histogram of the Whole_length9 of the generated sequences')
plt.savefig('Whole_length9.png', dpi=500)




fname = 'Whole_MFE_scores9.txt'

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
plt.title('Histogram of the Whole_MFE_scores9 of the generated sequences')
plt.savefig('Whole_MFE_scores9.png', dpi=500)