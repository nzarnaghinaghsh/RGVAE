import os
import numpy as np


# Directory containing the split files
directory = "./"  # Replace with your directory path if different





f1 = 'Whole_target_struct_dist1_tt.txt'
f2 = 'Whole_target_struct_dist2_tt.txt'
f3 = 'Whole_target_struct_dist3_tt.txt'
f4 = 'Whole_rnas_tt.txt'
f5 = 'Whole_gcs_scores_tt.txt'
f6 = 'Whole_forbidden_scores_tt.txt'

ft = "selected_results.txt"
i = 0

with open(f1) as f:
    Whole_target_struct_dist1_t = f.readlines()
for i in range(len(Whole_target_struct_dist1_t)):
    Whole_target_struct_dist1_t[ i ] = Whole_target_struct_dist1_t[ i ].strip()

with open(f2) as f:
    Whole_target_struct_dist2_t = f.readlines()
for i in range(len(Whole_target_struct_dist2_t)):
    Whole_target_struct_dist2_t[ i ] = Whole_target_struct_dist2_t[ i ].strip()

with open(f3) as f:
    Whole_target_struct_dist3_t = f.readlines()
for i in range(len(Whole_target_struct_dist3_t)):
    Whole_target_struct_dist3_t[ i ] = Whole_target_struct_dist3_t[ i ].strip()

with open(f4) as f:
    Whole_rnas = f.readlines()
for i in range(len(Whole_rnas)):
    Whole_rnas[ i ] = Whole_rnas[ i ].strip()
'''
with open(f5) as f:
    Whole_gcs_scores = f.readlines()
for i in range(len(Whole_gcs_scores)):
    Whole_gcs_scores[ i ] = Whole_gcs_scores[ i ].strip()

with open(f6) as f:
    Whole_forbidden_scores = f.readlines()
for i in range(len(Whole_forbidden_scores)):
    Whole_forbidden_scores[ i ] = Whole_forbidden_scores[ i ].strip()
'''

j1=0
j2=0
j3=0
min_struct_dist1 = 1000
min_struct_dist2 = 1000
min_struct_dist3= 1000
with open(ft, 'w') as outfile0:
    for i in range(len(Whole_target_struct_dist3_t)):
        if ((len(Whole_rnas[ i ])>1)):
            '''
                outfile0.write("RNA sequence: ")
                outfile0.write(Whole_rnas[ i ])
                outfile0.write("\n")
                outfile0.write("Whole_gcs[ i ]: ")
                outfile0.write(Whole_gcs[ i ])
                outfile0.write("\n")
                outfile0.write("Whole_motif_scores[ i ]: ")
                outfile0.write(Whole_motif_scores[ i ])
                outfile0.write("\n")
                outfile0.write("Whole_forbidden_scores[ i ]: ")
                outfile0.write(Whole_forbidden_scores[ i ])
                outfile0.write("\n")
                outfile0.write("Whole_MFE_scores[ i ]: ")
                outfile0.write(Whole_MFE_scores[ i ])
                outfile0.write("\n")
            '''
            if (min_struct_dist3 > float(Whole_target_struct_dist3_t[ i ])):
                min_struct_dist3 = float(Whole_target_struct_dist3_t[ i ])
                j3 = i
            if (min_struct_dist1 > float(Whole_target_struct_dist1_t[ i ])):
                min_struct_dist1 = float(Whole_target_struct_dist1_t[ i ])
                j1 = i
            if (min_struct_dist2 > float(Whole_target_struct_dist2_t[ i ])):
                min_struct_dist2 = float(Whole_target_struct_dist2_t[ i ])
                j2 = i
print("summation of the two scores:")
print("min struct_dist3 sequence:", Whole_rnas[ j3 ])
print("min_struct_dist3:", min_struct_dist3)

print("summation of the two scores:")
print("min struct_dist1 sequence:", Whole_rnas[ j1 ])
print("min_struct_dist1:", min_struct_dist1)

print("summation of the two scores:")
print("min struct_dist2 sequence:", Whole_rnas[ j2 ])
print("min_struct_dist2:", min_struct_dist2)

#print("shape", Whole_target_struct_dist3_t.shape())
'''
dist3_sorted = sorted(np.array(float(Whole_target_struct_dist3_t)), reverse=True)
dist2_sorted = sorted(np.array(float(Whole_target_struct_dist2_t)), reverse=True)
dist1_sorted = sorted(np.array(float(Whole_target_struct_dist1_t)), reverse=True)
'''

dist3_sorted = sorted([float(s) for s in Whole_target_struct_dist3_t])
dist2_sorted = sorted([float(s) for s in Whole_target_struct_dist2_t])
dist1_sorted = sorted([float(s) for s in Whole_target_struct_dist1_t])

print("dist3_sorted = ", dist3_sorted[0:10])
print("dist2_sorted = ", dist2_sorted[0:10])
print("dist1_sorted = ", dist1_sorted[0:10])




# Find indices where elements are equal
temp = [float(s) for s in Whole_target_struct_dist3_t]
print("temp =", temp[0])
print(np.array(temp).shape, np.array(dist3_sorted).shape)
print(np.array(temp).dtype, np.array(dist3_sorted).dtype)

for i in range(0,10):
    print("i = ", i)
    indices = np.where(np.array(temp) == np.array(dist3_sorted)[i])[0]

    print("indices",indices)  # Output: [0 2 4]
    print("RNA sequences: ", Whole_rnas[indices[0]])


#print("Whole_gcs:", Whole_gcs[ j ])
#print("Whole_motif_scores:", Whole_motif_scores[ j ])
#print("Whole_forbidden_scores:", Whole_forbidden_scores[ j ])
#print("Whole_MFE_scores:", Whole_MFE_scores[ j ])

#with open(f1, 'r') as infile1, open(f2, 'r') as infile2, open(f3, 'r') as infile3, open(f4, 'r') as infile4:
#    for line in infile1:
#    i = i+1
    
'''
with open(ft, 'w') as outfile0:
    with open(f1, 'r') as infile1:
        for line in infile1:
            outfile0.write(line)
    infile1.close()
    with open(f2, 'r') as infile2:
        for line in infile2:
            outfile0.write(line)
    infile2.close()
    with open(f3, 'r') as infile3:
        for line in infile3:
            outfile0.write(line)
    infile3.close()
    with open(f4, 'r') as infile4:
        for line in infile4:
            outfile0.write(line)
    infile4.close()
    with open(f5, 'r') as infile5:
        for line in infile5:
            outfile0.write(line)
    infile5.close()
    with open(f6, 'r') as infile6:
        for line in infile6:
            outfile0.write(line)
    infile6.close()
    with open(f7, 'r') as infile7:
        for line in infile7:
            outfile0.write(line)
    infile7.close()
    with open(f8, 'r') as infile8:
        for line in infile8:
            outfile0.write(line)
    infile8.close()
    with open(f9, 'r') as infile9:
        for line in infile9:
            outfile0.write(line)
    infile9.close()
    with open(f10, 'r') as infile10:
        for line in infile10:
            outfile0.write(line)
    infile10.close()
outfile0.close()

'''





print(f"All files combined")