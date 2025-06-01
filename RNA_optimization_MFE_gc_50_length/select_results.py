import os

# Directory containing the split files
directory = "./"  # Replace with your directory path if different






f1 = 'Whole_gcs_tt.txt'
f2 = 'Whole_length_tt.txt'
f3 = 'Whole_MFE_scores_tt.txt'
f4 = 'Whole_rnas_tt.txt'
f5 = 'Whole_gcs_scores_tt.txt'

ft = "selected_results.txt"
i = 0

with open(f1) as f:
    Whole_gcs = f.readlines()
for i in range(len(Whole_gcs)):
    Whole_gcs[ i ] = Whole_gcs[ i ].strip()

with open(f2) as f:
    Whole_length = f.readlines()
for i in range(len(Whole_length)):
    Whole_length[ i ] = Whole_length[ i ].strip()

with open(f3) as f:
    Whole_MFE_scores = f.readlines()
for i in range(len(Whole_MFE_scores)):
    Whole_MFE_scores[ i ] = Whole_MFE_scores[ i ].strip()

with open(f4) as f:
    Whole_rnas = f.readlines()
for i in range(len(Whole_rnas)):
    Whole_rnas[ i ] = Whole_rnas[ i ].strip()

with open(f5) as f:
    Whole_gcs_scores = f.readlines()
for i in range(len(Whole_gcs_scores)):
    Whole_gcs_scores[ i ] = Whole_gcs_scores[ i ].strip()

j=0
min_MFE = 10
with open(ft, 'w') as outfile0:
    for i in range(len(Whole_gcs_scores)):
        if ((float(Whole_gcs_scores[i])<0.02) and (len(Whole_rnas[ i ])>1) and (int(Whole_length[ i ])>99) and (int(Whole_length[ i ])<151)):
                outfile0.write("RNA sequence: ")
                outfile0.write(Whole_rnas[ i ])
                outfile0.write("\n")
                outfile0.write("Whole_gcs[ i ]: ")
                outfile0.write(Whole_gcs[ i ])
                outfile0.write("\n")
                outfile0.write("Whole_length[ i ]: ")
                outfile0.write(Whole_length[ i ])
                outfile0.write("\n")
                outfile0.write("Whole_MFE_scores[ i ]: ")
                outfile0.write(Whole_MFE_scores[ i ])
                outfile0.write("\n")
                if (min_MFE > float(Whole_MFE_scores[ i ])):
                    min_MFE = float(Whole_MFE_scores[ i ])
                    j = i

print("min MFE sequence with constraints:", Whole_rnas[ j ])
print("Whole_gcs:", Whole_gcs[ j ])
print("Whole_length:", Whole_length[ j ])
print("Whole_MFE_scores:", Whole_MFE_scores[ j ])
print("min_MFE:", min_MFE)
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