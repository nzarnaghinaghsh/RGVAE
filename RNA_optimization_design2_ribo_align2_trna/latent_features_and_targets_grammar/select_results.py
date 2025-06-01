import os

# Directory containing the split files
directory = "./"  # Replace with your directory path if different





f1 = 'target_structure_scores0.txt'
f2 = 'target_structure_scores02.txt'
f3 = 'MFE_scores_t.txt'
f4 = 'trna200maxlen_102k_mod.txt'
f5 = 'g_c_scores_t.txt'
f6 = 'forbidden_scores_t.txt'

ft = "selected_results.txt"
i = 0

with open(f1) as f:
    target_structure_scores0 = f.readlines()
for i in range(len(target_structure_scores0)):
    target_structure_scores0[ i ] = target_structure_scores0[ i ].strip()

target_structure_scores03 = []
with open(f2) as f:
    target_structure_scores02 = f.readlines()
for i in range(len(target_structure_scores02)):
    target_structure_scores02[ i ] = target_structure_scores02[ i ].strip()
    target_structure_scores03.append(float(target_structure_scores02[i]) + float(target_structure_scores0[i]))




'''
with open(f3) as f:
    Whole_MFE_scores = f.readlines()
for i in range(len(Whole_MFE_scores)):
    Whole_MFE_scores[ i ] = Whole_MFE_scores[ i ].strip()
'''
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
min_struct_dist3 = 1000
with open(ft, 'w') as outfile0:
    for i in range(len(target_structure_scores03)):
        if ((len(Whole_rnas[ i ])>1)):
            if (min_struct_dist3 > float(target_structure_scores03[ i ])):
                min_struct_dist3 = float(target_structure_scores03[ i ])
                j3 = i
            if (min_struct_dist2 > float(target_structure_scores02[ i ])):
                min_struct_dist2 = float(target_structure_scores02[ i ])
                j2 = i
            if (min_struct_dist1 > float(target_structure_scores0[ i ])):
                min_struct_dist1 = float(target_structure_scores0[ i ])
                j1 = i

print("min struct_dist1 sequence:", Whole_rnas[ j1 ])
print("min_struct_dist1:", min_struct_dist1)


print("min struct_dist2 sequence:", Whole_rnas[ j2 ])
print("min_struct_dist2:", min_struct_dist2)



print("min struct_dist3 sequence:", Whole_rnas[ j3 ])
print("min_struct_dist3:", min_struct_dist3)


'''


j=0
min_struct_dist = 1000
with open(ft, 'w') as outfile0:
    for i in range(len(Whole_gcs_scores)):
        temp_f = Whole_forbidden_scores[ i ]
        temp_m = Whole_motif_scores[ i ]
        if ((float(Whole_gcs_scores[i])<0.02) and (len(Whole_rnas[ i ])>1) and (int(temp_f[0])==0) and (int(temp_m[0])==0)):
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
                if (min_MFE > float(Whole_MFE_scores[ i ])):
                    min_MFE = float(Whole_MFE_scores[ i ])
                    j = i

print("min MFE sequence with constraints:", Whole_rnas[ j ])
print("Whole_gcs:", Whole_gcs[ j ])
print("Whole_gcs_scores:", Whole_gcs_scores[ j ])
print("Whole_motif_scores:", Whole_motif_scores[ j ])
print("Whole_forbidden_scores:", Whole_forbidden_scores[ j ])
print("Whole_MFE_scores:", Whole_MFE_scores[ j ])
print("min_MFE:", min_MFE)
#with open(f1, 'r') as infile1, open(f2, 'r') as infile2, open(f3, 'r') as infile3, open(f4, 'r') as infile4:
#    for line in infile1:
#    i = i+1
    

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