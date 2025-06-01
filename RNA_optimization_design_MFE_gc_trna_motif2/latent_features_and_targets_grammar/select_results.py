import os

# Directory containing the split files
directory = "./"  # Replace with your directory path if different





f1 = 'g_c.txt'
f2 = 'motif_scores1.txt'
f3 = 'MFE_scores_t.txt'
f4 = 'trna200maxlen_102k_mod.txt'
f5 = 'g_c_scores_t.txt'
f6 = 'target_structure_scores.txt'
f7 = 'motif_scores2.txt'
f8 = 'motif_scores3.txt'

ft = "selected_results.txt"
i = 0

with open(f1) as f:
    Whole_gcs = f.readlines()
for i in range(len(Whole_gcs)):
    Whole_gcs[ i ] = Whole_gcs[ i ].strip()

with open(f2) as f:
    motif_scores1 = f.readlines()
for i in range(len(motif_scores1)):
    motif_scores1[ i ] = motif_scores1[ i ].strip()

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

with open(f6) as f:
    target_structure_scores = f.readlines()
for i in range(len(target_structure_scores)):
    target_structure_scores[ i ] = target_structure_scores[ i ].strip()


with open(f7) as f:
    motif_scores2 = f.readlines()
for i in range(len(motif_scores2)):
    motif_scores2[ i ] = motif_scores2[ i ].strip()


with open(f8) as f:
    motif_scores3 = f.readlines()
for i in range(len(motif_scores3)):
    motif_scores3[ i ] = motif_scores3[ i ].strip()



j = 0
m = 0
min_MFE = 10
min_h = 100
with open(ft, 'w') as outfile0:
    for i in range(len(Whole_gcs_scores)):
        m1 = motif_scores1[ i ]
        m2 = motif_scores2[ i ]
        m3 = motif_scores3[ i ]
        if ((float(Whole_gcs_scores[i])<0.02) and (len(Whole_rnas[ i ])>1) and (int(m2[ 0 ])==0) and (int(m1[ 0 ])==0) and (int(m3[ 0 ])==0)):
                outfile0.write("RNA sequence: ")
                outfile0.write(Whole_rnas[ i ])
                outfile0.write("\n")
                outfile0.write("Whole_gcs[ i ]: ")
                outfile0.write(Whole_gcs[ i ])
                outfile0.write("\n")
                outfile0.write("motif_scores1[ i ]: ")
                outfile0.write(motif_scores1[ i ])
                outfile0.write("\n")
                outfile0.write("motif_scores2[ i ]: ")
                outfile0.write(motif_scores2[ i ])
                outfile0.write("\n")
                outfile0.write("motif_scores3[ i ]: ")
                outfile0.write(motif_scores3[ i ])
                outfile0.write("\n")
                outfile0.write("target_structure_scores[ i ]: ")
                outfile0.write(target_structure_scores[ i ])
                outfile0.write("\n")
                outfile0.write("Whole_MFE_scores[ i ]: ")
                outfile0.write(Whole_MFE_scores[ i ])
                outfile0.write("\n")
                if (min_h > float(target_structure_scores[ i ])):
                    min_h = float(target_structure_scores[ i ])
                    m = i
                if (min_h == float(target_structure_scores[ i ])):
                    print("yes")
                    if (min_MFE > float(Whole_MFE_scores[ i ])):
                        min_h = float(target_structure_scores[ i ])
                        m = i
                if (min_MFE > float(Whole_MFE_scores[ i ])):
                    min_MFE = float(Whole_MFE_scores[ i ])
                    j = i

print("min hamming distance sequence with constraints:", min_h)
print("Whole_gcs:", Whole_gcs[ m ])
print("motif_scores1:", motif_scores1[ m ])
print("motif_scores2:", motif_scores2[ m ])
print("motif_scores3:", motif_scores3[ m ])
print("Whole_h_dist:", target_structure_scores[ m ])
print("Whole_MFE_scores:", Whole_MFE_scores[ m ])
print("min_MFE:", min_MFE)





print(f"All files combined")