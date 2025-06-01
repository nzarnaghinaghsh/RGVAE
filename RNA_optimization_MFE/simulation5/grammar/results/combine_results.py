import os

# Directory containing the split files
directory = "./"  # Replace with your directory path if different

'''
Whole_rnas0.txt
Whole_target_struct_dist10.txt
Whole_target_struct_dist20.txt
'''
f = []
for i in range(0,15):
    f.append("Whole_MFE_scores" + str(i) + ".txt")

ft = "Whole_MFE_scores_t.txt"

with open(ft, 'w') as outfile0:
    for i in range(0,15):
        with open(f[i], 'r') as infile1:
            for line in infile1:
                outfile0.write(line)
        infile1.close()
outfile0.close()


'''
f = []
for i in range(0,50):
    f.append("Whole_target_struct_dist2" + str(i) + ".txt")

ft = "Whole_target_struct_dist2_t.txt"


with open(ft, 'w') as outfile0:
    for i in range(0,50):
        with open(f[i], 'r') as infile1:
            for line in infile1:
                outfile0.write(line)
        infile1.close()
outfile0.close()



f = []
for i in range(0,50):
    f.append("Whole_target_struct_dist3" + str(i) + ".txt")

ft = "Whole_target_struct_dist3_t.txt"


with open(ft, 'w') as outfile0:
    for i in range(0,50):
        with open(f[i], 'r') as infile1:
            for line in infile1:
                outfile0.write(line)
        infile1.close()
outfile0.close()



f = []
for i in range(0,50):
    f.append("Whole_rnas" + str(i) + ".txt")

ft = "Whole_rnas_t.txt"


with open(ft, 'w') as outfile0:
    for i in range(0,50):
        with open(f[i], 'r') as infile1:
            for line in infile1:
                outfile0.write(line)
        infile1.close()
outfile0.close()





f1 = "Whole_target_struct_dist1_t.txt"
f2 = "Whole_target_struct_dist2_t.txt"
f3 = "Whole_target_struct_dist3_t.txt"


with open(f1) as f:
    RNA_data1 = f.readlines()

for i in range(len(RNA_data1)):
    RNA_data1[ i ] = RNA_data1[ i ].strip()

with open(f2) as f:
    RNA_data2 = f.readlines()

RNA_data3 = []
for i in range(len(RNA_data2)):
    RNA_data2[ i ] = RNA_data2[ i ].strip()
    RNA_data3.append(float(RNA_data2[i]) + float(RNA_data1[i]))



with open(f3, 'w') as outfile0:
    for i in range(len(RNA_data1)):
        outfile0.write(str(RNA_data3[i]))
        outfile0.write('\n')
outfile0.close
'''

print(f"All files are combined.")