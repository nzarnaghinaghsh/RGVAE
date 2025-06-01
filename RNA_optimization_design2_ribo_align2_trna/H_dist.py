target_structure1 = "..((((((............((((.((((........((....)).((((......)))))))).))))..(((((.....)))))...)))))).."
target_structure2 = "..((((((....(((((...(((.....)))......)))))(((.((((......)))))))..(((...(((((.....)))))))))))))).."
    

len_t1 = len(target_structure1)
len_t2 = len(target_structure2)


subopt = []
with open("subopt_align.txt", "r") as file:
    for line in file:
        subopt.append(line.strip())




H1_dist_t = []
H2_dist_t = []


j1 = 0
j2 = 0
min_H1 = 100
min_H2 = 100
for j in range (0,6274):
    structure0 = subopt [3*j+2]
    structure = structure0[0:96]
    print("structure",structure)
    print("j",j)
    len_s = len(structure)
    min_len = min(len_t1,len_s) # Min length
    max_len = max(len_t1,len_s) # Max length
    i1 = 0
    i2 = 0
    for i in range(len_t1-len_s):
        temp1 = target_structure1[i:i+len_s]
        temp2 = target_structure2[i:i+len_s]
        print("i",i)
        print("temp1",temp1)
        H1_dist = 0 # Hamming distance
        H2_dist = 0 # Hamming distance
        for t in range(min_len):
            if (temp1[t] != structure[t]):
                H1_dist = H1_dist + 1
                #H1_dist = H1_dist + (max_len - min_len)
            
            if (temp2[t] != structure[t]):
                H2_dist = H2_dist + 1
                #H2_dist = H2_dist + (max_len - min_len)
        
        if (H1_dist < min_H1):
            min_H1 = H1_dist
            j1 = j
            i1 = i
        H1_dist_t.append(H1_dist)    
            
        if (H2_dist < min_H2):
            min_H2 = H2_dist
            j2 = j
            i2 = i
        H2_dist_t.append(H2_dist) 
                

print("structure number with min H-dist to the 1st sequence:", j1)
print("i1 =",i1)
print("min_H1 =",min_H1)  
print("H1_dist_t=",H1_dist_t)

print("structure number with min H-dist to the 2nd sequence:", j2)
print("i2 =",i2)
print("min_H2 =",min_H2) 
print("H2_dist_t=",H2_dist_t)

