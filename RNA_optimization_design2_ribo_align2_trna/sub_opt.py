import RNA

#sequence = "UGACAUUGCCUAACGCCUAGGCGUC"
#sequence = "AAUUGCUGCCAAUAAUCUUCGCUUUGUUCAAAGCGAAUUUACUGCGGAAAGGCUCCUUUCCGCAGAAAUCUAGGAAGAUAUAUUGGAGGCAGCAAU"
sequence = "AAUUGCUGCCAAUAAUCUUCACUUUGUUCAAAGUGAAUUUACUGCGGAAAGGCUCCUUUCCGCAGAAAUCUAGGAAGAUAUAUUGGAGGCAGCAAU"

# Set global switch for unique ML decomposition
RNA.cvar.uniq_ML = 5

subopt_data = { 'counter' : 1, 'sequence' : sequence }

# Print a subopt result as FASTA record
def print_subopt_result(structure, energy, data):
    if not structure == None:
        print(">subopt {:d}".format(data['counter']))
        print("{}\n{} [{:6.2f}]".format(data['sequence'], structure, energy))
        # increase structure counter
        data['counter'] = data['counter'] + 1

# Create a 'fold_compound' for our sequence
a = RNA.fold_compound(sequence)

# Enumerate all structures 500 dacal/mol = 5 kcal/mol arround
# the MFE and print each structure using the function above
a.subopt_cb(500, print_subopt_result, subopt_data);


