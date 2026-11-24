from Python.populationGeneticsBook.chapter1.questions.PopulationGenetics import PopulationGenetics

MM = 298
MN = 489
NN = 213

n = MM + MN + NN

observed = [MM, MN, NN]

# Calculate M allele frequency
M_count = 2 * MM + MN
p = PopulationGenetics.allele_frequency(M_count, 2 * n)

# Create model
population = PopulationGenetics(p)

print(f"M frequency: {population.p:.4f}")
print(f"N frequency: {population.q:.4f}")

# Expected genotype frequencies
P, Q, R = population.hardy_weinberg()


print(f"Observed counts: \n"
      f"MM: {MM}\n"
      f"MN: {MN}\n"
      f"NN: {NN}\n")


print(f"Expected counts: \n"
      f"MM: {population.expected_counts(n)[0]}\n"
      f"MN: {population.expected_counts(n)[1]}\n"
      f"NN: {population.expected_counts(n)[2]}\n")


chi_val = population.chi_square(observed)
deg_freedom = population.degrees_of_freedom(observed)

print(f"The P value is {population.p_value(observed)}")