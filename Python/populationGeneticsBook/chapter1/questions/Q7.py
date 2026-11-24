from PopulationGenetics import PopulationGenetics


SS = 99
Ss = 418
ss = 483

total_alleles = SS + Ss + ss
S_count = 2 * SS + Ss
observed = [SS, Ss, ss]

p = PopulationGenetics.allele_frequency(S_count, 2 * total_alleles)
population = PopulationGenetics(p)
P, Q, R = population.hardy_weinberg()

print(f"Frequency of S: {population.p:.4f}")
print(f"Frequency of s: {population.q:.4f}")

expected = population.expected_counts(total_alleles)
print(f"Observed counts: {observed}")
print(f"Expected counts: {expected}")

print(f"X^2: {population.chi_square(observed)})")
print(f"P value: {population.p_value(observed)}")