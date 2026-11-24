from PopulationGenetics import PopulationGenetics


rh_positive = 230
rh_negative = 170
n = rh_positive + rh_negative

q = PopulationGenetics.estimate_homozygous_allele_frequency(rh_negative, n)
p = 1 - q

population = PopulationGenetics(p, n)
population.run_results()
population.print_results()


print(population.get_heterozygous_to_recessive_ratio())
print(population.get_heterozygous_to_dominant_ratio())

