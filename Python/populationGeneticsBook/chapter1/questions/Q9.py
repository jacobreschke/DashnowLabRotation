from PopulationGenetics import PopulationGenetics


mel_pos = 87
mel_neg = 13
n = mel_pos + mel_neg

q = PopulationGenetics.estimate_homozygous_allele_frequency(mel_neg, n)
p = 1 - q

population = PopulationGenetics(p, n)

population.run_results()
population.print_results()
print(population.ratio_hetero_to_homo_zygote)