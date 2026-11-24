from MultiAllelePopulationGenetics import MultiAllelePopulationGenetics

# This requires a different class/calculations

observed = {
    ("A", "A"): 724,
    ("B", "B"): 110,
    ("O", "O"): 763,
    ("A", "B"): 20,
}
n = sum(observed.values())

alleles = MultiAllelePopulationGenetics.allele_frequencies(observed)

population = MultiAllelePopulationGenetics(alleles, n)

population.print_results()