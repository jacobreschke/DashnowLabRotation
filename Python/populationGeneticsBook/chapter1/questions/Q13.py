from MultiAllelePopulationGenetics import MultiAllelePopulationGenetics

alleles = {
    "Adh-1": 0.11,
    "Adh-2": 0.84,
    "Adh-3": 0.01,
    "Adh-4": 0.04,
}

population = MultiAllelePopulationGenetics(alleles)

for (a, b), frequency in population.hardy_weinberg().items():
    print(f"{a}/{b}: {frequency:.4f}")

