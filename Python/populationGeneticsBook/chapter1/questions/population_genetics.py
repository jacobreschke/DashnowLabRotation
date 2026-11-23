
import numpy as np
from scipy.stats import chi2

class PopulationGenetics:

    def __init__(self, p):
        # Allele frequencies
        self.p = p
        self.q = 1 - p

    # Calculate allele frequency from counts
    @staticmethod
    def allele_frequency(allele_count, total_alleles):
        return allele_count / total_alleles

    # Hardy-Weinberg genotype frequencies
    def hardy_weinberg(self):
        p = self.p
        q = self.q

        P = p ** 2
        Q = 2 * p * q
        R = q ** 2

        return P, Q, R

    # Expected genotype counts
    def expected_counts(self, population_size):
        P, Q, R = self.hardy_weinberg()

        return (
            P * population_size,
            Q * population_size,
            R * population_size
        )

    # Random mating frequencies
    def mating_frequencies(self):
        P, Q, R = self.hardy_weinberg()

        return {
            "AA x AA": P ** 2,
            "AA x Aa": 2 * P * Q,
            "AA x aa": 2 * P * R,
            "Aa x Aa": Q ** 2,
            "Aa x aa": 2 * Q * R,
            "aa x aa": R ** 2
        }

    # Genotype frequencies in next generation
    def next_generation(self):
        mating = self.mating_frequencies()

        P = (mating["AA x AA"]
             + mating["AA x Aa"] / 2
             + mating["Aa x Aa"] / 4)

        Q = (mating["AA x Aa"] / 2
             + mating["AA x aa"]
             + mating["Aa x Aa"] / 2
             + mating["Aa x aa"] / 2)

        R = (mating["Aa x Aa"] / 4
             + mating["Aa x aa"] / 2
             + mating["aa x aa"])

        return P, Q, R

    # Analyze restriction-site data
    @staticmethod
    def restriction_sites(data, m, j):
        # n = chromosomes sampled
        # m = total restriction sites
        # k = polymorphic restriction sites
        # j = nucleotides per restriction site

        n = data.shape[0]
        k = data.shape[1]

        # Count chromosomes cut at each variable site
        ci = np.sum(data, axis=0)

        # Count constant sites
        constant_sites = m - k

        # Total cuts
        c = np.sum(ci) + constant_sites * n

        # Proportion of polymorphic nucleotides
        p_nuc = (c - n * (m - k)) / (j * c)

        # Sum of squared counts
        sum_squared = np.sum(ci ** 2)
        sum_squared += constant_sites * n ** 2

        # Nucleotide heterozygosity
        h_nuc = (n * c - sum_squared) / (j * c * (n - 1))

        return p_nuc, h_nuc

    def chi_square(self, observed):
        expected = self.expected_counts(sum(observed))

        chi = 0

        for obs, exp in zip(observed, expected):
            chi += ((obs - exp)**2) / exp

        return chi

    def degrees_of_freedom(self, num_of_categories, estimated_parameters=1):
        return len(num_of_categories) - 1 - estimated_parameters


    def p_value(self, observed, estimated_parameters=1):
        statistic = self.chi_square(observed)
        deg_freedom = self.degrees_of_freedom(observed, estimated_parameters)

        return chi2.sf(statistic, deg_freedom)