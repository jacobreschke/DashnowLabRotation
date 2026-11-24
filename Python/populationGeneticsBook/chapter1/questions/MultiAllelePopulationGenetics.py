import numpy as np
from scipy.stats import chi2
from itertools import combinations_with_replacement
from collections import Counter


class MultiAllelePopulationGenetics:

    def __init__(self, alleles, n=None):
        if not np.isclose(sum(alleles.values()), 1):
            raise ValueError("Allele frequencies must sum to 1.")

        if any(p < 0 or p > 1 for p in alleles.values()):
            raise ValueError("Invalid allele frequency.")

        if n is not None and n <= 0:
            raise ValueError("Population size must be positive.")

        self.alleles = alleles.copy()
        self.n = n
        self.genotypes = {}

    # Calculate allele frequencies from observed genotype counts
    @staticmethod
    def allele_frequencies(observed):
        counts = Counter()

        for (a, b), count in observed.items():
            counts[a] += count
            counts[b] += count

        total = sum(counts.values())

        return {
            allele: count / total
            for allele, count in counts.items()
        }

    # Estimate allele frequency from homozygous genotype counts
    # Assumes Hardy-Weinberg equilibrium
    @staticmethod
    def estimate_homozygous_allele_frequency(
            homozygous_count, population_size):

        R = homozygous_count / population_size
        return np.sqrt(R)

    # Hardy-Weinberg genotype frequencies
    def hardy_weinberg(self):
        self.genotypes = {}

        for a, b in combinations_with_replacement(self.alleles, 2):

            p = self.alleles[a]
            q = self.alleles[b]

            if a == b:
                frequency = p ** 2
            else:
                frequency = 2 * p * q

            self.genotypes[(a, b)] = frequency

        return self.genotypes

    # Expected genotype counts
    def expected_counts(self):
        if self.n is None:
            raise ValueError("Sample size is required.")

        genotypes = self.hardy_weinberg()

        return {
            genotype: frequency * self.n
            for genotype, frequency in genotypes.items()
        }

    # Number of possible diploid genotypes
    def genotype_count(self):
        k = len(self.alleles)

        return k * (k + 1) // 2

    # Random mating frequencies
    def mating_frequencies(self):
        genotypes = self.hardy_weinberg()
        mating = {}

        for g1, g2 in combinations_with_replacement(genotypes, 2):

            f1 = genotypes[g1]
            f2 = genotypes[g2]

            if g1 == g2:
                frequency = f1 ** 2
            else:
                frequency = 2 * f1 * f2

            mating[(g1, g2)] = frequency

        return mating

    # Genotype frequencies in the next generation
    def next_generation(self):
        mating = self.mating_frequencies()

        offspring = {
            genotype: 0.0
            for genotype in self.genotypes
        }

        # Preserve the allele ordering used by hardy_weinberg()
        allele_order = {
            allele: i
            for i, allele in enumerate(self.alleles)
        }

        for (parent1, parent2), frequency in mating.items():

            for a in parent1:
                for b in parent2:

                    genotype = tuple(
                        sorted((a, b), key=allele_order.get)
                    )

                    offspring[genotype] += frequency / 4

        return offspring

    # Analyze restriction-site data
    @staticmethod
    def restriction_sites(data, m, j):

        n = data.shape[0]
        k = data.shape[1]

        ci = np.sum(data, axis=0)
        constant_sites = m - k
        c = np.sum(ci) + constant_sites * n
        p_nuc = (c - n * (m - k)) / (j * c)

        sum_squared = np.sum(ci ** 2)
        sum_squared += constant_sites * n ** 2
        h_nuc = (
            (n * c - sum_squared)
            / (j * c * (n - 1))
        )

        return p_nuc, h_nuc

    # Expected homozygosity
    def homozygosity(self):
        return sum(
            p ** 2
            for p in self.alleles.values()
        )

    # Expected heterozygosity
    def heterozygosity(self):
        return 1 - self.homozygosity()

    # Observed heterozygosity
    @staticmethod
    def observed_heterozygosity(observed):
        heterozygous = sum(
            count
            for (a, b), count in observed.items()
            if a != b
        )

        return heterozygous / sum(observed.values())

    # Observed homozygosity
    @staticmethod
    def observed_homozygosity(observed):
        return (
            1 -
            MultiAllelePopulationGenetics.observed_heterozygosity(
                observed
            )
        )

    # Ratio of heterozygous to homozygous individuals
    def heterozygous_to_homozygous_ratio(self):
        return self.heterozygosity() / self.homozygosity()

    # Chi-square statistic
    def chi_square(self, observed):
        expected = self.expected_counts()

        if sum(observed.values()) != self.n:
            raise ValueError(
                "Observed counts must sum to population size."
            )

        # Normalize genotype ordering
        allele_order = {
            allele: i
            for i, allele in enumerate(self.alleles)
        }

        normalized = Counter()

        for genotype, count in observed.items():

            if any(a not in self.alleles for a in genotype):
                raise ValueError("Unknown allele.")

            genotype = tuple(
                sorted(genotype, key=allele_order.get)
            )

            normalized[genotype] += count

        chi = 0

        for genotype, exp in expected.items():

            obs = normalized.get(genotype, 0)

            if exp == 0:
                if obs > 0:
                    return np.inf
                continue

            chi += ((obs - exp) ** 2) / exp

        return chi

    # Degrees of freedom
    def degrees_of_freedom(self, estimated_parameters=None):

        k = sum(
            p > 0
            for p in self.alleles.values()
        )

        categories = k * (k + 1) // 2

        if estimated_parameters is None:
            estimated_parameters = k - 1

        return categories - 1 - estimated_parameters

    # Chi-square p-value
    def p_value(self, observed, estimated_parameters=None):

        statistic = self.chi_square(observed)

        df = self.degrees_of_freedom(estimated_parameters)

        if df <= 0:
            raise ValueError("Insufficient degrees of freedom.")

        return chi2.sf(statistic, df)

    # Sampling variance from directly observed allele counts
    # Assumes independent allele copies under HWE
    def sampling_variance(self, allele):
        if self.n is None:
            raise ValueError("Sample size is required.")

        p = self.alleles[allele]

        return p * (1 - p) / (2 * self.n)

    # Standard error of directly measured allele frequency
    def standard_error(self, allele):
        return np.sqrt(self.sampling_variance(allele))

    # Sampling variance when estimated from homozygous counts
    def homozygous_sampling_variance(self, allele):
        if self.n is None:
            raise ValueError("Sample size is required.")

        R = self.alleles[allele] ** 2

        return (1 - R) / (4 * self.n)

    # Standard error when estimated from homozygous counts
    def homozygous_standard_error(self, allele):
        return np.sqrt(
            self.homozygous_sampling_variance(allele)
        )

    # Run population analysis
    def run_results(self):
        results = {
            "allele_frequencies": self.alleles.copy(),
            "genotype_frequencies": self.hardy_weinberg(),
            "homozygosity": self.homozygosity(),
            "heterozygosity": self.heterozygosity(),
            "heterozygous_to_homozygous_ratio":
                self.heterozygous_to_homozygous_ratio(),
            "degrees_of_freedom": self.degrees_of_freedom()
        }

        if self.n is not None:
            results["expected_counts"] = self.expected_counts()

        return results

    # Print population analysis
    def print_results(self):

        if self.n is not None:
            print(f"Sample size: {self.n}")

        print(f"Number of alleles: {len(self.alleles)}")
        print(f"Possible genotypes: {self.genotype_count()}")

        print("\nAllele frequencies:")

        for allele, frequency in self.alleles.items():
            print(f"{allele}: {frequency:.4f}")

        print("\nGenotype frequencies:")

        for (a, b), frequency in self.hardy_weinberg().items():
            print(f"{a}/{b}: {frequency:.4f}")

        print(f"\nHomozygosity: {self.homozygosity():.2%}")
        print(f"Heterozygosity: {self.heterozygosity():.2%}")

        print(
            f"Heterozygous/Homozygous ratio: "
            f"{self.heterozygous_to_homozygous_ratio():.4f}"
        )

        print(f"Degrees of freedom: {self.degrees_of_freedom()}")

        if self.n is not None:

            print("\nExpected genotype counts:")

            for (a, b), count in self.expected_counts().items():
                print(f"{a}/{b}: {count:.2f}")

            print("\nAllele statistics:")

            for allele in self.alleles:
                print(f"\nAllele {allele}")
                print(f"Variance: {self.sampling_variance(allele):.6f}")
                print(f"Standard error: {self.standard_error(allele):.6f}")