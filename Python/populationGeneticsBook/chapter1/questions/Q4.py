import numpy as np

# 1 = +
# 0 = -
# [0]BamH1 [1]Hindlll [2]Pstl [3]Xhol [4]Pstl [5]EcoRl [6]EcoRl

data = np.array([
    [1, 0, 0, 1, 1, 0, 0],
    [1, 0, 0, 0, 0, 1, 1],
    [0, 0, 1, 0, 0, 1, 0],
    [0, 1, 0, 1, 0, 1, 1],
    [0, 1, 0, 1, 0, 1, 1],
    [0, 1, 0, 1, 0, 1, 1],
    [0, 1, 0, 1, 0, 1, 1],
    [0, 1, 0, 1, 0, 1, 1],
    [0, 1, 0, 1, 0, 1, 1],
    [0, 0, 0, 1, 0, 1, 0]
])

# n = num of chromosomes
# m = total restriction sites examined
# k = Polymorphic restriction sites
# j = nucleotides per restriction site
# ci = num of chromosomes cut at site i
# c = total cuts across all restriction sites

n = data.shape[0]
m = 23
k = data.shape[1]
j = 6


# count how many chromosomes contain each restriction site
ci = np.sum(data, axis=0) # axis=0: sums down each column
print("Chromosomes cut at each site: ", ci)

constant_sites = m - k
c = np.sum(ci) + (constant_sites * n)
print("Total cuts: ", c)

# calculate the proportion of polymorphic nucleotides
p_nuc = (c - n * (m - k))/(j * c) * 100
print(f"percentage: {p_nuc:.4f}%")

# square the counts for polymorphic sites
sum_squared = np.sum(ci**2)

# Include the 16 constant sites
sum_squared += constant_sites * (n**2)

print("Sum of squared counts:", sum_squared)

heterozygosity = (n * c - sum_squared) / (j * c * (n - 1))
print(f"Heterozygosity = {heterozygosity:.4f}")

