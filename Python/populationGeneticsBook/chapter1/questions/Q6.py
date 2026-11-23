import math


def calculate_allele_frequency(num_of_alleles, total_alleles):
    return num_of_alleles/total_alleles


num_of_alleles = 29
total_alleles = 60


p = calculate_allele_frequency(num_of_alleles, total_alleles)
q = 1 - p

P = p ** 2
Q = 2 * p * q
R = q ** 2

print(f"BB: {P:.4f}")
print(f"Bb: {Q:.4f}")
print(f"bb: {R:.4f}")
print(f"Total: {P + Q + R:.4f}\n")


BB_BB = P ** 2
Bb_Bb = Q ** 2
bb_bb = R ** 2
BB_Bb = 2 * P * Q
Bb_bb = 2 * Q * R
BB_bb = 2 * P * R


print("Mating frequencies")
print(f"BBxBB | {BB_BB:.4f}")
print(f"BBxBb | {BB_Bb:.4f}")
print(f"BBxbb | {BB_bb:.4f}")
print(f"BbxBb | {Bb_Bb:.4f}")
print(f"Bbxbb | {Bb_bb:.4f}")
print(f"bbxbb | {bb_bb:.4f}")


# (P ** 2) + 2PQ/2 + Q**2/4
p_prime = BB_BB + (BB_Bb / 2) + (Bb_Bb / 4)
q_prime = (BB_Bb / 2) + BB_bb + (Bb_Bb / 2) + (Bb_bb / 2)
r_prime = (Bb_Bb / 4) + (Bb_bb / 2) + bb_bb

print(f"P': {p_prime:.4f}")
print(f"Q': {q_prime:.4f}")
print(f"R': {r_prime:.4f}")

n = 60
m = 1
k = 1
j = 4
c = 29 + ((m-k) * n)  # 29 + 0

p_nuc = (c - n * (m - k))/(j * c)
print(f"percentage: {p_nuc:.4f}")  # 1/4

sum_squared = c ** 2

h_nuc = (n * c - sum_squared) / (j * c * (n - 1))
print(f"Heterozygosity = {h_nuc:.4f}")
