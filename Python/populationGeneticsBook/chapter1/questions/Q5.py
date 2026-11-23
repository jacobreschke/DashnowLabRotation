
p = 0.3579
q = 0.6421

P = p ** 2  # AA
Q = 2 * p * q  # Aa
R = q ** 2  # aa


print(f"FF: {P:.4f}")
print(f"FS: {Q:.4f}")
print(f"SS: {R:.4f}")
print(f"Total: {P + Q + R:.4f}")


FF_FF = P ** 2
FS_FS = Q ** 2
SS_SS = R ** 2
FF_FS = 2 * P * Q
FS_SS = 2 * Q * R
FF_SS = 2 * P * R

print("Mating frequencies")
print(f"FFxFF | {FF_FF:.4f}")
print(f"FFxFS | {FF_FS:.4f}")
print(f"FFxSS | {FF_SS:.4f}")
print(f"FSxFS | {FS_FS:.4f}")
print(f"FSxSS | {FS_SS:.4f}")
print(f"SSxSS | {SS_SS:.4f}")


# (P ** 2) + 2PQ/2 + Q**2/4
p_prime = FF_FF + (FF_FS / 2) + (FS_FS / 4)
q_prime = (FF_FS / 2) + FF_SS + (FS_FS / 2) + (FS_SS / 2)
r_prime = (FS_FS / 4) + (FS_SS / 2) + SS_SS

print(f"P': {p_prime:.4f}")
print(f"Q': {q_prime:.4f}")
print(f"R': {r_prime:.4f}")
