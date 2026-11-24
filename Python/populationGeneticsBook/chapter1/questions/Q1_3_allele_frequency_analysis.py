import math


def calculate_allele_frequency(num_of_alleles, total_alleles):
    return num_of_alleles/total_alleles

def calculate_variance(allele_frequency, total_alleles):
    return (allele_frequency * (1-allele_frequency))/total_alleles


def calculate_standard_error(variance, standard_deviations):
    return standard_deviations * math.sqrt(variance)

def find_confidence_intervals(allele_frequency, variance):
    lower = allele_frequency - calculate_standard_error(variance, 1)
    upper = allele_frequency + calculate_standard_error(variance, 1)
    print(f"68% confidence interval = {lower:.2%} - {upper:.2%}")

    lower = allele_frequency - calculate_standard_error(variance, 2)
    upper = allele_frequency + calculate_standard_error(variance, 2)
    print(f"95% confidence interval = {lower:.2%} - {upper:.2%}")

    lower = allele_frequency - calculate_standard_error(variance, 3)
    upper = allele_frequency + calculate_standard_error(variance, 3)
    print(f"99.7% confidence interval = {lower:.2%} - {upper:.2%}")




def main():


    print("--------------------Q1--------------------")
    num_of_alleles = 2200
    total_alleles = 3200
    frequency_of_f =calculate_allele_frequency(num_of_alleles, total_alleles)
    variance = calculate_variance(frequency_of_f, total_alleles)
    find_confidence_intervals(frequency_of_f, variance)

    print("--------------------Q2--------------------")
    frequency_of_a = calculate_allele_frequency(2*4+32, 200)
    print(f"frequency of A allele is: {frequency_of_a: .2%}")

    print("--------------------Q3--------------------")
    print("\nE6f")
    total_allele_count = 1467
    e6f_count = 175
    freq_of_r6f = calculate_allele_frequency(e6f_count, total_allele_count)
    variance_of_r6f = calculate_variance(freq_of_r6f, total_allele_count)
    find_confidence_intervals(freq_of_r6f, variance_of_r6f)

    print("\nE6s")
    e6s_count = 314
    freq_of_e6s = calculate_allele_frequency(e6s_count, total_allele_count)
    variance_of_e6s = calculate_variance(freq_of_e6s, total_allele_count)
    find_confidence_intervals(freq_of_e6s, variance_of_e6s)

    print("\nECf")
    ecf_count = 436
    freq_of_ecf = calculate_allele_frequency(ecf_count, total_allele_count)
    variance_of_ecf = calculate_variance(freq_of_ecf, total_allele_count)
    find_confidence_intervals(freq_of_ecf, variance_of_ecf)

    print("\nECs")
    ecs_count = 53
    freq_of_ecs = calculate_allele_frequency(ecs_count, total_allele_count)
    variance_of_ecs = calculate_variance(freq_of_ecs, total_allele_count)
    find_confidence_intervals(freq_of_ecs, variance_of_ecs)

    print("\nODHf")
    odhf_count = 460
    freq_of_odhf = calculate_allele_frequency(odhf_count, total_allele_count)
    variance_of_odhf = calculate_variance(freq_of_odhf, total_allele_count)
    find_confidence_intervals(freq_of_odhf, variance_of_odhf)

    print("\nODHs")
    odhs_count = 29
    freq_of_odhs = calculate_allele_frequency(odhs_count, total_allele_count)
    variance_of_odhs = calculate_variance(freq_of_odhs, total_allele_count)
    find_confidence_intervals(freq_of_odhs, variance_of_odhs)


if __name__ == "__main__":
    main()