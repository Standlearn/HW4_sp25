import math
import random
from numericalMethods import Probability, GPDF


def truncated_lognormal(mu, sigma, Dmin, Dmax, N=100):
    """
    Generates a sample of N rock diameters following a truncated log-normal distribution.
    """
    samples = []
    while len(samples) < N:
        D = math.exp(random.gauss(mu, sigma))
        if Dmin <= D <= Dmax:
            samples.append(D)
    return samples


def compute_statistics(samples):
    """
    Computes the mean and variance of a sample.
    """
    mean = sum(samples) / len(samples)
    variance = sum((x - mean) ** 2 for x in samples) / (len(samples) - 1)
    return mean, variance


def t_test(sample_A, sample_B, alpha=0.05):
    """
    Performs a one-sided t-test comparing two suppliers' gravel.
    """
    mean_A, var_A = compute_statistics(sample_A)
    mean_B, var_B = compute_statistics(sample_B)

    n_A, n_B = len(sample_A), len(sample_B)
    pooled_std = math.sqrt((var_A / n_A) + (var_B / n_B))
    t_stat = (mean_A - mean_B) / pooled_std

    # Compute critical value
    c = Probability(GPDF, (0, 1), alpha, False)

    conclusion = (
        "Reject null hypothesis: samplingMeanA=samplingMeanB\n"
        "Accept alternative hypothesis: samplingMeanA>samplingMeanB"
        if t_stat > c
        else "Fail to reject null hypothesis: No significant difference."
    )

    return t_stat, c, conclusion


def main():
    """
    Runs the t-test comparison between two suppliers.
    """
    print("Enter parameters for the feedstock that is common to both company A and company B.")
    mu = float(input("Mean of ln(D) for the pre-sieved rocks? (ln(2.0)=0.693, where D is in inches): ") or 0.693)
    sigma = float(input("Standard deviation of ln(D) for the pre-sieved rocks? ") or 1.0)

    print("\nEnter parameters for company A sieving and sampling operations:")
    Dmax_A = float(input("Large aperture size? ") or 1.000)
    Dmin_A = float(input("Small aperture size? ") or 0.375)

    num_samples = int(input("How many samples? ") or 11)
    N = int(input("How many items in each sample? ") or 100)

    print("\nEnter parameters for company B sieving and sampling operations:")
    Dmax_B = float(input("Large aperture size? ") or 0.875)
    Dmin_B = float(input("Small aperture size? ") or 0.375)

    # Generate samples
    sample_A = truncated_lognormal(mu, sigma, Dmin_A, Dmax_A, N)
    sample_B = truncated_lognormal(mu, sigma, Dmin_B, Dmax_B, N)

    # Perform t-test
    t_stat, c, conclusion = t_test(sample_A, sample_B)

    print(f"\nc={c:.3f}, t_0={t_stat:.3f}")
    print(conclusion)

    print("\nGo again? (Y/N)")
    repeat = input(">> ").strip().lower()
    if repeat == "y":
        main()


if __name__ == "__main__":
    main()
