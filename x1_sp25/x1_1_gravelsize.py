import math
import random

def truncated_lognormal(mu, sigma, Dmin, Dmax, N=100):
    """
    Generates a sample of N rock diameters following a truncated log-normal distribution.

    :param mu: Mean of ln(D)
    :param sigma: Standard deviation of ln(D)
    :param Dmin: Minimum allowed rock diameter
    :param Dmax: Maximum allowed rock diameter
    :param N: Number of samples to generate (default=100)
    :return: List of rock diameters
    """
    samples = []
    while len(samples) < N:
        D = math.exp(random.gauss(mu, sigma))  # Generate log-normal value
        if Dmin <= D <= Dmax:
            samples.append(D)
    return samples

def compute_statistics(samples):
    """
    Computes the mean and variance of a sample.

    :param samples: List of sample values
    :return: (mean, variance)
    """
    mean = sum(samples) / len(samples)
    variance = sum((x - mean) ** 2 for x in samples) / (len(samples) - 1)
    return mean, variance

def main():
    """
    Runs the gravel simulation, generating rock samples and computing statistics.
    """
    print("Mean of ln(D) for the pre-sieved rocks? (ln(2.0)=0.693, where D is in inches): ")
    mu = float(input(">> ") or 0.693)

    print("Standard deviation of ln(D) for the pre-sieved rocks?")
    sigma = float(input(">> ") or 1.0)

    print("Large aperture size?")
    Dmax = float(input(">> ") or 1.000)

    print("Small aperture size?")
    Dmin = float(input(">> ") or 0.375)

    print("How many samples?")
    num_samples = int(input(">> ") or 11)

    print("How many items in each sample?")
    N = int(input(">> ") or 100)

    all_means = []
    all_variances = []

    print("\n")
    for i in range(num_samples):
        sample = truncated_lognormal(mu, sigma, Dmin, Dmax, N)
        mean, variance = compute_statistics(sample)
        all_means.append(mean)
        all_variances.append(variance)
        print(f"Sample {i}: mean = {mean:.3f}, var = {variance:.3f}")

    # Compute mean and variance of sample means
    mean_of_means = sum(all_means) / num_samples
    variance_of_means = sum((x - mean_of_means) ** 2 for x in all_means) / (num_samples - 1)

    print(f"\nMean of the sampling mean: {mean_of_means:.3f}")
    print(f"Variance of the sampling mean: {variance_of_means:.6f}")

    print("\nGo again? (Y/N)")
    repeat = input(">> ").strip().lower()
    if repeat == "y":
        main()

if __name__ == "__main__":
    main()