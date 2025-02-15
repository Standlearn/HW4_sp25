import random
import math


def truncated_lognormal_sample(mu, sigma, Dmin, Dmax, N=100):
    """
    Generates N samples from a truncated log-normal distribution
    with given mu (mean of ln(D)), sigma (std dev of ln(D)),
    Dmin (minimum diameter), and Dmax (maximum diameter).
    """
    samples = []
    while len(samples) < N:
        D = math.exp(random.gauss(mu, sigma))  # Generate log-normal sample
        if Dmin <= D <= Dmax:
            samples.append(D)
    return samples


def compute_statistics(samples):
    """
    Computes the mean and variance of a given list of samples.
    """
    N = len(samples)
    mean_D = sum(samples) / N
    variance_D = sum((x - mean_D) ** 2 for x in samples) / (N - 1)
    return mean_D, variance_D


def main():
    """
    Simulates an industrial gravel sieving process based on user inputs.
    The program generates 11 samples of 100 rocks each and computes
    the sample mean and variance for each, as well as the overall mean
    and variance of the sampling distribution.
    """
    # User inputs
    mu = float(input("Enter the mean of ln(D) (default: 0.0):") or 0.0)
    sigma = float(input("Enter the standard deviation of ln(D) (default: 1.0):") or 1.0)
    Dmax = float(input("Enter the maximum rock diameter (Dmax) (default: 2.0):") or 2.0)
    Dmin = float(input("Enter the minimum rock diameter (Dmin) (default: 0.5):") or 0.5)

    num_samples = 11  # Number of sample sets
    sample_size = 100  # Rocks per sample

    sample_means = []
    sample_variances = []

    for i in range(num_samples):
        sample = truncated_lognormal_sample(mu, sigma, Dmin, Dmax, sample_size)
        mean_D, variance_D = compute_statistics(sample)
        sample_means.append(mean_D)
        sample_variances.append(variance_D)
        print(f"Sample {i + 1}: Mean = {mean_D:.4f}, Variance = {variance_D:.4f}")

    # Compute statistics of the sampling distribution
    overall_mean = sum(sample_means) / num_samples
    overall_variance = sum((x - overall_mean) ** 2 for x in sample_means) / (num_samples - 1)

    print(f"\nOverall Mean of Sample Means: {overall_mean:.4f}")
    print(f"Variance of the Sample Means: {overall_variance:.4f}")


if __name__ == "__main__":
    main()
