import math
import numpy as np
from random import random as rnd
from scipy.integrate import quad
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


def ln_PDF(args):
    D, mu, sig = args
    if D == 0.0:
        return 0.0
    p = 1 / (D * sig * math.sqrt(2 * math.pi))
    _exp = -((math.log(D) - mu) ** 2) / (2 * sig ** 2)
    return p * math.exp(_exp)


def tln_PDF(args):
    D, mu, sig, F_DMin, F_DMax = args
    return ln_PDF((D, mu, sig)) / (F_DMax - F_DMin)


def F_tlnpdf(args):
    mu, sig, D_Min, D_Max, D, F_DMax, F_DMin = args
    if D > D_Max or D < D_Min:
        return 0
    P, _ = quad(lambda D: tln_PDF((D, mu, sig, F_DMin, F_DMax)), D_Min, D)
    return P


def makeSample(args, N=100):
    ln_Mean, ln_sig, D_Min, D_Max, F_DMax, F_DMin = args
    probs = [rnd() for _ in range(N)]

    def root_function(D, prob):
        return F_tlnpdf((ln_Mean, ln_sig, D_Min, D_Max, D, F_DMax, F_DMin)) - prob

    d_s = [fsolve(root_function, x0=(D_Min + D_Max) / 2, args=(p,))[0] for p in probs]
    return d_s


def sampleStats(D, doPrint=False):
    N = len(D)
    mean = sum(D) / N
    var = sum((d - mean) ** 2 for d in D) / (N - 1)
    if doPrint:
        print(f"mean = {mean:.3f}, var = {var:.3f}")
    return mean, var


def getFDMaxFDMin(args):
    mean_ln, sig_ln, D_Min, D_Max = args
    F_DMax, _ = quad(lambda D: ln_PDF((D, mean_ln, sig_ln)), 0, D_Max)
    F_DMin, _ = quad(lambda D: ln_PDF((D, mean_ln, sig_ln)), 0, D_Min)
    return F_DMin, F_DMax


def makeSamples(args):
    mean_ln, sig_ln, D_Min, D_Max, F_DMax, F_DMin, N_sampleSize, N_samples, doPrint = args
    Samples, Means = [], []
    for n in range(N_samples):
        sample = makeSample((mean_ln, sig_ln, D_Min, D_Max, F_DMax, F_DMin), N=N_sampleSize)
        Samples.append(sample)
        sample_Stats = sampleStats(sample)
        Means.append(sample_Stats[0])
        if doPrint:
            print(f"Sample {n}: mean = {sample_Stats[0]:.3f}, var = {sample_Stats[1]:.3f}")
    return Samples, Means


def main():
    mean_ln = float(input("Enter mean of ln(D) or press enter for default = 0.693: ") or 0.693)
    sig_ln = float(input("Enter standard deviation of ln(D) or press enter for default = 1.0: ") or 1.0)
    D_Min = float(input("Enter minimum sieve aperture size or press enter for default = .375: ") or .375)
    D_Max = float(input("Enter maximum sieve aperture size or press enter for default = 1.0: ") or 1.0)
    N_samples = int(input("Enter number of samples or press enter for default = 11: ") or 11)
    N_sampleSize = int(input("Enter sample size or press enter for default = 100: ") or 100)

    F_DMin, F_DMax = getFDMaxFDMin((mean_ln, sig_ln, D_Min, D_Max))
    Samples, Means = makeSamples((mean_ln, sig_ln, D_Min, D_Max, F_DMax, F_DMin, N_sampleSize, N_samples, True))
    stats_of_Means = sampleStats(Means)
    print(f"Mean of the sampling mean:  {stats_of_Means[0]:.3f}")
    print(f"Variance of the sampling mean:  {stats_of_Means[1]:.6f}")


if __name__ == '__main__':
    main()
