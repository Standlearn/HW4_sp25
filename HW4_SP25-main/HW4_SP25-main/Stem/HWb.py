import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import lognorm
from scipy.integrate import quad


def truncated_lognormal_plot():
    # Solicit user input for Log-Normal parameters
    mu = float(input("Enter the mean (mu) of the log-normal distribution (default =.693): ") or .693)  # .693 as default
    sigma = float(input("Enter the standard deviation (sigma) of the log-normal distribution (default =1.0): ") or 1.0)  # 1 as default
    d_min = float(input("Enter the minimum D (D_min) of the truncated log-normal distribution (default =.375): ") or .375)  # .375 as default
    d_max = float(input("Enter the maximum D (D_max) of the truncated log-normal distribution (default =1.0): ") or 1.0)  # 1.0

    # Define the scaling factor for the log-normal distribution
    scale = np.exp(mu)

    # Define x values for the range of the truncated log-normal distribution
    x = np.linspace(d_min, d_max, 500)

    # Compute the unnormalized PDF
    pdf = lognorm.pdf(x, sigma, scale=scale)  # log normal density

    # Compute normalization constant by integrating over the truncated range
    normalization_factor, _ = quad(lambda t: lognorm.pdf(t, sigma, scale=scale), d_min, d_max)

    # Normalize the PDF
    pdf /= normalization_factor

    # Compute the normalized CDF
    cdf = np.cumsum(pdf) * (x[1] - x[0])
    cdf /= cdf[-1]  # Ensure CDF ends at 1

    # Set the upper limit of integration for shading
    p_threshold = 0.75
    d_threshold = d_min + (d_max - d_min) * p_threshold
    p_value = np.interp(d_threshold, x, cdf)  # Interpolate from computed CDF

    # Plot the PDF with shaded region
    fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axs[0].plot(x, pdf, label='Normalized PDF of Truncated Log-Normal')
    axs[0].fill_between(x, pdf, where=(x <= d_threshold), color='gray', alpha=0.5)
    axs[0].set_ylabel("f(D)")
    axs[0].set_title("Truncated Log-Normal Distribution")
    axs[0].annotate(f'P(D<{d_threshold:.2f}|TLN({mu:.2f},{sigma:.2f}))={p_value:.2f}',
                    xy=(d_threshold, lognorm.pdf(d_threshold, sigma, scale=scale) / normalization_factor),
                    xytext=(d_min + (d_max - d_min) * 0.3, max(pdf) * 0.8),
                    arrowprops=dict(arrowstyle='->', connectionstyle="arc3"))

    # Add text annotation to the PDF graph with mathematical notation
    text_x = d_min + (d_max - d_min) * 0.1
    text_y = max(pdf) * 0.9
    axs[0].text(text_x, text_y,
                r'$f(D) = \frac{1}{D \sigma \sqrt{2\pi}} e^{-\frac{1}{2} \left( \frac{\ln(D) - \mu}{\sigma} \right)^2}$',
                fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

    # Plot the normalized CDF with probability threshold
    axs[1].plot(x, cdf, label='Normalized CDF of Truncated Log-Normal')
    axs[1].hlines(p_value, d_min, d_threshold, color='black', linewidth=1)
    axs[1].vlines(d_threshold, 0, p_value, color='black', linewidth=1)
    axs[1].plot(d_threshold, p_value, 'ro', markerfacecolor='white', markeredgecolor='red')
    axs[1].set_xlabel("x")
    axs[1].set_ylabel(r'$\Phi(x) = \int_{D_{min}}^{p} f(D) dD$')

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    truncated_lognormal_plot()
# endregion
