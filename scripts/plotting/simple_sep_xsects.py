import numpy as np
import matplotlib.pyplot as plt

# Define the epsilon range (x-axis)
epsilon = np.linspace(0, 0.9, 100)

# Slope and intercept with uncertainties
slope = 0.61
intercept = 0.30
slope_err = 0.05     # <-- update with your actual error
intercept_err = 0.04 # <-- update with your actual error

# Calculate sigma (y-axis)
sigma = slope * epsilon + intercept

# Epsilon values for red dots and their errors
epsilon_dots = np.array([0.292, 0.779])
sigma_dots = slope * epsilon_dots + intercept
sigma_errors = np.array([0.02, 0.03])  # Error bars for each point

# Plot
plt.figure(figsize=(5, 5))
plt.plot(epsilon, sigma, label=r'$\sigma = 0.61\,\epsilon + 0.30$', color='blue', linewidth=2)
plt.errorbar(epsilon_dots, sigma_dots, yerr=sigma_errors, fmt='o', color='red',
             markersize=8, capsize=5, label='Data Points', zorder=5)
plt.xlabel(r'$\epsilon$', fontsize=14)
plt.ylabel(r'$d\sigma / dt\, d\phi$ ($\mu$b/GeV$^2$)', fontsize=14)
plt.title('Sigma vs Epsilon', fontsize=16)
plt.xlim(0, 0.9)
plt.ylim(0.1, None)

# Legend box on top left
plt.legend(loc="upper left", fontsize=12)

# Sigma_L and Sigma_T values at lower right (axes fraction)
textstr = (
    r'$\sigma_L = %.2f \pm %.2f$' % (slope, slope_err) + '\n' +
    r'$\sigma_T = %.2f \pm %.2f$' % (intercept, intercept_err)
)
plt.text(
    0.65, 0.08, textstr,
    transform=plt.gca().transAxes,
    fontsize=13, va='bottom', ha='left',
    bbox=dict(facecolor='white', alpha=0.7, edgecolor='gray')
)

plt.tight_layout()

# Save the figure as PNG
plt.savefig('sigma_vs_epsilon.png', dpi=300)
plt.close()