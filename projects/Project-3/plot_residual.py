import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt("residual_log_coarse_p3_q3_hlle.txt")

# Split columns
iterations = data[:, 0]
residuals = data[:, 1]

# Plot
plt.figure()
plt.plot(iterations, residuals)

# Log scale (important for residuals)
plt.yscale('log')

plt.rcParams['font.size'] = 14  # Increase font size for better readability
plt.rcParams['font.family'] = 'Times New Roman'  # Use a serif font for a more professional look
# Labels
plt.xlabel("Iteration")
plt.ylabel("L1 Residual Norm")

# Grid (nice for log plots)
plt.grid(True, which="both", linestyle="--")

plt.tight_layout()
plt.savefig("residual_plot_coarse_p3_q3_hlle.png", dpi=300)