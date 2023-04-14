import numpy as np
import matplotlib.pyplot as plt

# Generate example data
x = np.linspace(0.01, 1, 1000)
y = np.exp(-5 * x) + 0.02 * np.random.rand(1000)

fig, ax1 = plt.subplots()

# Plot with logarithmic scale
ax1.set_title("Logarithmic (blue) vs. Square Root (orange) Scale")
ax1.set_xscale("log")

ax1.plot(x, y, '.', label="Logarithmic")
ax1.set_xlabel("X-axis")
ax1.set_ylabel("Y-axis (logarithmic scale)")

# Create a second y-axis with square root scale
ax2 = ax1.twinx()
ax1.set_yscale("log")
ax2.set_yscale('function', functions=(np.sqrt, np.square))
ax2.plot(x, y, '.', color='orange', label="Square Root")
ax2.set_ylabel("Y-axis (square root scale)")

# Combine legends
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.tight_layout()
plt.show()
