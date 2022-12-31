import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Generate some data for the plot
x = np.linspace(-5, 5, 100)
y = np.sin(x)

# Create the plot using the `plot` function
plt.plot(x, y)

# Define a custom set of colors for the color bar
colors = ["#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"]
cmap = ListedColormap(colors)

# Add a color bar to the plot using the `colorbar` function
plt.colorbar(cmap=cmap)

# Show the plot
plt.show()
