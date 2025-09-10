import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from skimage import color
from matplotlib import cm

def create_average_colormap(cmap1, cmap2, n_colors=256):
    """
    Create a perceptually uniform colormap by averaging two perceptually uniform colormaps.

    Args:
    cmap1 (str): The name of the first colormap.
    cmap2 (str): The name of the second colormap.
    n_colors (int, optional): The number of colors in the resulting colormap. Default is 256.

    Returns:
    LinearSegmentedColormap: The custom colormap.
    """

    colors1 = cmap1(np.linspace(0, 1, n_colors))[:, :3]
    colors2 = cmap2(np.linspace(0, 1, n_colors))[:, :3]

    lab_colors1 = color.rgb2lab(colors1)
    lab_colors2 = color.rgb2lab(colors2)

    averaged_lab_colors = (lab_colors1 + lab_colors2) / 2
    averaged_rgb_colors = color.lab2rgb(averaged_lab_colors)

    custom_cmap = LinearSegmentedColormap.from_list("average_colormap", averaged_rgb_colors)

    return custom_cmap


if __name__ == "__main__":
    cmap = cm.viridis
    print("Number of points in colormap:", cmap.N)  # Output: Number of points in colormap: 256


    # Example usage:
    cmap1 = cm.viridis
    cmap2 = cm.plasma

    average_cmap = create_average_colormap(cmap1, cmap2)

    # Create sample data
    data = np.random.rand(10, 10)

    # Plot the data with the custom colormap
    plt.imshow(data, cmap=average_cmap)
    plt.colorbar()
    plt.show()
