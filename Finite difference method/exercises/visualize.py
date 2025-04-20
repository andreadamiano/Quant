import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read the data
data = np.loadtxt("solution.csv")

# Create grid
N, M = data.shape
x = np.linspace(0, 1, N)
y = np.linspace(0, 1, M)
X, Y = np.meshgrid(y, x)  # Note: swapped for correct orientation

# 2D Heatmap plot
plt.figure(figsize=(10, 8))
plt.imshow(data, extent=[0, 10, 0, 1], origin='lower', cmap='viridis')
# plt.colorbar(label='Temperature')
plt.title("Solution to Laplace Equation (2D Heatmap)")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()

# 3D Surface plot
# fig = plt.figure(figsize=(12, 8))
# ax = fig.add_subplot(111, projection='3d')
# surf = ax.plot_surface(X, Y, data, cmap='viridis')
# fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
# ax.set_title("Solution to Laplace Equation (3D Surface)")
# ax.set_xlabel("X")
# ax.set_ylabel("Y")
# ax.set_zlabel("Temperature")
# plt.show()