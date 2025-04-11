import numpy as np
import matplotlib.pyplot as plt

# Load the solution from file
U = np.loadtxt("solution")

# Extract dimensions
M_plus_1, N_plus_1 = U.shape
x = np.linspace(0, 1, M_plus_1)  # Space grid
t = np.linspace(0, 1, N_plus_1)  # Time grid

# Create a heatmap
plt.figure(figsize=(10, 6))
plt.imshow(U, extent=[0, 1, 0, 1], origin="lower", aspect="auto", cmap="hot")
plt.colorbar(label="Temperature (u)")
plt.xlabel("Time (t)")
plt.ylabel("Space (x)")
plt.title("Heat Equation Solution: u(x,t)")
plt.show()