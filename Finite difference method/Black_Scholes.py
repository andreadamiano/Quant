import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load data from file (assuming your C++ program saved it)
try:
    data = np.loadtxt("solution")
    M = 100  # Should match your C++ M value
    N = 1000  # Should match your C++ N value
    U = data.reshape((M+1, N+1))  # Reshape to original dimensions
    
    # Create grids
    S_max = 200.0
    T = 1.0
    S = np.linspace(0, S_max, M+1)
    t = np.linspace(0, T, N+1)
    T_grid, S_grid = np.meshgrid(t, S)
    
    # Plot 1: Option price surface
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(S_grid, T_grid, U, cmap='viridis')
    ax.set_xlabel('Stock Price (S)')
    ax.set_ylabel('Time to Maturity (T-t)')
    ax.set_zlabel('Option Price (V)')
    ax.set_title('Black-Scholes Option Price Surface')
    
    # Plot 2: Final option prices at t=0
    plt.figure(figsize=(10, 5))
    plt.plot(S, U[:, -1], 'b-', linewidth=2)
    plt.xlabel('Stock Price (S)')
    plt.ylabel('Option Price (V)')
    plt.title('Option Prices at t=0 (Maturity)')
    plt.grid(True)
    
    # Plot 3: Option price evolution for specific stock prices
    plt.figure(figsize=(10, 5))
    for S_val in [80, 100, 120]:
        idx = int(S_val / (S_max/M))  # Find closest index
        plt.plot(t, U[idx, :], label=f'S={S_val}')
    plt.xlabel('Time (t)')
    plt.ylabel('Option Price (V)')
    plt.title('Option Price Evolution for Different Stock Prices')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.show()

except FileNotFoundError:
    print("Error: Could not find solution file. Make sure your C++ program ran successfully.")
except Exception as e:
    print(f"Error loading data: {str(e)}")