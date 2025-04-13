import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def load_solution(filename):
    with open(filename, 'r') as f:
        # Read all lines
        data = []
        for line in f:
            row = list(map(float, line.strip().split()))
            if row:  # skip empty lines
                data.append(row)
        
        U = np.array(data)
        M, N = U.shape[1]-1, U.shape[0]  # space points, time steps
        x = np.linspace(0, 1, M+1)
        t = np.linspace(0, 1, N)
        return x, t, U

def animate_wave(x, U):
    fig, ax = plt.subplots(figsize=(10, 6))
    line, = ax.plot(x, U[0], 'b-', linewidth=2)
    
    ax.set_xlim(0, 1)
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlabel('Position x')
    ax.set_ylabel('Displacement u(x,t)')
    ax.set_title('Wave Propagation')
    ax.grid(True)
    
    def update(frame):
        line.set_ydata(U[frame])
        ax.set_title(f'Wave at t = {frame/len(U):.2f}')
        return line,
    
    ani = FuncAnimation(fig, update, frames=len(U), interval=50, blit=True)
    plt.show()
    return ani

# Load and visualize
x, t, U = load_solution('solution.txt')

# 1. Animation
animate_wave(x, U)

# 2. Heatmap plot
plt.figure(figsize=(10, 6))
plt.imshow(U, extent=[0, 1, 0, 1], aspect='auto', cmap='RdBu')
plt.colorbar(label='Displacement')
plt.xlabel('Position x')
plt.ylabel('Time t')
plt.title('Wave Equation Solution')
plt.show()