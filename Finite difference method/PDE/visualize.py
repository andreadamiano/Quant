# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# # Read the data
# data = np.loadtxt("solution.csv")

# # Create grid
# N, M = data.shape
# x = np.linspace(0, 1, N)
# y = np.linspace(0, 1, M)
# X, Y = np.meshgrid(y, x)  # Note: swapped for correct orientation

# # 2D Heatmap plot
# plt.figure(figsize=(10, 8))
# plt.imshow(data, extent=[0, 10, 0, 1], origin='lower', cmap='viridis')
# # plt.colorbar(label='Temperature')
# plt.title("Solution to Laplace Equation (2D Heatmap)")
# plt.xlabel("X")
# plt.ylabel("Y")
# plt.show()

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


# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

# #load data
# U = np.loadtxt('solution.csv')  

# #plot setting
# fig, ax = plt.subplots()
# line, = ax.plot(U[0], lw=2) #plot the intial state (line, unpack the first element of the list)
# ax.set_ylim(-110, 110)
# ax.set_title("1D Wave Propagation")
# ax.set_xlabel("Space")
# ax.set_ylabel("Displacement")

# #animation function 
# def animate(i):
#     line.set_ydata(U[i]) #set new y data 
#     ax.set_text(f"Time Step {i}") #update title 
#     return line, title

# ani = animation.FuncAnimation(fig, animate, 
# frames=range(0, U.shape[0], 1), #n of frames 
# interval=50, #deay between frames 
# blit=True) #update only the return value of the animate function 

# plt.show()


# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

# #parameters
# Nx, Ny = 200+1, 200+1 #matrix dimensions
# matrices = []

# #load data
# data = np.loadtxt('solution.csv', delimiter = ",")

# for i in range (0, len(data), Ny):
#     matrix = data[i:i+Ny]
#     matrices.append(matrix)

# # print(matrices[1])


# #set up the figure and axis
# fig, ax = plt.subplots()
# img = ax.imshow(matrices[0], cmap='viridis', origin='lower', vmin=-1, vmax=1)
# fig.colorbar(img, ax=ax)
# ax.set_title("2D Wave Propagation")

# #animation function
# def animate(i):
#     img.set_data(matrices[i])
#     ax.set_title(f"2D Wave Propagation - Time step {i}")
#     return [img]

# #animate
# # ani = animation.FuncAnimation(fig, animate, frames=range(0, len(matrices), 1), interval=50, blit=True)
# ani = animation.FuncAnimation(fig, animate, frames=range(0, len(matrices), 1), interval=50)
# plt.show()

# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
# from mpl_toolkits.mplot3d import Axes3D

# # Parameters
# Nx, Ny = 200 + 1, 200 + 1  # Spatial resolution
# matrices = []

# # Load data
# data = np.loadtxt('solution.csv', delimiter=",")
# for i in range(0, len(data), Ny):
#     matrices.append(data[i:i + Ny])

# # Create meshgrid
# x = np.linspace(0, 1, Nx)  # Spatial domain [0, 1]
# y = np.linspace(0, 1, Ny)
# X, Y = np.meshgrid(x, y)

# # Set up the 3D plot
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(111, projection='3d')

# # Initialize surface plot with fixed Z-limits and symmetric colors
# Z = matrices[0]
# surf = ax.plot_surface(X, Y, Z, cmap='viridis', vmin=-1, vmax=1, rstride=2, cstride=2)

# # Add a reference plane at Z=0 (optional)
# ax.plot_surface(X, Y, np.zeros_like(Z), color='gray', alpha=0.3)

# # Configure axes
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# ax.set_zlim(-1.2, 1.2)  # Fixed limits to show oscillation
# ax.set_title('Standing Wave Oscillation (Z-axis)')

# # Animation function
# def animate(i):
#     ax.clear()
#     Z = matrices[i]
    
#     # Replot surface with consistent colors
#     surf = ax.plot_surface(X, Y, Z, cmap='viridis', vmin=-1, vmax=1, rstride=2, cstride=2)
#     ax.plot_surface(X, Y, np.zeros_like(Z), color='gray', alpha=0.3)  # Reference plane
    
#     # Keep axes consistent
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
#     ax.set_zlim(-1.2, 1.2)
#     ax.set_title(f'Time Step {i}')
    
#     # Adjust view for better perspective
#     # ax.view_init(elev=30, azim=45)
    
#     return [surf]

# # Animate
# ani = animation.FuncAnimation(fig, animate, frames=range(0, len(matrices), 2), interval=50, blit=False)
# plt.show()




import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

#parameters
Nx, Ny = 200 + 1, 200 + 1  #space grid 

#load data 
data = np.loadtxt('solution.csv', delimiter=",")
matrices = [data[i:i+Ny] for i in range(0, len(data), Ny)] #list comprehension

#meshgrid
x = np.linspace(0, 1, Nx)
y = np.linspace(0, 1, Ny)
X, Y = np.meshgrid(x, y)

#set up the figure
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d') #111 means a single subplot in a 1-row, 1-column grid at the first position , projection='3d' creates a 3D axis for 3D plotting.

#initial plot
Z = matrices[0] 
surf = ax.plot_surface(
    X, Y, Z, #plot a surface using x, y, z coordinates 
    cmap='viridis', #color map
    vmin=-1, vmax=1, #color limits for the surface 
    rstride=4, cstride=4, #control how may data points are used to plot the surface 
    antialiased=False #anti-aliasing is a technique used to smooth the edges of objects
)
ax.view_init(
    elev=10, #vertical angle 
    azim=20) #horizonatal angle 

# Reference plane
# ax.plot_surface(X, Y, np.zeros_like(Z), color='gray', alpha=0.2)
# ax.set(xlabel='X', ylabel='Y', zlabel='Z', zlim=(-1.2, 1.2))

#animation function
def update(i):
    Z = matrices[i] #get new matrix 
    ax.clear()

    #replot surface 
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', vmin=-1, vmax=1, rstride=4, cstride=4, antialiased=False)
    # ax.plot_surface(X, Y, np.zeros_like(Z), color='gray', alpha=0.2)

    #set view and labels 
    ax.view_init(elev=10, azim=20)
    ax.set(xlabel='X', ylabel='Y', zlabel='Z', zlim=(-1.2, 1.2), title=f'Time Step {i}') #update axis labels 
    return [surf]

#animate
ani = animation.FuncAnimation(
    fig, #figure object to animate
    update, #update function
    frames=range(0, len(matrices), 5), #which frames to animate
    interval=50, #time interval between frames 
    blit=False, #cannot be used for 3d plots 
    repeat=False) #the animation wil repeat after the last frame 

#visualize
plt.tight_layout(rect=[0, 0, 1, 0.9]) 
plt.show()

# #save animation 
# ani.save('wave_animation.gif', writer='pillow', fps=15)

