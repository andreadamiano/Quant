import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("result.csv")
S = data[:,0]
# v = data[:,1]
n_timesteps = len(data)
time = np.linspace(0,1, n_timesteps)

#visualize
plt.figure(figsize=(10,8))
plt.plot(time, S)
plt.title("Stock price")
plt.show()

# plt.figure(figsize=(10,8))
# plt.plot(time, v)
# plt.title("Volatility")
# plt.show()
