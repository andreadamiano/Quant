import numpy as np
import matplotlib.pyplot as plt

#read the data from CSV
data = np.loadtxt('BM.csv')
n_simulations = data.shape[0]
T = 1.0
n_timepoints = data.shape[1]
time = np.linspace(0, T, n_timepoints)

#plot data 
plt.figure(figsize=(10,8))

for i in range(n_simulations):
    plt.plot(time, data[i,:], alpha=0.6, linewidth=0.8)
    
#mean path
mean_path = np.mean(data, axis=0)
plt.plot(time, mean_path, 'k-', linewidth=2, label='Mean path')
    
#95% confidence interval
std_dev = np.std(data, axis=0)
plt.fill_between(time, mean_path - 1.96*std_dev, mean_path + 1.96*std_dev, color='gray', alpha=0.2, label='95% CI')


plt.title(f'Brownian Motion Simulations (n={n_simulations})')
plt.xlabel('Time')
plt.ylabel('W(t)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
