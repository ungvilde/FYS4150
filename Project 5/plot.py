from common import *


data = load("data/problem7partA.bin") 
data_conj = np.conj(data) # complex conjugate of U_n
p_data = np.multiply(data, data_conj) # probabilities

total_grid_prob = []
for i in range(320):
    total_grid_prob.append(np.real(np.sum(p_data[i,:,:]))) # sum of probabilities in grid

error = (np.array(total_grid_prob) - 1) # difference between sum of probs and 1.0
T = 0.008
dt = 2.5 * 10**(-5)
time = np.arange(dt, T+dt, step=dt)
plt.plot(time, error, label = "$v_0 = 0$")

data = load("data/problem7partB.bin")
data_conj = np.conj(data)
p_data = np.multiply(data, data_conj)

total_grid_prob = []
for i in range(320):
    total_grid_prob.append(np.real(np.sum(p_data[i,:,:]))) # sum of probabilities in grid

probs_diff = (np.array(total_grid_prob) - 1)
plt.plot(time, probs_diff, label = "$v_0 = 10^{10}$")
plt.xlabel("Time")
plt.ylabel("Error")
plt.legend()
plt.tight_layout()
plt.show()

