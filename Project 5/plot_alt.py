from common import *

h = 0.005
T = 0.008
dt = 2.5 * 10**(-5)
N_timesteps = int(T / dt)
n = int(0.001/dt)
data = pa.cx_cube()

x_idx = int(0.85 / h)   # x index where the detector is located

user_input = input("Compile animations? (y/n): ")

cm = 1/2.54
params = {
    'legend.fontsize': 10,
    'font.size': 11,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'lines.markersize': 3.0
    }

plt.rcParams.update(params)

def initalize_data(data):
    data = load(data)
    print(data.shape)
    data_conj = np.conj(data)   # complex conjugate of U_n
    p_data = np.multiply(data, data_conj) # probabilities
    return p_data, data

def run_animations(p_data, label):
    animate(np.real(p_data), label + '.mp4', h=h, T=T, dt=dt)

def error_plot(data, plot_label, anim_label):
    p_data, data = initalize_data(data)
    if user_input == "y":
        run_animations(p_data, anim_label)
    
    total_grid_prob = []
    for i in range(N_timesteps):
        total_grid_prob.append(np.real(np.sum(p_data[i,:,:])))  # sum of probabilities in grid
    
    error = np.array(total_grid_prob) - 1
    time = np.arange(0, T, step=dt)
    plt.plot(time, error, label = plot_label, linestyle='dotted')
    # plt.hlines(y=0, xmin=0, xmax=T, linestyles='dashed', colors='k')

def prob_detectorsheet(data, x_idx, save_label):
    data_x = data[-1, x_idx, :]
    p_data_x = np.real(np.multiply(data_x, np.conj(data_x)))

    p_data_x /= np.sum(p_data_x)
    plt.figure(figsize=(9*cm,7*cm))
    plt.plot(p_data_x)
    plt.xticks([0, 99, 199], ["0", "0.5", "1"])
    plt.xlabel("$y$-values"), plt.ylabel("$P(y\, | \, x = 0.8, t = 0.002)$")
    plt.tight_layout(), plt.grid()
    plt.savefig(save_label)

# data = pa.cx_mat()
# data.load("data/potential_1_slit.bin", pa.arma_binary)
# data = np.array(data, dtype=np.cdouble)
# data = np.real(data)
# colormap(data, fig_label="Potential", save=True, save_label="potenial_1_slits")


# Plot the total probability error
plt.figure(figsize=(9*cm,7*cm))
error_plot("data/problem7partA.bin", "No barrier", "anim/no_barrier")
error_plot("data/problem7partB.bin", "Double slit barrier", "anim/double_slit_barrier")
plt.xlabel("Time"), plt.ylabel("Error")
plt.legend(), plt.tight_layout(), plt.grid()
plt.savefig("figs/prob_error.pdf")

# Plot the double slit experiment at t = 0, 0.001 and 0.002
p_data, data = initalize_data("data/problem8.bin")
colormap(np.real(p_data[0]), fig_label="Probability", save_label="prob_T0.0")
colormap(np.real(p_data[n]), fig_label="Probability", save_label="prob_T0.001")
colormap(np.real(p_data[-1]), fig_label="Probability", save_label="prob_T0.002")

colormap(np.real(data[0]), fig_label="Re($U^n$)", save_label="real_U_T0.0")
colormap(np.real(data[n]), fig_label="Re($U^n$)", save_label="real_U_T0.001")
colormap(np.real(data[-1]), fig_label="Re($U^n$)", save_label="real_U_T0.002")

colormap(np.imag(data[0]), fig_label="Im($U^n$)", save_label="imag_U_T0.0")
colormap(np.imag(data[n]), fig_label="Im($U^n$)", save_label="imag_U_T0.001")
colormap(np.imag(data[-1]), fig_label="Im($U^n$)", save_label="imag_U_T0.002")

# Plot probability at the detector sheet placed at x = 0.8 and animate the output
prob_detectorsheet(data, x_idx, "figs/cond_prob_x0.8.pdf")
if user_input == "y":
    run_animations(p_data, label="anim/problem_9_twoslits")

p_data, data = initalize_data("data/problem9_3_slits.bin")
prob_detectorsheet(data, x_idx, "figs/cond_prob_x0.8_3slits.pdf")
if user_input == "y":
    run_animations(p_data, label="anim/problem_9_threeslits")

p_data, data = initalize_data("data/problem9_1_slit.bin")
prob_detectorsheet(data, x_idx, "figs/cond_prob_x0.8_1slit.pdf")
if user_input == "y":
    run_animations(p_data, label="anim/problem_9_oneslit")