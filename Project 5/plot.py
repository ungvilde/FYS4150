from common import *

h = 0.005
T = 0.008
dt = 2.5 * 10**(-5)
N_timesteps = int(T / dt)
n = int(0.001/dt)

x_idx = int(0.8 / h)   # x index where the detector is located

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


# Plot the total probability error
plt.figure(figsize=(9*cm,7*cm))
error_plot("problem7partA.bin", "No barrier", "no_barrier")
error_plot("problem7partB.bin", "Double slit barrier", "double_slit_barrier")
plt.xlabel("Time"), plt.ylabel("Error")
plt.legend(), plt.tight_layout(), plt.grid()
plt.savefig("prob_error.pdf")


# Plot the double slit experiment at t = 0, 0.001 and 0.002
p_data, data = initalize_data("problem8.bin")
colormap(np.real(p_data[0]), fig_label="Probability", save_label="prob_T0.0")
colormap(np.real(p_data[n]), fig_label="Probability", save_label="prob_T0.001")
colormap(np.real(p_data[-1]), fig_label="Probability", save_label="prob_T0.002")

colormap(np.sqrt(abs(np.real(data[0]))), fig_label="$\sqrt{|{Re(U^n)}|}$", save_label="real_U_T0.0")
colormap(np.sqrt(abs(np.real(data[n]))), fig_label="$\sqrt{|{Re(U^n)}|}$", save_label="real_U_T0.001")
colormap(np.sqrt(abs(np.real(data[-1]))), fig_label="$\sqrt{|{Re(U^n)}|}$", save_label="real_U_T0.002")

colormap(np.sqrt(abs(np.imag(data[0]))), fig_label="$\sqrt{|{Im(U^n)}|}$", save_label="imag_U_T0.0")
colormap(np.sqrt(abs(np.imag(data[n]))), fig_label="$\sqrt{|{Im(U^n)}|}$", save_label="imag_U_T0.001")
colormap(np.sqrt(abs(np.imag(data[-1]))), fig_label="$\sqrt{|{Im(U^n)}|}$", save_label="imag_U_T0.002")


# Plot probability at the detector sheet placed at x = 0.8 and animate the output
prob_detectorsheet(data, x_idx, "cond_prob_x0.8.pdf")
if user_input == "y":
    run_animations(p_data, label="problem_8")

p_data, data = initalize_data("problem9_3_slits.bin")
prob_detectorsheet(data, x_idx, "cond_prob_x0.8_3slits.pdf")
if user_input == "y":
    run_animations(p_data, label="problem_9_threeslits")

p_data, data = initalize_data("problem9_1_slit.bin")
prob_detectorsheet(data, x_idx, "cond_prob_x0.8_1slit.pdf")
if user_input == "y":
    run_animations(p_data, label="problem_9_oneslit")
