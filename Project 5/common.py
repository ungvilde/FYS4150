import matplotlib
from matplotlib.animation import FuncAnimation
import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa
import sys

h = 0.005
T = 0.008
dt = 2.5 * 10**(-5)
N_timesteps = int(T / dt)
n = int(0.001/dt)

cm = 1/2.54

if len(sys.argv) > 1:
    user_input = sys.argv[1]    # (y/any input) input from user
else:
    user_input = None

def load(filename):
    data = pa.cx_cube()
    data.load(filename, pa.arma_binary)
    data = np.array(data, dtype=np.cdouble)
    return data

def colormap(data, fig_label, save_label, h=0.005, save=True):
    x_points = np.arange(0, 1+h, h)
    y_points = np.arange(0, 1+h, h)
    x_min, x_max = x_points[0], x_points[-1]
    y_min, y_max = y_points[0], y_points[-1]

    fig = plt.figure(figsize=(9*cm, 10*cm))
    ax = plt.gca()

    # Create a colour scale normalization according to the max z value in the first frame
    norm = matplotlib.cm.colors.Normalize(vmin=0, vmax=np.max(data))

    # Plot the first frame
    img = ax.imshow(data.T, extent=[x_min, x_max, y_min, y_max], cmap=plt.get_cmap("viridis"), norm=norm)
    # Axis labels
    fontsize = 12
    plt.xlabel("x", fontsize=fontsize)
    plt.ylabel("y", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    # Add a colourbar
    cbar = fig.colorbar(img, ax=ax, location='top', shrink=.8)
    cbar.set_label(fig_label, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    plt.tight_layout()
    if save:
        plt.savefig("colormap_" + save_label + ".pdf")


def animate(data, anim_label, h, dt, T):
    # Set up a 2D xy grid
    x_points = np.arange(0, 1+h, h)
    y_points = np.arange(0, 1+h, h)
    x, y = np.meshgrid(x_points, y_points, sparse=True)

    # Array of time points
    t_points = np.arange(0, T+dt, dt)

    # Some settings
    fontsize = 12
    t_min = t_points[0]
    x_min, x_max = x_points[0], x_points[-1]
    y_min, y_max = y_points[0], y_points[-1]

    # Create figure
    fig = plt.figure()
    ax = plt.gca()

    # Create a colour scale normalization according to the max z value in the first frame
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(data.T[0]))

    # Plot the first frame
    img = ax.imshow(data[0].T, extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

    # Axis labels
    plt.xlabel("x", fontsize=fontsize)
    plt.ylabel("y", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    # Add a colourbar
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label("probability", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

    # Add a text element showing the time
    time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white", 
                        horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

    # Function that takes care of updating the z data and other things for each frame
    def animation(i):
        # Normalize the colour scale to the current frame?
        norm = matplotlib.cm.colors.Normalize(vmin=np.min(data[i].T), vmax=np.max(data[i].T))
        img.set_norm(norm)

        # Update z data
        img.set_data(data[i].T)

        # Update the time label
        current_time = t_min + i * dt
        time_txt.set_text("t = {:.3e}".format(current_time))

        return img

    # Use matplotlib.animation.FuncAnimation to put it all together
    anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(data), 2), repeat=False, blit=0)

    # Run the animation!
    # plt.show()

    # # Save the animation
    anim.save(anim_label, writer="ffmpeg", bitrate=-1, fps=30)

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
