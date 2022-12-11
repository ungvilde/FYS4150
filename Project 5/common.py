import matplotlib
from matplotlib.animation import FuncAnimation
import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

def load(filename):
    data = pa.cx_cube()
    data.load(filename, pa.arma_binary)
    data = np.array(data, dtype=np.cdouble)
    return data

def animate(data, h, dt, T):

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
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(data[0]))

    # Plot the first frame
    img = ax.imshow(data[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

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
        norm = matplotlib.cm.colors.Normalize(vmin=np.min(data[i]), vmax=np.max(data[i]))
        img.set_norm(norm)

        # Update z data
        img.set_data(data[i])

        # Update the time label
        current_time = t_min + i * dt
        time_txt.set_text("t = {:.3e}".format(current_time))

        return img

    # Use matplotlib.animation.FuncAnimation to put it all together
    anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(data), 2), repeat=False, blit=0)

    # Run the animation!
    plt.show()

    # # Save the animation
    # anim.save('./animation.mp4', writer="ffmpeg", bitrate=-1, fps=30)
