import numpy as np
import matplotlib.pyplot as plt
# PyArmadillo is a Python module that provides a Python interface to the Armadillo C++ linear algebra library.
import pyarma as pya
# Animation module from matplotlib
import matplotlib.animation as animation

# Example of declaring a pyarma matrix
A = pya.mat()
A.load(file)


# Example of making an animation, where fig is from "fig, ax = plt.subplots()", animate declares what to animate
ani = animation.FuncAnimation(fig, animate, range(N), interval=1000)
