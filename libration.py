"""

https://nickcharlton.net/posts/drawing-animating-shapes-matplotlib.html
"""


import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, sin

import matplotlib.animation as animation





e = 1         # Ellipticité
R = 380000.   # Rayon (km)
t_stop =28       # Time range of the simulation
T = 28        # Periode (days)


def position(t):
    """Return position
    """

    x = R * cos(t/T * 2 * np.pi) 
    y = R * sin(t/T * 2 * np.pi) 
    theta = t/T * 2 * np.pi

    return x, y, theta



# create a time array from 0..t_stop in days
dt = 0.5
t = np.arange(0, t_stop, dt)

y = np.empty((len(t), 3))

for i, _t in enumerate(t):
    y[i,0], y[i,1], y[i,2] = position(_t)


fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(autoscale_on=False, xlim=(-2*R, 2*R), ylim=(-2*R, 2*R))
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
circle = plt.Circle((0,0), R/10, color='r')
ax.add_patch(circle)
time_template = 'time = {} days'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)



def init():
    circle.center = (0., 0.)
    # ax.add_patch(circle)
    return line, circle, time_text

def animate(i):
    thisx = [0, y[i,0]]
    thisy = [0, y[i,1]]

    line.set_data(thisx, thisy)
    circle.center
    circle.center = (y[i,0], y[i,1])
    # circle.radius = R/10
    time_text.set_text(time_template.format(i*dt))
    return line, circle, time_text


ani = animation.FuncAnimation(
    fig, animate, 
#    init_func=init, 
    frames=
    len(y), interval=dt*100, blit=True)
plt.show()
