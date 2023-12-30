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



class Libration():
    """Calcule les paramètres de l'orbite lunaire pour mettre en évidence la libration
       - orbite elliptique
       - inclinaison de l'orbite par rapport à l'écliptique et inclinaison de l'axe de rotation lunaire
       - calcul du paralaxe entre observation à l'aube ou au crépuscule

       On peut ainsi mettre en évidence, indépendamment et en exagérant les paramètres :
       - la libration en longitude
       - la libration en latitude
       - la libration de paralaxe

       par contre, on ne mettra pas en évidence la libration physique


    """
    def __init__(self):
        self.time = 0.

        # Paramètres du problème
        self.ellipticity = 0.05490*10               # Ellipticité de l'orbite lunaire
        self.a = 384399                     # demi grand axe de l'orbite lunaire (km)
        self.orbit_inclination = 5.145        # Inclinaison de l'orbite lunaire par rapport au plan de l'écliptique (°)
        self.earth_axis_inclination = 27    # Inclinaison de l'axe de rotation de la Terre par rapport au plan de l'écliptique
        self.moon_axis_inclination = 6.68    # Inclinaison de l'axe de rotation de la Lune par rapport au plan de son orbite autour de la Terre
        self.T = 27.321                       # Période de rotation (sidérale ?) de la Lune
        self.earth_radius = 12000/2*5           # Rayon de la Terre (km)
        self.moon_radius = 1736.*50            # Rayon de la Lune (km)

        # Paramètres d'état
        self._self_rotation_angle = 0. # Angle de rotation propre
        self._theta = 0.               # coordonnée angulaire du référentiel géocentrique
        self._omega = 0.               # paramètre de rotation propre de la Lune
        self._r = 0.                   # coordonnée de distance du référenciel géocentrique



    def set_time(self,time):
        self.time = time

    def get_distance(self, theta=None):
        """Renvoi la distance Terre-Lune (r)
        
        see https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion

        """

        p = self.a*(1-self.ellipticity**2)

        if theta==None:
            theta = self.time/self.T *2 * np.pi
        self._r = p/(1+self.ellipticity * np.cos(theta))

        return self._r

    def get_orbit_XYZ(self, theta=None):
        """Renvoi les coordonnées cartésiennes de la lune dans le référentiel geocentrique.

        note : le référentiel géocentrique -> e_x et et e_y sont dans le plan de l'écliptique 

        """


        if theta == None :
            theta = self.time/self.T * 2 * np.pi

        r = self.get_distance(theta)
        x = r * cos(theta) 
        y = r * sin(theta) 
        return x, y, None

    
    


class LibrationPlot(Libration):
    """Complète la classe Libration avec des fonction qui tracent les éléments à afficher
    """

    def __init__(self, axis):
        Libration.__init__(self)
        self.top_orbit_axis, self.moon_face_axis = axis

        # Disque lunaire
        self.top_moon_disc = plt.Circle((self.a, 0), self.moon_radius, color='y')
        self.top_orbit_axis.add_patch(self.top_moon_disc)
        self.plot_moon_orbit_lines()
    
    def plot_top_moon_disc(self):
        """Trace le disque lunaire (à updater)
        """
        x, y, z = self.get_orbit_XYZ()
        self.top_moon_disc.center = (x, y)
        return [self.top_moon_disc]

    def plot_top_moon_meridians(self):
        """ Trace les méridiens de la Lune vu de dessus (à updater)

        pour l'instant j'en trace 6 (3 lignes)
        """

        x, y, z = self.get_orbit_XYZ()

        pix2 = 2 * np.pi
        meridian_angles = [0, -45, 45, 90]
        meridian_angles = [a/180*np.pi for a in meridian_angles ]

        try:
            self.top_moon_meridians
        except AttributeError:
            self.top_moon_meridians = []
            for i, pos_i in enumerate(meridian_angles):
                moon_meridian, = self.top_orbit_axis.plot([], [], '-', lw=2., color='k')
                self.top_moon_meridians.append(moon_meridian)


        for i, angle in enumerate(meridian_angles):
            print(x,y, angle)
            pt_1 = [
                    -1*self.moon_radius*np.cos(angle + self._omega) + x, 
                       self.moon_radius*np.cos(angle + self._omega) + x,    
                    ]
            pt_2 = [
                    -1*self.moon_radius*np.sin(angle + self._omega) + y,
                       self.moon_radius*np.sin(angle + self._omega) + y
                       ]
            self.top_moon_meridians[i].set_data(pt_1, pt_2) 
            

        return self.top_moon_meridians

    def plot_moon_orbit_lines(self):
        """Trace l'orbite elliptique de la Lune (statique)
        """
        self.moon_orbit_lines, = self.top_orbit_axis.plot([], [], '--', lw=0.5, color='k')
        pos = np.empty((1000, 2))
        for i, pos_i in enumerate(pos):
            x, y, z = self.get_orbit_XYZ(i*np.pi*2/1000)
            pos[i,0] = x
            pos[i,1] = y
        self.moon_orbit_lines.set_data(pos[:,0], pos[:,1])


    def plot_moon_face(self):
        """Trace la lune vu de face : disque et méridiens
        """

        try:
            self.moon_face_disc
        except AttributeError:
            self.moon_face_disc = plt.Circle((0, 0), 1, color='y')
            self.moon_face_axis.add_patch(self.moon_face_disc)

        self.get_distance()
        self.moon_face_disc.radius = (0.5*self._r/self.a)
        return [self.moon_face_disc]
    
        


    def animate(self, time):

        self.set_time(time)

        # Plot de l'orbite vu de dessus
        self.plot_top_moon_disc()
        self.plot_top_moon_meridians()
        self.plot_moon_face()
        # return top_moon_disc
        # Must return an iterable
        # print(self.top_moon_disc + self.top_moon_meridians)
        return [self.top_moon_disc] + self.top_moon_meridians + [self.moon_face_disc]


libration_parameters = Libration()

fig = plt.figure(figsize=(5, 4))

# fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
ax1 = fig.add_subplot(221, autoscale_on=False, xlim=(-2*libration_parameters.a, 2*libration_parameters.a), ylim=(-2*libration_parameters.a, 2*libration_parameters.a))
ax2 = fig.add_subplot(222, autoscale_on=False, xlim=(-1, 1), ylim=(-1, 1))
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax1.grid()

top_earth_disc = plt.Circle((0, 0), libration_parameters.earth_radius, color='b')
ax1.add_patch(top_earth_disc)


libration = LibrationPlot((ax1, ax2))
# create a time array from 0..t_stop in days
dt = 0.0005
# t = np.arange(0, libration_parameters.T, dt)
# t = np.arange(0, libration_parameters.T/1000, dt)
t = np.arange(0, 360, dt)

print(f'nb points = {len(t)}')


ani = animation.FuncAnimation(
    fig, libration.animate, frames= len(t), interval=dt*100000, blit=True)
plt.show()
