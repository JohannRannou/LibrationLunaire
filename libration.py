"""

https://nickcharlton.net/posts/drawing-animating-shapes-matplotlib.html
"""


import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.animation as animation

from numpy import cos, sin
from numpy import cos
from numpy import sin
from numpy import pi



TWOPI = 2 * pi



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
        self.b = self.a*np.sqrt(1-self.ellipticity**2)                     # demi petit axe de l'orbite lunaire (km)
        self.c = np.sqrt(self.a**2 - self.b**2)                            # distance centre-foyer (km)
        self.earth_axis_inclination = 23.44    # Obliquité : inclinaison de l'axe de rotation de la Terre par rapport au plan de l'écliptique
        self.orbit_inclination = 5.145       # Inclinaison de l'orbite lunaire par rapport au plan de l'écliptique (°)
        self.moon_axis_inclination = 6.68   # Inclinaison de l'axe de rotation de la Lune par rapport au plan de son orbite autour de la Terre
        self.alpha = (self.orbit_inclination - self.moon_axis_inclination )  # angle en degrés entre l'équateur de la Lune et le plan de l'écliptique
        self.T = 27.321                       # Période de rotation (sidérale ?) de la Lune
        self.earth_radius = 12000/2*5           # Rayon de la Terre (km)
        self.moon_radius = 1736.*50            # Rayon de la Lune (km)

        # Paramètres d'état
        self._theta = 0.               # coordonnée angulaire du référentiel géocentrique
        self._omega = 0.               # paramètre de rotation propre de la Lune
        self._r = 0.                   # coordonnée de distance du référenciel géocentrique
        self._x = 0.                   # x,y,z : 3 paramètres cartésien dans le référentiel géocentrique
        self._y = 0.
        self._z = 0.
        self._face_radius = 1.

        self._modify()

    def _modify(self):
        """Modify parameters to highlight a specific phenomenon
        """
        self.orbit_inclination = 10       # Inclinaison de l'orbite lunaire par rapport au plan de l'écliptique (°)
        self.moon_axis_inclination = 20   # Inclinaison de l'axe de rotation de la Lune par rapport au plan de son orbite autour de la Terre
        self.T = 10.
        # self.orbit_inclination = 20
        # Recalculé ici au cas ou les paramètres seraient modifiés ici
        self.alpha = (self.orbit_inclination - self.moon_axis_inclination)
    
    def set_time(self,time):
        """
        time in days
        """
        self.time = time


    def _compute_rotation_matrices(self):
        """Calcule les matrices de rotation
        
        Voir notes manuscripts.
        """

        self.P_l_lp = np.eye(3)
        self.P_lp_o = np.eye(3)
        self.P_o_theta = np.eye(3)

        self.P_l_lp[1,1] = cos(-self._omega + pi/2)
        self.P_l_lp[2,1] = sin(-self._omega + pi/2)
        self.P_l_lp[1,2] = -sin(-self._omega + pi/2)
        self.P_l_lp[2,2] = cos(-self._omega + pi/2)

        self.P_lp_o[0,0] = cos(-self.alpha)
        self.P_lp_o[1,0] = sin(-self.alpha)
        self.P_lp_o[0,1] = -sin(-self.alpha)
        self.P_lp_o[1,1] = cos(-self.alpha)

        self.P_o_theta[1,1] = cos(self._theta)
        self.P_o_theta[2,1] = sin(self._theta)
        self.P_o_theta[1,2] = -sin(self._theta)
        self.P_o_theta[2,2] = cos(self._theta)

    def update(self, time):
        """
        time in days
        """
        self.set_time(time)
        self._theta = self.time / self.T * 2 * np.pi
        self.get_distance()
        self.compute_orbit_XYZ()
        self.compute_omega()
        self.compute_beta()
        self.compute_face_radius()
        self._compute_rotation_matrices()
        print(f'at day {self.time} : theta = {self._theta*360/TWOPI}, omega = {self._omega*360/TWOPI}, r = {self._r}')

    def get_distance(self, theta=None):
        """Renvoi la distance Terre-Lune (r)
        
        see https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion

        """

        p = self.a*(1-self.ellipticity**2)

        if theta==None:
            theta = self.time/self.T *2 * np.pi
        self._r = p/(1+self.ellipticity * np.cos(theta))

        return self._r

    def compute_face_radius(self):
        """Retourne le diamètre apparent de la Lune en °
        """
        # self._face_radius = self.moon_radius*self.a / self._r*360/TWOPI
        self._face_radius = 0.25 * self.a / self._r

    def compute_omega(self):
        # Pour l'instant c'est faux, ça pointe vers le centre de l'ellipse

        self._omega = np.arctan2(self._y , (self._x + self.a*self.ellipticity))

    def compute_beta(self):

        self._omega = np.arctan2(self._y , (self._x + self.a*self.ellipticity))
        


    def compute_orbit_XYZ(self, theta=None):
        """Renvoi les coordonnées cartésiennes de la lune dans le référentiel geocentrique.

        note : le référentiel géocentrique -> e_x et et e_y sont dans le plan de l'écliptique 

        """


        if theta == None :
            theta = self.time/self.T * 2 * np.pi

        self.get_distance(theta)
        self._x = self._r * cos(theta) 
        self._y = self._r * sin(theta) 
        self._z = self._x * np.tan(np.radians(self.orbit_inclination)) 

        return self._x, self._y, self._z 

    
    


class LibrationPlot(Libration):
    """Complète la classe Libration avec des fonction qui tracent les éléments à afficher


    Il y a plusieurs vues :

    

        * vue de face


        * vue de dessus






    """

    def __init__(self, axis):
        Libration.__init__(self)
        self.top_orbit_axis, self.moon_face_axis, self.front_orbit_axis = axis

        # Disque lunaire
        self.top_moon_disc = plt.Circle((self.a, 0), self.moon_radius, color='gold', alpha=0.75)
        self.top_orbit_axis.add_patch(self.top_moon_disc)
        self.plot_moon_orbit_segments()
    
        self.meridian_angles = np.array([180, 180-45, 90]) * np.pi/180
        self.meridian_colors = ['r', 'b', 'k']

        self.parallel_angles = np.array([180, 180-45, 90]) * np.pi/180
        self.parallel_colors = ['r', 'b', 'k']


    def plot_top_moon_disc(self):
        """Trace le disque lunaire (à updater)
        """
        x, y, z = self.compute_orbit_XYZ()
        self.top_moon_disc.center = (x, y)
        return [self.top_moon_disc]


    def plot_front_moon_disc(self):
        """Trace le disque lunaire (à updater)
        """
        try:
            self.front_moon_disc
        except:
            self.front_moon_disc = plt.Circle((self.a, 0), self.moon_radius, color='gold', alpha=0.75)
            self.front_orbit_axis.add_patch(self.front_moon_disc)

        self.front_moon_disc.center = (self._x, self._z)


  
    def plot_moon_orbit_segments(self):
        """Trace l'orbite elliptique de la Lune (statique)
        sur les vues de face et de dessus
        """
        self.moon_orbit_segments, = self.top_orbit_axis.plot([], [], '--', lw=0.5, color='k')
        N=100
        pos = np.empty((N, 2))
        for i, pos_i in enumerate(pos):
            x, y, z = self.compute_orbit_XYZ(i*np.pi*2/N)
            pos[i,0] = x
            pos[i,1] = y
        self.moon_orbit_segments.set_data(pos[:,0], pos[:,1])



        self.front_moon_orbit_segments, = self.front_orbit_axis.plot([-(self.a+self.c), (self.a-self.c)], [-(self.a+self.c)*np.tan(np.radians(self.orbit_inclination)), +(self.a-self.c)*np.tan(np.radians(self.orbit_inclination))], '--', lw=0.5, color='k')

    def plot_moon_face(self):
        """Trace la lune vu de face : disque et méridiens
        """

        try:
            self.moon_face_disc
        except AttributeError:
            self.moon_face_disc = plt.Circle((0, 0), self._face_radius, color='gold', alpha=0.75)
            self.moon_face_axis.add_patch(self.moon_face_disc)

        self.moon_face_disc.radius = (self._face_radius)
        return [self.moon_face_disc]
    

    def plot_moon_face_meridians(self):
        """Trace les méridiens de la Lune vu de la Terre (à updater)

        pour l'instant j'en trace 6 (3 lignes)
        """

        TWOPI = 2 * np.pi

        try:
            self.moon_face_meridians
        except AttributeError:
            self.moon_face_meridians = []
            for i, pos_i in enumerate(self.meridian_angles):
               moon_face_meridian = mpatches.Arc((0, 0), 2 * self._face_radius, 0., color=self.meridian_colors[i])
               self.moon_face_axis.add_patch(moon_face_meridian)
               self.moon_face_meridians.append(moon_face_meridian)
        for i, meridian in enumerate(self.moon_face_meridians):
            meridian.width = -2 *self._face_radius*np.sin(self._omega - self._theta + self.meridian_angles[i])
            meridian.height = 2 * self._face_radius
            # meridian.theta1 = 0 + meridian_angles[i]*360/TWOPI
            # meridian.theta2 = 180 + meridian_angles[i]*360/TWOPI

            if np.sin(self._omega - self._theta + self.meridian_angles[i]) < 0 :
                meridian.theta1 = -90 + 1e-6 
                meridian.theta2 =  90 - 1e-6
            else:
                meridian.theta1 = 90 - 1e-6 
                meridian.theta2 = 270 - 1e-6

            # Ce qui suit semble ne pas fonctionner
            if np.cos(self._omega - self._theta + self.meridian_angles[i]) < 0 :
                meridian.linestyle = '--'
                meridian.linewidth = 2
            else:
                meridian.linestyle = '--'
                meridian.linewidth = 0.2

            # meridian.theta1 = -90 - 0.000001 
            # meridian.theta2 =  90 -0.000001


        return self.moon_face_meridians



    def plot_moon_face_equator(self):
        """Trace l'équateur de la Lune vu de la Terre (à updater)

        La projection du cercle de l'équateur est assez compliqué : c'est un ellipse qui change d'axe entre deux extrèmes :
            - une ellipse la plus ouverte "horizontale"
            - une ellipse fermée (un segment droit) incliné de l'angle alpha

        """
        try: 
            self.equator_segments
        except:
            self.equator_segments, = self.moon_face_axis.plot([], [], '-', lw=2., color='c')
            # On dessine N segments
        N = 30
        t_list = np.linspace( 0, 2 * np.pi, N)
        pos = np.empty((N, 2))
        for i, t in enumerate(t_list):
            # pos[i,0] = self._face_radius * (np.cos(self._theta-np.pi/2) * np.cos(self.alpha) * np.cos(t) + np.sin(self._theta-np.pi/2)*np.sin(t))
            # pos[i,1] = self._face_radius * (-np.sin(self.alpha) * np.cos(t))

            X_l = self._face_radius * np.array([0, cos(t), sin(t)])

            X_theta = self.P_o_theta.T @ self.P_lp_o.T @ self.P_l_lp.T @ X_l

            pos[i,0] = X_theta[2]
            pos[i,1] = X_theta[0]

        self.equator_segments.set_data(pos[:,0], pos[:,1])




    def plot_top_moon_meridians(self):
        """ Trace les méridiens de la Lune vu de dessus (à updater)

        pour l'instant j'en trace 6 (3 lignes)
        """
 
        try:
            self.top_moon_meridians
        except AttributeError:
            self.top_moon_meridians = []
            for i, pos_i in enumerate(self.meridian_angles):
                moon_meridian, = self.top_orbit_axis.plot([], [], '-', lw=1, color=self.meridian_colors[i])
                self.top_moon_meridians.append(moon_meridian)

        # for i, angle in enumerate(self.meridian_angles):
        #     pt_1 = [
        #                 x,    
        #                 self.moon_radius*np.cos(angle + self._omega) + x, 
        #             ]
        #     pt_2 = [
        #                 y,
        #                self.moon_radius*np.sin(angle + self._omega) + y,
        #                ]
        #     self.top_moon_meridians[i].set_data(pt_1, pt_2) 
            

        return self.top_moon_meridians




    def plot_front_moon_parallels(self):
        """ Trace les parallèles de la Lune vu de dessus (à updater)

        pour l'instant je ne trace que l'équateur et ça devrait suffire)

        Je rajouter l'axe de rotation aussi, ce n'est pas un parallèle certes ...
        """

        x, y, z = self.compute_orbit_XYZ()
       
        try:
            self.front_moon_parallels
        except:
            self.front_moon_parallels = []
            moon_equator, = self.front_orbit_axis.plot([], [], '-', lw=1, color = 'k')
            self.front_moon_parallels.append(moon_equator)
            
            moon_rotation_axis, = self.front_orbit_axis.plot([], [], '--', lw=0.5
                                                             , color = 'k')
            self.front_moon_parallels.append(moon_rotation_axis)


        print(self.alpha)
        pt_1 = [x - self.moon_radius * np.cos(np.radians(self.alpha)),
                x + self.moon_radius * np.cos(np.radians(self.alpha)),
                ]
        
        pt_2 = [z - self.moon_radius * np.sin(np.radians(self.alpha)),
                z + self.moon_radius * np.sin(np.radians(self.alpha)),
                ]
        self.front_moon_parallels[0].set_data(pt_1, pt_2)

        pt_3 = [x - self.moon_radius * 1.4 * np.sin(np.radians(self.alpha)),
                x + self.moon_radius * 1.4 * np.sin(np.radians(self.alpha)),
                ]
        
        pt_4 = [z + self.moon_radius * 1.4 * np.cos(np.radians(self.alpha)),
                z - self.moon_radius * 1.4 * np.cos(np.radians(self.alpha)),
                ]
        self.front_moon_parallels[1].set_data(pt_3, pt_4)



        return self.front_moon_parallels





    def animate(self, time):

        self.update(time)
        # Plot de l'orbite vu de dessus
        self.plot_top_moon_disc()
        self.plot_front_moon_disc()
        self.plot_top_moon_meridians()
        self.plot_moon_face()
        self.plot_moon_face_meridians()
        self.plot_moon_face_equator()

        self.plot_front_moon_parallels()
        self.plot_moon_face_equator()

        # Must return an iterable
        return [self.top_moon_disc] +[self.front_moon_disc] + self.top_moon_meridians + [self.moon_face_disc] + self.moon_face_meridians + self.front_moon_parallels + [self.equator_segments]


libration_parameters = Libration()

fig = plt.figure(figsize=(10, 10))

# fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
ax1 = fig.add_subplot(223, autoscale_on=False, xlim=(-2*libration_parameters.a, 2*libration_parameters.a), ylim=(-2*libration_parameters.a, 2*libration_parameters.a))
ax2 = fig.add_subplot(224, autoscale_on=False, xlim=(-1, 1), ylim=(-1, 1))
ax3 = fig.add_subplot(221, autoscale_on=False, xlim=(-2*libration_parameters.a, 2*libration_parameters.a), ylim=(-2*libration_parameters.a, 2*libration_parameters.a))
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax1.grid()
ax2.grid()
ax3.grid()


top_earth_disc = plt.Circle((0, 0), libration_parameters.earth_radius, color='b')
ax1.add_patch(top_earth_disc)
front_earth_disc = plt.Circle((0, 0), libration_parameters.earth_radius, color='b')
ax3.add_patch(front_earth_disc)


libration = LibrationPlot((ax1, ax2, ax3))
# create a time array from 0..t_stop in days
dt = 0.1
t = np.arange(0, libration_parameters.T, dt)
# t=[7.5]
print(t)
print(f'nb points = {len(t)}')


ani = animation.FuncAnimation(
    fig, libration.animate, frames= t, interval=50, blit=True)
plt.show()
