import numpy as np
import matplotlib.pyplot as plt

from constants import *


class SpurGear:

    def __init__(self, N, Pd, phi, gear_type, ts=None, a=None, b=None):
        # Values  at the standard pitch cicle
        self.N = N                  # Number of teeth
        self.Pd = Pd                # Diametral pitch
        self.phi = phi              # Pressure angle
        self.gear_type = gear_type  # Gear type (int/ext)
        if ts is not None:
            self.ts = ts            # Tooth thickness
        else:
            self.ts = self.ps() / 2
        if a is not None:
            self.a = a              # Addendum
        else:
            self.a = 1 / Pd
        if b is not None:
            self.b = b              # Dedendum
        else:
            self.b = 1.25 / Pd
        # Operating pitch radius and pressure angle
        # Any associated operating values may be calculated using them
        self.R_op = None
        self.phi_op = None

    def __repr__(self):
        if self.gear_type == EXTERNAL:
            g = 'External'
        else:
            g = 'Internal'
        r = '{} Spur Gear\n'.format(g)
        r += '\tNumber of Teeth: {}\n'.format(self.N)
        r += '\tStandard Pitch Diameter: {}\n'.format(2*self.Rs())
        r += '\tBase Diameter: {}\n'.format(2*self.Rb())
        r += '\tDiametral Pitch: {}\n'.format(self.Pd)
        r += '\tModule: {}\n'.format(1/self.Pd)
        r += '\tPressure Angle: {}\n'.format(np.degrees(self.phi))
        r += '\tTooth Thickness: {}\n'.format(self.ts)
        r += '\tAddendum Diameter: {}\n'.format(2*self.Ra())
        r += '\tDedendum Diameter: {}\n'.format(2*self.Rd())
        if self.R_op is not None:
            r += '\tOperating Pitch Diameter: {}\n'.format(2*self.R_op)
        return r

    def set_operating_point(self, R_op, phi_op):
        self.R_op = R_op
        self.phi_op = phi_op
        
    def Rs(self):
        # Standard pitch radius
        return self.N / (2*self.Pd)

    def Rb(self):
        # Base radius
        return self.Rs() * np.cos(self.phi)

    def ap(self):
        # Angular pitch
        return 2*np.pi / self.N

    def ps(self):
        # Standard pitch
        return 2*np.pi*self.Rs() / self.N

    def pb(self):
        # Base pitch
        return 2*np.pi*self.Rb() / self.N
		
    def rho_p(self):
        # Radius of curvature to pitch point
        return self.Rb()*np.tan(self.phi_op)
            
    def rho_a(self):
        # Radius of curvature to addendum diameter
        return np.sqrt(self.Ra()**2 - self.Rb()**2)

    def phi_R(self, R):
        # Pressure angle at specified radius
        return np.arccos(self.Rb() / R)

    def R_phi(self, phi):
        # Radius at specified pressure angle
        return self.Rb() / np.cos(phi)

    @staticmethod
    def inv(phi):
        # Involute function
        return np.tan(phi) - phi

    @staticmethod
    def invinv(x):
        # Inverse involute function
        q = x**(2/3)
        p = 1 + 1.04004*q + .32451*q**2 - .00321*q**3 - \
            .00894*q**4 + .00319*q**5 - .00048*q**6
        return np.arccos(1/p)


class ExternalSpurGear(SpurGear):

    def __init__(self, N, Pd, phi, ts=None, a=None, b=None):
        super().__init__(N, Pd, phi, EXTERNAL, ts, a, b)

    def Ra(self):
        # Addendum (tip/outer/major) radius
        return self.Rs() + self.a

    def Rd(self):
        # Dedendum (root/minor) radius
        return self.Rs() - self.b

    def a_op(self):
        return self.a + self.Rs() - self.R_op

    def b_op(self):
        return self.b + self.R_op - self.Rs()

    def t_R(self, R):
        # Tooth thickness at specified radius
        inv_phiS = self.inv(self.phi)
        inv_phiR = self.inv(self.phi_R(R))
        return R*(self.ts/self.Rs() + 2*(inv_phiS - inv_phiR))

    def th_R(self, R):
        # Polar angle at specified radius
        return self.t_R(R) / (2*R)

    def R_th(self, th):
        # Radius at specified polar angle
        inv_phiS = self.inv(self.phi)
        inv_phiR = self.ts / (2*self.Rs()) + inv_phiS - th
        phiR = self.invinv(inv_phiR)
        return self.R_phi(phiR)

    def animate(self, xC=0, yC=0, theta0=0):
        # Animate one rotation of gear in specified position
        nframes = 100
        fig, ax = plt.subplots()
        for i in range(nframes):
            theta = theta0 + 2*np.pi * (i/nframes)
            self.draw(ax, xC, yC, theta)
            plt.pause(.0001)
            ax.clear()
        plt.show()
        plt.close(fig)

    def draw(self, ax, xC=0, yC=0, theta0=0):
        # Draw gear in specified position
        xg, yg = self._gear_coordinates(xC, yC, theta0)
        xL, yL = self._centerline(xC, yC, theta0)
        ax.plot(xg, yg, color='black')
        ax.plot(xL, yL, color='blue')
        ax.set_xlim([-1.25*self.Ra() + xC, 1.25*self.Ra() + xC])
        ax.set_ylim([-1.25*self.Ra() + yC, 1.25*self.Ra() + yC])
        ax.set_aspect('equal')

    def _centerline(self, xC, yC, theta0):
        npoints = 10
        xL = np.linspace(xC, xC + self.Ra()*np.cos(theta0), npoints)
        yL = np.zeros(npoints)
        for i in range(npoints):
            yL[i] = yC + (xL[i] - xC)*np.tan(theta0)
        return xL, yL

    def _gear_coordinates(self, xC, yC, theta0):
        # Draw gear on screen
        npoints = 50
        x0, y0 = self._tooth0_coordinates()
        xg, yg = [], []
        for i in range(self.N):
            theta = theta0 + 2*np.pi * (i / self.N)
            xt, yt = self._tooth_coordinates(x0, y0, theta)
            for j in range(npoints):
                xg.append(xC + xt[j])
                yg.append(yC + yt[j])
        return xg, yg

    def _tooth_coordinates(self, x0, y0, theta):
        # Take 0deg tooth coordinates and rotate to tooth position
        npoints = 50
        c = np.cos(theta)
        s = np.sin(theta)
        xt, yt = np.zeros(npoints), np.zeros(npoints)
        for i in range(npoints):
            xt[i] = c*x0[i] - s*y0[i]
            yt[i] = s*x0[i] + c*y0[i]
        return xt, yt 

    def _tooth0_coordinates(self):
        # Compute basic coordinates for one tooth at 0deg position
        npoints = 50
        theta = np.linspace(-self.ap()/2, self.ap()/2, npoints)
        x0, y0 = np.zeros(npoints), np.zeros(npoints)
        tha = self.th_R(self.Ra())
        if self.Rd() < self.Rb():
            thd = self.th_R(self.Rb())
        else:
            thd = self.th_R(self.Rd())
        # Positive polar angle points
        for i in range(npoints // 2, npoints):
            th = theta[i]
            if th < tha:
                R = self.Ra()
            elif th < thd:
                R = self.R_th(th)
            else:
                R = self.Rd()
            x0[i] = R * np.cos(th)
            y0[i] = R * np.sin(th)
        # Mirror points
        for i in range(0, npoints // 2):
            x0[i] = x0[npoints - i - 1]
            y0[i] = -y0[npoints - i - 1]
        return x0, y0


class Tests:

    def test1(self):
        g = ExternalSpurGear(Pd=4, N=16, phi=np.radians(25), a=.375, ts=.5091)
        x0, y0 = g._tooth0_coordinates()
        plt.plot(x0, y0)
        plt.xlim([-1.25*g.Ra(), 1.25*g.Ra()])
        plt.ylim([-1.25*g.Ra(), 1.25*g.Ra()])
        plt.show()

    def test2(self):
        g = ExternalSpurGear(Pd=4, N=16, phi=np.radians(20))
        fig, ax = plt.subplots()
        g.draw(ax, xC=2, yC=3, theta0=np.radians(5))
        plt.show()

    def test3(self):
        g = ExternalSpurGear(Pd=4, N=16, phi=np.radians(20))
        g.animate()

if __name__ == '__main__':
    t = Tests()
    t.test2()
