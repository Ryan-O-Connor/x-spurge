import matplotlib.pyplot as plt
import numpy as np


class RackCutter:

    def __init__(self, Pd, phi, a, r):
        self.Pd = Pd
        self.phi = phi
        self.a = a
        self.r = r
        self.p = np.pi / Pd
        self.t = self.p / 2
        self.h = a - r*(1 - np.sin(phi))

    def __repr__(self):
        r = 'Rack Cutter:\n'
        r += '\tDiametral Pitch: {}\n'.format(self.Pd)
        r += '\tPressure Angle: {}\n'.format(self.phi)
        r += '\tAddendum: {}\n'.format(self.a)
        r += '\tFillet Radius: {}\n'.format(self.r)
        r += '\tFillet Height: {}\n'.format(self.h)
        r += '\tTooth Thickness: {}\n'.format(self.t)
        return r

    def rack_coordinates(self):
        npoints = 200
        t, h, a, r, phi = self.t, self.h, self.a, self.r, self.phi
        xr = np.linspace(-t, t, npoints)
        yr = np.zeros(npoints)
        for i, x in enumerate(xr):
            if x < -t/2:
                # Rack flank, below reference line
                dx = -t/2 - x
                yr[i] = -dx/np.tan(phi)
            elif x < -t/2 + h*np.tan(phi):
                # Rack flank, above reference line
                dx = t/2 + x
                yr[i] = dx / np.tan(phi)
            elif x < -t/2 + h*np.tan(phi) + r*np.cos(phi):
                # Rack fillet
                xC = -t/2 + h*np.tan(phi) + r*np.cos(phi)
                yC = a - r
                dx = xC - x
                dy = np.sqrt(r**2 - dx**2)
                yr[i] = yC + dy
            elif x <= 0:
                # Rack addendum
                yr[i] = a
            else:
                # Reflect coordinates about y axis
                yr[i] = yr[npoints - i - 1]
        return xr, yr

    def cut_coordinates(self, Rg):
        npoints = 200
        xr, yr = self.rack_coordinates()
        xg, yg = np.zeros(npoints), np.zeros(npoints)
        for i in range(npoints):
            x = yr[i]
            y = xr[i]
            xi = x
            eta = xi / np.tan(self.phi)
            ur = eta - y
            beta = (ur - self.p / 2) / Rg
            R = np.sqrt((Rg + xi)**2 + eta**2)
            th = np.arctan2(eta, Rg + xi) - beta
            xg[i] = R*np.cos(th)
            yg[i] = R*np.sin(th)
        return xg, yg


if __name__ == '__main__':
    Pd = 16
    phir = np.radians(14.5)
    ar = 1.25 / Pd
    rr = 0.30 / Pd
    rack = RackCutter(Pd, phir, ar, rr)
    print(rack)
    fig, ax = plt.subplots()
    # xr, yr = rack.rack_coordinates()
    Rg = 36 / 2 / Pd
    xg, yg = rack.cut_coordinates(Rg)
    ax.plot(xg, yg)
    plt.grid()
    plt.show()


