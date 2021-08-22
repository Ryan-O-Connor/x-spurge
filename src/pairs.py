import numpy as np
import matplotlib.pyplot as plt

import gears


class GearPair:

    def __init__(self, pinion, gear, C):
        self.pinion = pinion
        self.gear = gear
        if type(C) == type('STANDARD') and C.upper() == 'STANDARD':
            self.C = self.Cs()
        else:
            self.C = C
        self.set_operating_points()

    def __repr__(self):
        r = 'Gear Pair:\n'
        r += '***Pinion***\n'
        r += self.pinion.__repr__()
        r += '***Gear***\n'
        r += self.gear.__repr__()
        r += '***Mesh Geometry***\n'
        r += '\tStandard Center Distance: {}\n'.format(self.Cs())
        r += '\tOperating Center Distance: {}\n'.format(self.C)
        r += '\tOperating Pressure Angle: {}\n'.format(np.degrees(self.phi_op()))
        r += '\tBacklash: {}\n'.format(self.backlash())
        r += '\tContact Ratio: {}\n'.format(self.mc())
        return r

    def checks(self):
        # Check against each type of interference
        self.center_distance_check()
        self.root_interference_check()
        self.fillet_interference_check()
        self.working_depth_check()
        self.tip_thickness_check()


class ExternalSpurGearPair(GearPair):

    def animate(self, area):
        # Animate one pinion revolution
        nframes = 100
        fig, ax = plt.subplots()
        for i in range(nframes):
            theta_p = 2*np.pi * (i / nframes)
            self.pinion.draw(ax, xC=0, yC=0, theta0=theta_p)
            theta_g = self.beta_2(theta_p) - np.pi
            self.gear.draw(ax, xC=self.C, yC=0, theta0=theta_g)
            if area.upper() == 'MESH':
                ax.set_xlim([0.5*self.pinion.R_op, 1.25*self.pinion.R_op])
                ax.set_ylim([-3*self.pinion.ps(), 3*self.pinion.ps()])
            else:
                ax.set_xlim([-1.25*self.pinion.Ra(), 1.25*self.gear.Ra()+self.C])
                ax.set_ylim([-1.25*self.gear.Ra(), 1.25*self.gear.Ra()])
            plt.pause(.0001)
            ax.clear()
        plt.show()
        plt.close(fig)

    def set_operating_points(self):
        R_opp = self.C*self.pinion.N / (self.pinion.N + self.gear.N)
        self.pinion.set_operating_point(R_opp, self.phi_op())
        R_opg = self.C*self.gear.N / (self.pinion.N + self.gear.N)
        self.gear.set_operating_point(R_opg, self.phi_op())

    def Cs(self):
	# Standard center distance
        return self.gear.Rs() + self.pinion.Rs()

    def phi_op(self):
        # Operating pressure angle
        x = (self.pinion.Rb() + self.gear.Rb()) / self.C
        return np.arccos(x)

    def pp(self):
        # Operating pitch
        return 2*np.pi*self.C/(self.pinion.N + self.gear.N)

    def B(self):
        # Pitch minus operating tooth thicknesses
        t_opp = self.pinion.t_R(self.pinion.R_op)
        t_opg = self.gear.t_R(self.gear.R_op)
        b = self.pp() - t_opp - t_opg
        return b

    def beta_2(self, beta_1):
        # Position of gear given position of pinion
        R_op1 = self.pinion.R_op
        R_op2 = self.gear.R_op
        t_op1 = self.pinion.t_R(R_op1)
        t_op2 = self.gear.t_R(R_op2)
        return (-0.5*(t_op1 + t_op2) - R_op1*beta_1) / R_op2

    def mc(self):
        # Profile contact ratio
        Z = self.pinion.rho_a() + self.gear.rho_a() - (self.pinion.rho_p() + self.gear.rho_p())
        mc = Z / self.pinion.pb()
        return mc

    def center_distance_check(self):
        print('***Center Distance Check***')
        if self.Cs() <= self.C <= self.Cs() + np.sqrt(self.C*self.B()*np.tan(self.phi_op())):
            print('PASS')
        else:
            print('FAIL')

    def root_clearance_check(self):
        print('***Root Clearance Check***')
        c_p = self.pinion.b_op() - self.gear.a_op()
        c_g = self.gear.b_op() - self.pinion.a_op()
        print('\tPinion root clearance: {}'.format(c_p))
        print('\tGear root clearance: {}'.format(c_g))
        print('\tStandard clearance: {}'.format(0.25/self.pinion.Pd))

    def base_clearance_check(self):
        # Check if contact zone interferes with base circle
        print('***Base Interference Check***')
        if np.sqrt(self.gear.Ra()**2 - self.gear.Rb()**2) < (self.pinion.Rb() + self.gear.Rb())*np.tan(self.phi_op()):
            if np.sqrt(self.pinion.Ra()**2 - self.pinion.Rb()**2) < (self.pinion.Rb() + self.gear.Rb())*np.tan(self.phi_op()): 
                print('PASS')
            else:
                print('FAIL: Interference at gear base circle')
        else:
            print('FAIL: Interference at pinion base circle')

    def fillet_interference_check(self):
        # Check if contact zone interferes with fillet
        print('***Fillter Interference Check***')
        RLp = np.sqrt(self.pinion.Rb()**2 + \
                ((self.pinion.Rb() + self.gear.Rb())*np.tan(self.phi_op()) \
                - np.sqrt(self.gear.Ra()**2 - self.gear.Rb()**2)))
        Rfp = self.pinion.Rf()
        RLg = np.sqrt(self.gear.Rb()**2 + \
                ((self.pinion.Rb() + self.gear.Rb())*np.tan(self.phi_op()) \
                - np.sqrt(self.pinion.Ra()**2 - self.pinion.Rb()**2)))
        Rfg = self.gear.Rf()

    def undercut_check(self):
        pass

    def working_depth_check(self):
        # Check for adequate working depth
        print('***Working Depth Check***')
        print('\tWorking depth: {}'.format(self.pinion.a_op() + self.gear.a_op()))
        print('\tStandard working depth: {}'.format(2/self.pinion.Pd))

    def tip_thickness_check(self):
        # Check if tooth tips are sharp
        print('***Tooth Tip Thickness Check***')
        t_Ap = self.pinion.t_R(self.pinion.Ra())
        if t_Ap < 0:
            print('FAIL: Sharp pinion tooth tips')
        t_Ag = self.pinion.t_R(self.gear.Ra())
        if t_Ag < 0:
            print('FAIL: Sharp gear tooth tips')


class InternalSpurGearPair(GearPair):

    def Cs(self):
	# Standard center distance
        return self.gear.Rs() - self.pinion.Rs()
       
    def phi_op(self):
	# Operating pressure angle
        x = (self.gear.Rb() - self.pinion.Rb()) / self.C
        return np.arccos(x)

    def mc(self):
        # Profile contact ratio
        Z = self.pinion.rho_a() - self.gear.rho_a() - (self.pinion.rho_p() - self.gear.rho_p())
        mc = Z / self.pinion.pb()
        return mc


class Tests:

    def test1(self):
        p =  gears.ExternalSpurGear(N=16, Pd=1/8, phi=np.radians(20), ts=14.90, a=11.2)
        g =  gears.ExternalSpurGear(N=95, Pd=1/8, phi=np.radians(20), ts=16.76, a=13.8)
        pair = ExternalSpurGearPair(p, g, C=453)
        print(pair)

        
    def test2(self):
        p =  gears.ExternalSpurGear(N=16, Pd=1/8, phi=np.radians(20), ts=14.90, a=11.2)
        g =  gears.ExternalSpurGear(N=95, Pd=1/8, phi=np.radians(20), ts=16.76, a=13.8)
        pair = ExternalSpurGearPair(p, g, C=453)
        pair.animate(area='Mesh')


if __name__ == '__main__':
    t = Tests()
    t.test2()




