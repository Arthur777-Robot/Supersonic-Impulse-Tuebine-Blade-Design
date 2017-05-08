# -*- coding: utf-8 -*-
'''
The MIT License (MIT)
Copyright (c) 2017 Seiji Arther Murakami

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

This program was based on below papers
1 - NACA RM L52B06, APPLICATION OF SUPERSONIC VORTEX-FLOW THEORY TO THE DESIGN OF SUPERSONIC IMPULSE COMPRESSOROR TURBINE-BLADE SECTIONS
2 - DESIGN OF TURBINE BLADES SUITABLE FOR SUPERSONIC RELATIVE INLET VELOCITIES AND THE INVESTIGATION OF THEIR PERFORMANCE IN CASCADES:	PART I-THEORY AND DESIGN
3 - DESIGN OF TURBINE BLADES SUITABLE FOR SUPERSONIC RELATIVE INLET VELOCITIES AND THE INVESTIGATION OF THEIR PERFORMANCE IN CASCADES:	PART 11-EXPERIMENTS, RESULTS AND DISCUSSION
4 - NASA TN D-4421, ANALYTICAL INVESTIGATION OF SUPERSONIC TURBOMACHINERY BLADING: I - Computer Program for Blading Design
5 - NASA TN D-4422, ANALYTICAL INVESTIGATION OF SUPERSONIC TURBOMACHINERY BLADING: II - Analysis of Impulse Turbine-Blade Sections
'''

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy import integrate
from scipy import optimize

plt.close('all')


# Characteristic line function

class Blade():
    def __init__(self, gamma, mach):
        print('Initializing parameter...')
        self.gamma = gamma
        self.Rstar_min = math.sqrt((self.gamma - 1)/(self.gamma + 1))
        self.const = self.chara_line(1)
        self.ve = int(round(self.get_Pr(mach)))
        print("Inlet Mach = ", mach, "Inlet Prandtle meyer angle = ", self.ve)

    # if Rstar value is less than self.Rstar_min, if will give math error
    def chara_line(self, Rstar):

        fai = 0.5 * (np.sqrt((self.gamma + 1)/(self.gamma - 1)) * np.arcsin((self.gamma - 1) / Rstar**2 - self.gamma) + np.arcsin((self.gamma + 1) * Rstar**2 - self.gamma))

        return fai

    def chara_x(self, Rstar, theta):
        return Rstar * np.sin(theta)

    def chara_y(self, Rstar, theta):
        return Rstar * np.cos(theta)

    # draws 2 characteristic lines. the angle starts from v1(compression wave) and v2(expansion wave).
    def draw_lines(self, v1, v2):

        v1 = np.deg2rad(v1)
        v2 = np.deg2rad(v2)

        counter = 1000
        i = np.arange(counter)
        Rstar = self.Rstar_min + i / counter
        Rstar = Rstar[Rstar < 1]
        theta1 = (self.chara_line(Rstar) - self.const) + v1
        theta2 = -(self.chara_line(Rstar) - self.const) + v2
        x0 = self.chara_x(Rstar, theta1)
        y0 = self.chara_y(Rstar, theta1)
        x1 = self.chara_x(Rstar, theta2)
        y1 = self.chara_y(Rstar, theta2)

        num = 360
        theta = np.arange(0, num)
        R_x = - np.sin(np.deg2rad(theta))
        R_y = np.cos(np.deg2rad(theta))
        R_x_min = - np.sin(np.deg2rad(theta)) * self.Rstar_min
        R_y_min = np.cos(np.deg2rad(theta)) * self.Rstar_min

        plt.plot(x0, y0)
        plt.plot(x1, y1)
        plt.plot(R_x, R_y)
        plt.plot(R_x_min, R_y_min)
        plt.xlim(-1, 1)
        plt.ylim(-1, 1)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    def plot_chara(self, angle, step):

        counter = 1000
        i = np.arange(counter)
        Rstar = self.Rstar_min + i / counter
        Rstar = Rstar[Rstar < 1]
        fai = self.chara_line(Rstar)
        for j in range(0, angle, step):
            x1 = self.chara_x(Rstar, fai - self.const + np.deg2rad(j))
            y1 = self.chara_y(Rstar, fai - self.const + np.deg2rad(j))
            x2 = self.chara_x(Rstar, - (fai - self.const - np.deg2rad(j)))
            y2 = self.chara_y(Rstar, - (fai - self.const - np.deg2rad(j)))

            plt.plot(x1, y1, "r")
            plt.plot(x2, y2, "k")
        plt.xlim(-1, 1)
        plt.ylim(-1, 1)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    # top arc angle is 0
    # v1 must be smaller than v2
    def get_R(self, v1, v2):

        v1 = math.radians(v1)
        v2 = math.radians(v2)
        Rstar = 1.0	 # initial Rstar
        param = []
        temp_x = 1
        temp_y = 1
        myu_check = math.radians(90)

        # this for is very stupid. intersecting angle is always half the delta
        for num in range(1, 10):  # decimal precision of Rstar
            while(1):
                theta1 = (self.chara_line(Rstar) - self.const) + v1
                theta2 = -(self.chara_line(Rstar) - self.const) + v2
                x0 = self.chara_x(Rstar, theta1)
                y0 = self.chara_y(Rstar, theta1)
                x1 = self.chara_x(Rstar, theta2)
                y1 = self.chara_y(Rstar, theta2)
                # print(Rstar, math.degrees(theta1), math.degrees(theta2), x0, x1)
                # print("R* = %.5f,\ttheta1 = %.5f,\ttheta2 = %.5f,\tx0 = %.5f,\tx1 = %.5f" % (Rstar, math.degrees(theta1), math.degrees(theta2), x0, x1))

                if (x0 == x1) and (y0 == y1):
                    myu_check = 0
                    break
                elif x1 - x0 < 0:
                    myu_check = np.arctan2(temp_y - y1, temp_x - x1)
                    # if myu_check == 0:
                    #     print(theta1, theta2)
                    temp_x = x1
                    temp_y = y1
                    Rstar += 1.0/(10**num)

                    break
                else:
                    Rstar -= 1.0/(10**num)

                    if Rstar < self.Rstar_min:
                        Rstar += 1.0/(10**num)
                        break

            if (x0 == x1) and (y0 == y1):
                break

        param.append(Rstar)
        param.append(x0)
        param.append(y0)
        param.append(myu_check)
        return param

    def get_myu(self, M):
        return np.arcsin(1/M)

    def get_Mstar(self, Rstar):
        return 1/Rstar

    def get_mach(self, Mstar):
        return np.sqrt((2 * Mstar**2)/((self.gamma + 1) - (self.gamma - 1) * Mstar**2))

    def get_Ru(self, vu):
        return self.get_R(-vu, vu)[0]

    # In the papaer, theta are defines around 90-120deg
    def get_upper_arc(self, Rstar, theta_in, theta_out, vu, vout):

        alpha_in = 90 - theta_in - (vu - self.ve)
        alpha_out = 90 - theta_out - (vu - vout)

        print("Upper arc alpha_in = %.1f deg, alpha_out = %.1f deg" % (alpha_in, alpha_out))

        theta = np.arange(-alpha_in, alpha_out, 0.1)
        x = Rstar * np.sin(np.deg2rad(theta))
        y = Rstar * np.cos(np.deg2rad(theta))

        # print(beta_in,beta_out)
        # print(x[0],y[0])
        # print(x[-1],y[-1])
        plt.plot(x, y)
        plt.xlim(-3, 3)
        plt.ylim(-3, 3)
        plt.gca().set_aspect('equal', adjustable='box')
        # plt.show()

    def get_lower_arc(self, Rstar, theta_in, theta_out, vl, vout, shift=0):

        alpha_in = 90 - theta_in - (self.ve - vl)
        alpha_out = 90 - theta_out - (vout - vl)

        print("Lower arc alpha_in = %.1f deg, alpha_out = %.1f deg" % (alpha_in, alpha_out))

        theta = np.arange(-alpha_in, alpha_out, 0.1)
        x = Rstar * np.sin(np.deg2rad(theta))
        y = shift + Rstar * np.cos(np.deg2rad(theta))

        # print(alpha_in,alpha_out)
        # print(x[0],y[0])
        # print(x[-1],y[-1])

        plt.plot(x, y)
        plt.xlim(-3, 3)
        plt.ylim(-3, 3)
        plt.gca().set_aspect('equal', adjustable='box')
        # plt.show()

    # defines upper convex arc coordinates
    # ve is the entry angle
    # vu is the value to define the upper convex radius
    def upper_convex(self, vu,vend,theta):

        x, y = [], []
        Xstar, Ystar = [], []

        xtmp = 0
        ytmp = self.get_Ru(vu)

        for num in range(0, int(round(vend*2))):
            Xstar_b = xtmp
            Ystar_b = ytmp
            Rstar, Xstar_a, Ystar_a, myu_check = self.get_R(-vu, vu-num)
            myu = self.get_myu(self.get_mach(self.get_Mstar(Rstar)))
            a1 = math.tan(myu+math.radians(num/2.0))
            # a1 = math.tan(myu_check)
            b1 = Ystar_a - a1 * Xstar_a
            a2 = math.tan(math.radians(num/2.0))
            # a2 = math.tan(myu_check-myu)
            b2 = Ystar_b - a2 * Xstar_b

            # print(Rstar,Xstar_a,Ystar_a,math.degrees(myu_check),math.degrees(myu))

            xtmp = ((b2 - b1) / (a1 - a2))
            ytmp = xtmp * a2 + b2

            # print(num/2,Rstar,math.degrees(myu_check),math.degrees(myu),xtmp,ytmp)

            rotx, roty = self.rotate(xtmp, ytmp, theta)

            x += [(rotx)]
            y += [(roty)]
            Xstar += [(Xstar_a)]
            Ystar += [(Ystar_a)]

            xcomp = np.arange(-1, 0, 0.001)
            ycomp = a1*xcomp + b1

            # plt.plot(xcomp,ycomp)

        theta90 = np.arange(0, 90)
        R_x = - np.sin(np.deg2rad(theta90))
        R_y = np.cos(np.deg2rad(theta90))
        R_x_min = - np.sin(np.deg2rad(theta90)) * self.get_Ru(vu)
        R_y_min = np.cos(np.deg2rad(theta90)) * self.get_Ru(vu)

        plt.plot(x, y)
        # plt.plot(Xstar, Ystar)
        # plt.plot(R_x, R_y)
        # plt.plot(R_x_min, R_y_min)
        plt.xlim(-3, 3)
        plt.ylim(-3, 3)
        plt.gca().set_aspect('equal', adjustable='box')
        # plt.show()

        return x[-1], y[-1]

    # defines lower concave arc coordinates
    # Precaution. there are still calculation uncertainties.
    def lower_concave(self, v1, theta, shift=0):

        xinit = -math.sin(math.radians(v1))
        yinit = math.cos(math.radians(v1))

        xtmp = xinit
        ytmp = yinit

        x, y = [], []
        Xstar, Ystar = [], []

        for num in range(v1, self.ve*2):

            Xstar_b = xtmp
            Ystar_b = ytmp

            Rstar, Xstar_a, Ystar_a, myu_check = self.get_R(-num, v1)
            myu = self.get_myu(self.get_mach(self.get_Mstar(Rstar)))
            a1 = math.tan(-myu+math.radians(num/2.0))
            b1 = Ystar_a - a1 * Xstar_a
            a2 = math.tan(math.radians(num/2.0))
            b2 = Ystar_b - a2 * Xstar_b

            xtmp = ((b2 - b1) / (a1 - a2))
            ytmp = xtmp * a2 + b2

            # print(num/2.0,Rstar,math.degrees(myu_check),math.degrees(myu),xtmp,ytmp,Xstar_a,Ystar_a)

            rotx, roty = self.rotate(xtmp, ytmp, theta)

            x += [(rotx)]
            y += [(roty + shift)]
            Xstar += [(Xstar_a)]
            Ystar += [(Ystar_a)]

            xcomp = np.arange(-1, 0, 0.001)
            ycomp = a1*xcomp + b1

            # plt.plot(xcomp,ycomp)

        theta90 = np.arange(0, 90)
        R_x = - np.sin(np.deg2rad(theta90))
        R_y = np.cos(np.deg2rad(theta90))

        plt.plot(x, y)
        # plt.plot(Xstar, Ystar)
        # plt.plot(R_x, R_y)
        plt.xlim(-3, 3)
        plt.ylim(-3, 3)
        plt.gca().set_aspect('equal', adjustable='box')
        # plt.show()

        return x[-1], y[-1]

    def get_Q(self, Rlstar, Rustar):

        Mlstar = 1 / Rlstar
        Mustar = 1 / Rustar

        Mach = lambda Mach: (((self.gamma+1) / 2 - (self.gamma-1) * Mach**2 / 2)**(1 / (self.gamma - 1)))/Mach

        Q1 = Mlstar * Mustar / (Mustar - Mlstar)
        Q2 = integrate.quad(Mach, Mlstar, Mustar)

        Q = Q1 * Q2[0]

        return Q

    def get_Gstar(self, vl, vu, theta_in):

        Ae = 0
        Rlstar = self.get_Ru(vl)
        Rustar = self.get_Ru(vu)
        Q = self.get_Q(Rlstar, Rustar)

        Astar = (Rlstar - Rustar) / Q

        print(Astar)

        # Gstar = Ae/Astar * Q * (Rl - Ru) / math.cos(theta_in)
        print("daradara")

    def rotate(self, x, y, angle):

        theta = np.deg2rad(angle)
        a = np.array(((np.cos(theta), -np.sin(theta)),
                      (np.sin(theta), np.cos(theta))))
        b = np.array((x, y))
        return np.dot(a, b)

    def straight_line(self, theta, x, y, targetx):

        targety = np.tan(np.deg2rad(theta)) * targetx + y - np.tan(np.deg2rad(theta)) * x

        return targety

    def get_Pr(self, Mach):

        Mstar = (((self.gamma + 1) / 2 * Mach**2) / (1 + (self.gamma - 1) / 2 * Mach**2))**0.5
        tmp1 = math.pi/4 * (math.sqrt((self.gamma + 1)/(self.gamma - 1)) - 1)
        tmp2 = self.chara_line(1/Mstar)

        Pr = math.degrees(tmp1 + tmp2)

        return Pr

    def get_mach_from_prandtle_meyer(self, v1):
        mach0 = 1.0

        def func(mach, v1):
            return self.get_Pr(mach) - v1

        sol = optimize.root(func, mach0, args=(v1))
        mach = sol.x[0]
        return mach

    def valuables_limit(self, vl, vu):

        vlmin = 0
        vlmax = self.ve
        vumin = self.ve
        vumax = math.degrees((math.pi/2) * (math.sqrt((self.gamma + 1)/(self.gamma - 1)) - 1))

        print("upper max = %.2f, upper min = %.2f, lower max = %.2f, lower min = %.2f" % (vumax, vumin, vlmax, vlmin))

if __name__ == "__main__":

    gamma = 1.4
    mach_in = 2.5
    vout = 30
    vl = 0
    vu = 59

    total_turn_ang = 130
    theta_in = 90 - total_turn_ang / 2
    theta_out = 90 - total_turn_ang / 2

    print("Design Supersonic Turbine")

    f = Blade(gamma, mach_in)

    print("upper mach number = %.3f" % (f.get_mach_from_prandtle_meyer(vu)), "vu = ", vu)
    print("lower mach number = %.3f" % (f.get_mach_from_prandtle_meyer(vl)), "vl = ", vl)

    f.valuables_limit(vl, vu)

    # f.get_upper_arc(f.get_Ru(vu),theta_in_upper,theta_out_upper,vu,vout)
    # f.get_lower_arc(f.get_Ru(vl),theta_in_lower,theta_out_lower,vl,vout)
    # plt.show()
    #
    # a = f.rotate(1,0,45)
    # f.get_Gstar(vl,vu,theta_in_upper)

    Ualpha_in = 90 - theta_in - (vu - f.ve)
    Lalpha_in = 90 - theta_in - (f.ve - vl)
    
    print(Ualpha_in, Lalpha_in,Ualpha_in - theta_in)

    xlow, ylow = f.lower_concave(vl, Lalpha_in)
    xup, yup = f.upper_convex(vu,Ualpha_in - theta_in, Ualpha_in)

    f.get_lower_arc(f.get_Ru(vl), theta_in, theta_out, vl, vout)
    f.get_upper_arc(f.get_Ru(vu), theta_in, theta_out, vu, vout)

    newy = f.straight_line(total_turn_ang/2, xup, yup, xlow)
    plt.plot([xlow, xup], [newy, yup])
    shift = -abs(ylow - newy)

    f.lower_concave(vl, Lalpha_in, shift)
    f.get_lower_arc(f.get_Ru(vl), theta_in, theta_out, vl, vout, shift)
    plt.show()

    # f.draw_lines(-10, 20)
    # f.lower_concave(0, 30)
    # f.upper_convex(40, 35)
    # f.plot_chara(360, 5)
    print("finish")
