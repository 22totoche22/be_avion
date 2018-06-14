#-*- coding: utf-8 -*-
import dynamic, utils
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint

import numpy as np
import matplotlib.pyplot as plt
import dynamic as dyn
import utils as ut
import scipy.integrate

import pdb


def points_trim():
    h = [3000, 10000]  # altitude en m
    Ma = [0.4, 0.9]  # intervalel nombre de mach
    ms = [-0.2, 0.5]  # marge statique
    km = [0.1, 1]  # coefficient de réglage de la masse
    points = []
    for hi in h:
        for Mai in Ma:
            for msi in ms:
                for kmi in km:
                    points.append([hi, Mai, msi, kmi])
    return points




point = points_trim()[15]
alt, Ma, sm, km = point


def nonLinearModel(aircraft, alt, Ma, sm, km, wh, time):
    aircraft.set_mass_and_static_margin(km, sm)
    va = dyn.va_of_mach(Ma, alt)
    xtrim, uTrim = dyn.trim(aircraft, {'va':va, 'h':alt})
    xtrim[dyn.s_a] = xtrim[dyn.s_a] + np.arctan(wh/xtrim[dyn.s_va])
    X = scipy.integrate.odeint(dyn.dyn, xtrim, time, args=(uTrim, aircraft))
    return X



''' Modèle linéaire '''


def dynlin(dX, t, A, B, dU):
    dxdot = np.dot(A,dX) + np.dot(B,dU)
    return dxdot





def LinearModel(aircraft, alt, Ma, sm, km, wh, time):

    aircraft.set_mass_and_static_margin(km, sm)
    va = dyn.va_of_mach(Ma, alt)
    xtrim, utrim = dyn.trim(aircraft, {'va':va,'gamm':0, 'h':alt})
    xtrim[dyn.s_a] = xtrim[dyn.s_a] + np.arctan(wh / xtrim[dyn.s_va])
    wing = [0] * 6
    A, B = utils.num_jacobian(xtrim, utrim, aircraft, dyn.dyn)
    du = np.zeros((4,))
    wing[dyn.s_a] += np.arctan(wh / xtrim[dyn.s_va])
    dX = scipy.integrate.odeint(dynlin, wing, time, args=(A, B, du))
    X = dX+xtrim
    for i in range(len(X)):
        X[i][dyn.s_y] = X[i][dyn.s_va]*time[i]
    return X


''' Affichages des trajectoires '''

aircraft = dyn.Param_A321()
wh = 2 #m/s
time = np.arange(0, 240, 0.1)
delta_t = (time[1]-time[0])/len(time)



XnonLin = nonLinearModel(aircraft, alt, Ma, sm, km, wh, time)
fignonLin = plt.figure()
dyn.plot(time, XnonLin, figure=fignonLin)




XLin = LinearModel(aircraft, alt, Ma, sm, km, wh, time)
figLin = plt.figure()
dyn.plot(time, XLin, figure=figLin)



plt.show()


figDiff = plt.figure()
dyn.plot(time, abs(XLin - XnonLin), figure=figDiff)


plt.show()



print ("le point de trim choisis est : \n\n"
       "h = "+str(alt)+"\n"
       "Ma = "+str(Ma)+"\n"
       "ms = "+str(sm)+"\n"
       "km = "+str(km))










