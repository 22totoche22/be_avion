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
    h = [3000.0, 10000.0]  # altitude en m
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



choix = np.random.random_integers(0, len(points_trim())-1)
point = points_trim()[choix]
alt, Ma, sm, km = point


def nonLinearModel(aircraft, alt, Ma, sm, km, wh, time):
    aircraft.set_mass_and_static_margin(km, sm)
    va = dyn.va_of_mach(Ma, alt)
    xtrim, uTrim = dyn.trim(aircraft, {'va':va, 'h':alt})
    xtrim[dyn.s_a] = xtrim[dyn.s_a] + np.arctan(wh/xtrim[dyn.s_va])
    X = scipy.integrate.odeint(dyn.dyn, xtrim, time, args=(uTrim, aircraft))
    return X




aircraft = dyn.Param_A321()
wh = 2 #m/s

time = np.arange(0., 240, 0.1)
XnonLin = nonLinearModel(aircraft, alt, Ma, sm, km, wh, time)
fig = plt.figure()
dyn.plot(time, XnonLin, figure = fig)

plt.show()

''' Modèle linéaire '''


def dynlin(X, t, A, B, U):
    xdot = np.dot(A,X) + np.dot(B,U)
    return xdot





def LinearModel(aircraft, alt, Ma, sm, km, wh, time):
    aircraft.set_mass_and_static_margin(km, sm)
    va = dyn.va_of_mach(Ma, alt)
    xtrim, utrim = dyn.trim(aircraft, {'va':va, 'h':alt})
    xtrim[dyn.s_a] = xtrim[dyn.s_a] + np.arctan(wh/xtrim[dyn.s_va])
    A, B = utils.num_jacobian(xtrim, utrim, aircraft, dyn.dyn)
    X = scipy.integrate.odeint(dynlin, xtrim, time, args=(A, B, utrim))
    return X


XLin = LinearModel(aircraft, alt, Ma, sm, km, wh, time)
fig = plt.figure()
dyn.plot(time, XLin, figure = fig)
plt.show()



print ("le point de trim choisis est : \n\n"
       "h = "+str(alt)+"\n"
       "Ma = "+str(Ma)+"\n"
       "ms = "+str(sm)+"\n"
       "km = "+str(km))










