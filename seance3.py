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



point_1 = points_trim()[15]
point_2 = points_trim()[14]



def nonLinearModel(aircraft, point_trim, wh, time):

    alt, Ma, sm, km = point_trim
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



def LinearModel(aircraft, point_trim, wh, time):

    alt, Ma, sm, km = point_trim
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
time = np.arange(0, 250, 0.1)



XnonLin = nonLinearModel(aircraft, point_2, wh, time)
fignonLin = plt.figure()
dyn.plot(time, XnonLin, figure=fignonLin)



XLin = LinearModel(aircraft, point_1, wh, time)
figLin = plt.figure()
dyn.plot(time, XLin, figure=figLin)



figDiff = plt.figure()
dyn.plot(time, abs(XLin - XnonLin), figure=figDiff)



def modal_form(aircraft, point_trim):

    alt, Ma, sm, km = point_trim
    aircraft.set_mass_and_static_margin(km, sm)
    va = dyn.va_of_mach(Ma, alt)
    xtrim, utrim = dyn.trim(aircraft, {'va': va, 'gamm': 0, 'h': alt})
    xtrim_4 = xtrim[2:]
    A, B = utils.num_jacobian(xtrim, utrim, aircraft, dyn.dyn)
    A_4, B_4 = A[2:,2:], B[2:,:2]
    valeursp, M = np.linalg.eig(A_4)
    M_inv = np.linalg.inv(M)
    Am_4 = np.diag(valeursp)
    Bm_4 = np.dot(M_inv, B_4)
    return '{}\n{}'.format(Am_4, Bm_4)




def stability(aircraf, point_trim):

    alt, Ma, sm, km = point_trim
    aircraft.set_mass_and_static_margin(km, sm)
    va = dyn.va_of_mach(Ma, alt)
    xtrim, utrim = dyn.trim(aircraft, {'va': va, 'gamm': 0, 'h': alt})
    A, B = utils.num_jacobian(xtrim, utrim, aircraft, dyn.dyn)
    A_4 = A[2:,2:]
    valeursp, M = np.linalg.eig(A_4)
    V_reals = []
    for V in valeursp:
        V_reals.append(V.real)
    return V_reals



def controllability(aircraf, point_trim):

    alt, Ma, sm, km = point_trim
    aircraft.set_mass_and_static_margin(km, sm)
    va = dyn.va_of_mach(Ma, alt)
    xtrim, utrim = dyn.trim(aircraft, {'va': va, 'gamm': 0, 'h': alt})
    xtrim_4 = xtrim[2:]
    A, B = utils.num_jacobian(xtrim, utrim, aircraft, dyn.dyn)
    A_4, B_4 = A[2:, 2:], B[2:, :2]
    Q = np.zeros((4, 4*2))
    for i in range(3):
        Q[:,2*i:2*(i+1)] = np.dot(np.linalg.matrix_power(A_4, i),B_4)
    return Q



def function_transfere(aircraft, point_trim):

    alt, Ma, sm, km = point_trim
    aircraft.set_mass_and_static_margin(km, sm)
    va = dyn.va_of_mach(Ma, alt)
    xtrim, utrim = dyn.trim(aircraft, {'va': va, 'gamm': 0, 'h': alt})
    A, B = utils.num_jacobian(xtrim, utrim, aircraft, dyn.dyn)
    A_4, B_4 = A[2:,2:], B[2:,:2][:,0].reshape((4,1))
    C_4 = np.zeros((1,4))
    C_4[0][2] = 1
    Acc_4 = np.zeros((4, 4))
    Bcc_4 = np.zeros((4, 1))
    Bcc_4[3,0] = 1
    valeursp, M = np.linalg.eig(A_4)
    coeff = np.poly(valeursp)
    for i in range(3):
        Acc_4[i,i+1] = 1
        Acc_4[3,3-i] = -coeff[i+1]
    Acc_4[3,0] = -coeff[4]
    Qccc = np.zeros((4, 4))
    Q = np.zeros((4, 4))
    for i in range(4):
        Qccc[:,i:i+1] = np.dot(np.linalg.matrix_power(Acc_4, i), Bcc_4)
        Q[:,i:i+1] = np.dot(np.linalg.matrix_power(A_4, i), B_4)
    Mcc = np.dot(Q, np.linalg.inv(Qccc))
    Ccc_4 = np.dot(C_4, Mcc)
    num = list(Ccc_4)
    num.reverse()
    den = list(-Acc_4[3, :])
    den.append(1.)
    den.reverse()
    return num, den, valeursp



def pade_reduction(aircraft, point_trim):

    num, den, valeursp = function_transfere(aircraft, point_trim)
    poles = sorted(valeursp, key=lambda x: abs(x))
    pade_poles = poles[0:2]
    pade_den = np.poly(pade_poles)
    y_0, y_1 = pade_poles
    y_01 = y_0 + y_1
    pade_num = [0,0]
    a_0, c_0, a_1, c_1 = den[-1], num[0][-1], den[-2], num[0][-2]
    pade_num[1] = c_0/a_0*y_0*y_1
    pade_num[0] = 1/a_0**2*(c_1*a_0*y_0*y_1 - c_0*a_1*y_0*y_1 - c_0*a_0*y_01)
    return pade_num, pade_den



print(function_transfere(aircraft, point_1))

print(pade_reduction(aircraft, point_1))