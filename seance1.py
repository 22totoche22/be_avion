#-*- coding: utf-8 -*-
import dynamic, utils
import numpy as np
import matplotlib.pyplot as plt
import math


avion = dynamic.Param_A321()
rho0 = utils.isa(0)[1]


''''''''''''

def get_propulsion_M(h,m,nb): # altitude h(m), m is an interval for the mach, nb number of points
    rho = utils.isa(h)[1]
    M = np.linspace(m[0], m[1], nb)
    F = np.array([avion.F0*math.pow(rho/rho0, 0.6)*(0.568+0.25*math.pow(1.2-mach, 3)) for mach in M])
    return F,M


interv_M = [0.4,0.9]
h1 = 3000
h2 = 10000

F1,M = get_propulsion_M(h1,interv_M,100)
F2,M = get_propulsion_M(h2,interv_M,100)

plt.plot(M,F1,label = "h=3000m")
plt.plot(M,F2,label = "h=10000m")

plt.xlabel("Mach")
plt.ylabel("Poussee maximale en N")
plt.title("Poussee maximale en fonction \n du nombre de mach ")
plt.legend()
plt.show()




''''''''''''''''''
va = 100

def get_CL_alpha(va,q,dphr,P,alpha,nb):
    angle_alpha = np.linspace(alpha[0],alpha[1],nb)
    CL = np.array([dynamic.get_aero_coefs(va, alpha, q, dphr, P)[0] for alpha in angle_alpha])
    return angle_alpha, CL


interv_alpha = [-10*math.pi/180,20*math.pi/180]

alpha, CL = get_CL_alpha(va,0,-30*math.pi/180,avion,interv_alpha,100)
alpha, CL1 = get_CL_alpha(va,0,+20*math.pi/180,avion,interv_alpha,100)

plt.plot(alpha,CL,label = "dphr=-30*math.pi/180 rad")
plt.plot(alpha,CL1,label = "dphr=+20*math.pi/180 rad")

plt.xlabel("Alpha en radian")
plt.ylabel("CL")
plt.title("CL en fonction \n d'alpha ")
plt.legend()
plt.show()

''''''''

def get_Cm_alpha(va,q,dphr,P,alpha,nb,ms):
    angle_alpha = np.linspace(alpha[0], alpha[1], nb)
    Cm = np.array([ P.Cm0 - ms * P.CLa * (k - P.a0) + P.Cmq * P.lt / va * q + P.Cmd * dphr for k in angle_alpha])
    return angle_alpha, Cm

ms=[-0.1,0,0.2,1]


for i in range(4):
    alpha, vars()['{0}{1}'.format('Cm', i)] = get_Cm_alpha(va,0,0,avion,interv_alpha,100,ms[i])

plt.plot(alpha,Cm0,label = "ms=-0.1")
plt.plot(alpha,Cm1,label = "ms=0")
plt.plot(alpha,Cm2,label = "ms=0.2")
plt.plot(alpha,Cm3,label = "ms=1")

plt.xlabel("Alpha en radian")
plt.ylabel("Cm")
plt.title("Cm en fonction \n d'alpha ")
plt.legend()
plt.show()

''''''''

def dpha_alpha(va,q,P,alpha,nb,ms):
    angle_alpha = np.linspace(alpha[0], alpha[1], nb)
    P.set_mass_and_static_margin(P.m_k, ms)
    dphr = np.array([-(P.Cm0 - ms * P.CLa * (k - P.a0) + P.Cmq * P.lt / va * q )/ P.Cmd  for k in angle_alpha])
    return angle_alpha, dphr

def dpha_vt(va,q,P,alpha,nb,ms,coeff):
    angle_alpha = np.linspace(alpha[0], alpha[1], nb)
    P.set_mass_and_static_margin(P.m_k, ms)
    dphr = np.array([-(P.Cm0 - ms * P.CLa * (k - P.a0) + P.Cmq * P.lt / va * q )/ (P.Cmd*coeff)  for k in angle_alpha])
    return angle_alpha, dphr

coeff = [0.9,1,1.1]
for i in range(len(coeff)):
    alpha, vars()['{0}{1}'.format('vt', i)] = dpha_vt(va, 0, avion, interv_alpha, 100, ms[0],coeff[i])


for i in range(4):
    alpha, vars()['{0}{1}'.format('dphr', i)] = dpha_alpha(va,0,avion,interv_alpha,100,ms[i])

plt.plot(alpha,dphr0,label = "ms=-0.1")
plt.plot(alpha,dphr1,label = "ms=0")
plt.plot(alpha,dphr2,label = "ms=0.2")
plt.plot(alpha,dphr3,label = "ms=1")

plt.xlabel("Alpha en radian")
plt.ylabel("dphr")
plt.title("dphr en fonction \n d'alpha ")
plt.legend()
plt.show()

plt.plot(alpha,vt0,label = "vt=-0.1")
plt.plot(alpha,vt1,label = "vt=0")
plt.plot(alpha,vt2,label = "vt=0.2")

plt.xlabel("Alpha en radian")
plt.ylabel("dphr")
plt.title("dphr en fonction \n d'alpha ")
plt.legend()
plt.show()



#je sais pas pour vt, proportionnel
''''''''''''

def CLE_alpha(va,q,P,alpha,nb,ms):
    angle_alpha = np.linspace(alpha[0],alpha[1],nb)
    dphr = np.array([-(P.Cm0 - ms * P.CLa * (k - P.a0) + P.Cmq * P.lt / va * q) / P.Cmd for k in angle_alpha])
    CLE = np.array([dynamic.get_aero_coefs(va, alpha, q, dphr[i], P)[0] for i,alpha in enumerate(angle_alpha)])
    return angle_alpha, CLE

alpha, CLE0 = CLE_alpha(va,0,avion,interv_alpha,100,ms[2])
alpha, CLE1 = CLE_alpha(va,0,avion,interv_alpha,100,ms[3])

plt.plot(alpha,CLE0,label = "ms=0.2")
plt.plot(alpha,CLE1,label = "ms=1")

plt.xlabel("Alpha en radian")
plt.ylabel("CLE")
plt.title("CLE en fonction \n d'alpha ")
plt.legend()
plt.show()

''''''''''''


def polaire(avion,ms,q=0,va=100,alpha=interv_alpha,):

    CLE = CLE_alpha(va,q,avion,alpha,100,ms)[1]
    CD = np.array([avion.CD0 + avion.ki * CL ** 2 for CL in CLE])
    return CLE,CD

CD,CLE0 = polaire(avion,ms[2])
CD,CLE1 = polaire(avion,ms[3])

plt.plot(CLE0,CD,label = "ms=0.2")
plt.plot(CLE1,CD,label = "ms=1")

plt.xlabel("CD")
plt.ylabel("CL")
plt.title("polaire")
plt.legend()
plt.show()