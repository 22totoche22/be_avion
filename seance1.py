#-*- coding: utf-8 -*-
import dynamic, utils
import numpy as np
import matplotlib.pyplot as plt
import math

''''''''''''

avion = dynamic.Param_A321() #avion choisi
rho0 = utils.isa(0)[1] #rho0
nb_points = 100 #nombre de points pour les courbes
q=0 #angle de braquage de l'empennage horozontal en rad

''''''''''''

def get_propulsion_M(h,m,nb): #calcule la poussée maximale en focntion du nombre de mach, (altitude m, intervalle pour le mach, nombre de points)
    rho = utils.isa(h)[1]
    M = np.linspace(m[0], m[1], nb)
    F = np.array([avion.F0*math.pow(rho/rho0, 0.6)*(0.568+0.25*math.pow(1.2-mach, 3)) for mach in M])
    return F,M

interv_M = [0.4,0.9] #intervalle du mach
h = [3000, 10000] #altitude en m

prop = plt.figure()
for i in h:
    F,M = get_propulsion_M(i,interv_M,nb_points)
    plt.plot(M,F,label = "h="+str(i)+"m",figure=prop)
    plt.xlabel("Mach")
    plt.ylabel("Poussee maximale en N")
    plt.title("Poussee maximale en fonction \n du nombre de mach ")
    plt.legend()
plt.show()

''''''''''''''''''
va = 100 #vitesse de l'avion en m/s

def get_CL_alpha(va,q,dphr,P,alpha,nb): # calcule CL en fonction de l'angle d'incidence alpha , (vitesse m/s,vitesse de tangage rad/sec, angle de braquage de l'empennage horozontal en rad, intervalle d'incidnce alpha  rad, nombre de points)
    angle_alpha = np.linspace(alpha[0],alpha[1],nb)
    CL = np.array([dynamic.get_aero_coefs(va, alpha, q, dphr, P)[0] for alpha in angle_alpha])
    return utils.deg_of_rad(angle_alpha), CL #angle_alpha en degré

interv_alpha = [-10*math.pi/180,20*math.pi/180] #intervalle d'incidence alpha en rad
dphr = [-30*math.pi/180,+20*math.pi/180] #valeurs de dphr en rad

cl = plt.figure()
for i in dphr:
    alpha, CL = get_CL_alpha(va,q,i,avion,interv_alpha,nb_points)
    plt.plot(alpha,CL,label = "dphr="+str(utils.deg_of_rad(i))+" degre",figure=cl)
    plt.xlabel("Angle d'incidence alpha en degre")
    plt.ylabel("Coefficient de portance CL")
    plt.title("CL en fonction \n d'alpha ")
    plt.legend()
plt.show()

''''''''

def get_Cm_alpha(va,q,dphr,P,alpha,nb,ms): # calcule Cm en fonction de l'angle d'incidence alpha, (vitesse m/s, vitesse de tangage rad/sec, angle de braquage de l'empennage horizontal en rad, intervalle d'incidnce alpha  rad, nombre de points)
    angle_alpha = np.linspace(alpha[0], alpha[1], nb)
    Cm = np.array([ P.Cm0 - ms * P.CLa * (k - P.a0) + P.Cmq * P.lt / va * q + P.Cmd * dphr for k in angle_alpha])
    return utils.deg_of_rad(angle_alpha), Cm #angle alpha en degré

ms=[-0.1,0,0.2,1] #valeurs de marge statique
dphr0 = 0

cm = plt.figure()
for i in ms:
    alpha, Cm = get_Cm_alpha(va,q,dphr0,avion,interv_alpha,nb_points,i)
    plt.plot(alpha,Cm,label = "ms="+str(i)+"",figure=cm)
    plt.xlabel("Angle d'incidence alpha en degre")
    plt.ylabel("Coefficient du moment de tangage Cm")
    plt.title("Cm en fonction \n d'alpha ")
    plt.legend()
plt.show()

''''''''
qe =0 #vol stabilié
Cme =0 #  moment de tangage nul

def dphre_alpha(va,q,P,alpha,nb,ms): #dphre en fonction de de l'incidence alpha e ,(idem que precedemment)
    angle_alpha = np.linspace(alpha[0], alpha[1], nb)
    P.set_mass_and_static_margin(P.m_k, ms)
    dphre = np.array([-(P.Cm0 - ms * P.CLa * (k - P.a0) + P.Cmq * P.lt / va * q )/ P.Cmd  for k in angle_alpha])
    return utils.deg_of_rad(angle_alpha), utils.deg_of_rad(dphre) #angle alpha et dphre en degré

dphrealpha = plt.figure()
for i in ms:
    alpha,dphre = dphre_alpha(va, qe, avion, interv_alpha, nb_points, i)
    plt.plot(alpha,dphre,label = "ms="+str(i)+"",figure=dphrealpha)
    plt.xlabel("Angle d'incidence alpha en degre")
    plt.ylabel("angle de braquage de l'empennage horizontal dphre en degre")
    plt.title("dphre en fonction \n d'alpha pour vitesse de tangage q=0 et moment de tangage Cm=0")
    plt.legend()
plt.show()


coeff_vt = [0.9,1,1.1] #coefficient pour volume d'empennage Vt

def dphre_vt(va,q,P,alpha,nb,ms,coeff): #dphre en fonction de de l'incidence alpha e ,(idem que precedemment, coefficient pour volume d'empennage Vt)
    angle_alpha = np.linspace(alpha[0], alpha[1], nb)
    P.set_mass_and_static_margin(P.m_k, ms)
    dphre= np.array([-(P.Cm0 - ms * P.CLa * (k - P.a0) + P.Cmq * P.lt / va * q )/ (P.Cmd*coeff)  for k in angle_alpha])
    return utils.deg_of_rad(angle_alpha), utils.deg_of_rad(dphre) #angle alpha et dphre en degré

dphrevt = plt.figure()
ms_vt = ms[2]
for i in coeff_vt:
    alpha, dphre = dphre_vt(va, qe, avion, interv_alpha, nb_points, ms_vt,i)
    plt.plot(alpha,dphre,label = "vt="+str(-avion.Cmd*i/avion.CLat)+"",figure=dphrevt)
    plt.xlabel("Angle d'incidence alpha en degre")
    plt.ylabel("angle de braquage de l'empennage horizontal dphre en degre")
    plt.title("dphre en fonction \n d'alpha pour ms ="+str(ms_vt)+"")
    plt.legend()
plt.show()

''''''''''''

def CLE_alpha(va,q,P,alpha,nb,ms): #coefficident de portance equilibrée Cle lorsque dphr=dphre, (idem)
    angle_alpha = np.linspace(alpha[0],alpha[1],nb)
    dphr = np.array([-(P.Cm0 - ms * P.CLa * (k - P.a0) + P.Cmq * P.lt / va * q) / P.Cmd for k in angle_alpha])
    CLE = np.array([dynamic.get_aero_coefs(va, alpha, q, dphr[i], P)[0] for i,alpha in enumerate(angle_alpha)])
    return utils.deg_of_rad(angle_alpha), CLE

clealpha = plt.figure()
for i in range(2,len(ms)):
    alpha, CLE = CLE_alpha(va,qe,avion,interv_alpha,nb_points,ms[i])
    plt.plot(alpha,CLE,label = "ms="+str(ms[i])+"",figure=clealpha)
    plt.xlabel("Angle d'incidence alpha en degre")
    plt.ylabel("Coefficient de portance equilibree CLE")
    plt.title("CLE en fonction \n d'alpha ")
    plt.legend()
plt.show()

''''''''''''

def polaire(avion,ms,qe,va,alpha,nb): #polaire, idem
    CLE = CLE_alpha(va,qe,avion,alpha,nb,ms)[1]
    CDE = np.array([avion.CD0 + avion.ki * CL ** 2 for CL in CLE])
    return CLE,CDE

polair=plt.figure()
for i in range(2,len(ms)):
    CLE,CDE = polaire(avion,ms[i],qe,va,interv_alpha,nb_points)
    plt.plot(CDE,CLE,label = "ms="+str(ms[i])+"",figure=polair)
    plt.xlabel("CDe")
    plt.ylabel("CLe")
    plt.title("Polaire equilibree")
    plt.legend()
plt.show()

fmax = max(CLE / CDE)
print ("finesse_max :"+str(fmax))