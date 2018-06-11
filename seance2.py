#-*- coding: utf-8 -*-
import dynamic, utils
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint

h = [3000.0,10000.0] #altitude en m
Ma = [0.4,0.9] #intervalel nombre de mach
ms = [-0.2,0.5] #marge statique
km = [0.1,1] #coefficient de réglage de la masse
nb = 100 #nb de points pour les courbes

avion = dynamic.Param_A321()
rho0 = utils.isa(0)[1]


''''''''''''

#calcul des points de Trim

alphatrim=np.zeros((len(h),len(Ma), len(ms),len(km)))
dphrtrim= np.zeros((len(h),len(Ma), len(ms),len(km)))
dthtrim = np.zeros((len(h),len(Ma), len(ms),len(km)))
prop =  np.zeros((len(h),len(Ma), len(ms),len(km)))

for ib,i in enumerate(ms):
    for jb,j in enumerate(km):
        avion.set_mass_and_static_margin(j, i)
        for kb,k in enumerate(Ma):
            for lb,l in enumerate(h):
                va_h0= dynamic.va_of_mach(k, l, k=1.4, Rs=287.05)
                X,U = dynamic.trim(avion, {'va':va_h0,'gamm': 0.,'h': l})
                alphatrim[ib,jb,kb,lb] = X[3]
                dphrtrim[ib,jb,kb,lb] = U[0]
                dthtrim[ib,jb,kb,lb] = U[1]
                prop[ib,jb,kb,lb] = dynamic.propulsion_model(X, U, avion)

def affiche():
    for i in range(len(ms)):
        for j in range(len(km)):
            for k in range(len(Ma)):
                        plt.suptitle('ms=' + str(ms[i]) + ', km=' + str((km[j])))
                        plt.subplot(1,4,1)
                        plt.title("alpha")
                        plt.plot(h,utils.deg_of_rad(alphatrim[i,j,k,:]),label = "Mach = "+str(Ma[k])+"")
                        plt.ylabel('alpha en degre')
                        plt.xlabel('altitude')
                        plt.legend()
                        plt.subplot(1, 4, 2)
                        plt.title("dphr")
                        plt.plot(h, utils.deg_of_rad(dphrtrim[i, j, k, :]), label = "Mach = "+str(Ma[k])+"")
                        plt.ylabel('dphr en degre')
                        plt.xlabel('altitude')
                        plt.legend()
                        plt.subplot(1, 4, 3)
                        plt.title("dth")
                        plt.plot(h, utils.deg_of_rad(dthtrim[i, j, k, :]), label = "Mach = "+str(Ma[k])+"")
                        plt.ylabel('dth en degre')
                        plt.xlabel('altitude')
                        plt.legend()
            for l in range(len(h)):
                        plt.subplot(1, 4, 4)
                        plt.title("poussee")
                        plt.plot(Ma, prop[i, j,: ,l ], label="altitude = " + str(h[l]) + "")
                        plt.ylabel('poussee')
                        plt.xlabel('Mach')
                        plt.legend()
            plt.show()

affiche()

''''''''''''''''''
#
valeur_trim = alphatrim[0,0,0,0], dphrtrim[0,0,0,0], dthtrim[0,0,0,0]
point_trim = [3000.0, 0.4, -0.2, 0.1]
#
h_trim, Ma_trim, ms_trim, km_trim = point_trim
# p_trim, rho_trim, T_trim = utils.isa(h_trim)
#
#
# avion.set_mass_and_static_margin(km_trim, ms_trim)
# va_trim = dynamic.va_of_mach(Ma_trim, h_trim, k = 1.4, Rs = 287.05)
# Cl = 2*avion.m*avion.g/(rho_trim*avion.S*va_trim**2)
#
# interv_alpha = [-10*math.pi/180,20*math.pi/180]
#
# def CLE_alpha(va,q,P,alpha,nb,ms):
#     angle_alpha = np.linspace(interv_alpha[0],interv_alpha[1],nb)
#     dphr = np.array([-(P.Cm0 - ms * P.CLa * (k - P.a0) + P.Cmq * P.lt / va * q) / P.Cmd for k in angle_alpha])
#     CLE = np.array([dynamic.get_aero_coefs(va, alpha, q, dphr[i], P)[0] for i,alpha in enumerate(angle_alpha)])
#     return angle_alpha, CLE
#
# alpha, CLE = CLE_alpha(va_trim,0,avion,interv_alpha,100,ms_trim)
#
#
# alphatrim_0 = valeur_trim[0]
# alphatrim_0 = float(int(alphatrim_0*1000 + 0.5))/1000
# Cl = float(int(Cl*1000 +0.5))/1000
#
# plt.plot(alpha,CLE,label = "ms = "+str(ms_trim))
# plt.xlabel("Alpha (rad)")
# plt.ylabel("CLE")
# plt.title("CLE en fonction \n d'alpha ")
# ligne_cl = [Cl]*nb
# idx = np.argwhere(np.diff(np.sign(CLE - ligne_cl)) != 0).reshape(-1)
# idx = idx[0]
# x = alpha[:idx]
# alpha[idx] = float(int(alpha[idx]*1000 + 0.5))/1000
# colonne_alpha = [alpha[idx]]*2
# y = [CLE[0], Cl]
# colonne_alphatrim = [alphatrim_0]*2
# plt.plot(x, ligne_cl[:idx], '--', color = 'red')
# plt.plot(colonne_alpha, y, '--', color = 'red')
# plt.plot(colonne_alphatrim, [CLE[0], CLE[nb-1]], '--', color = 'red')
# plt.annotate("alphae = "+str(alpha[idx]), xy = (alpha[idx], CLE[0]), xytext = (0.175, -0.5),
#              arrowprops = {"facecolor" : "black", "width" : 0, "headwidth" : 7}, color = 'red')
# plt.annotate("cl = "+str(Cl), xy = (alpha[0], Cl), xytext = (-0.1, 1.0),
#              arrowprops = {"facecolor" : "black", "width" : 0, "headwidth" : 7}, color = 'red')
# plt.annotate("alpha = "+str(alphatrim_0), xy = (valeur_trim[0], CLE[0]), xytext = (-0.1, -0.5),
#              arrowprops = {"facecolor" : "black", "width" : 0, "headwidth" : 7}, color = 'red')
# plt.xlim(alpha[0], alpha[nb-1])
# plt.ylim(CLE[0], CLE[nb-1])
# plt.legend()
# plt.show()
#
#
# def dpha_alpha(va,q,P,alpha,nb,ms):
#     angle_alpha = np.linspace(alpha[0], alpha[1], nb)
#     P.set_mass_and_static_margin(P.m_k, ms)
#     dphr = np.array([-(P.Cm0 - ms * P.CLa * (k - P.a0) + P.Cmq * P.lt / va * q )/ P.Cmd  for k in angle_alpha])
#     return angle_alpha, dphr
#
# alpha, dphr = dpha_alpha(va_trim, 0, avion, interv_alpha, 100, ms_trim)
#
#
# dphrtrim_1 = valeur_trim[1]
# dphrtrim_1 = float(int(dphrtrim_1*1000 + 0.5))/1000
#
#
# plt.plot(alpha, dphr, label = "ms = "+str(ms))
# plt.xlabel("alpha (rad)")
# plt.ylabel("dphr (rad)")
# plt.title("dphr en fonction \n d'alpha")
# ligne_dphr = [dphr[idx]]*idx
# dphr[idx] = float(int(dphr[idx]*1000 + 0.5))/1000
# ligne_dphrtrim = [dphrtrim_1]*nb
# alpha[idx] = float(int(alpha[idx]*1000 + 0.5))/1000
# plt.plot(colonne_alphatrim, [dphr[0], valeur_trim[1]], '--', color = 'red')
# plt.plot(x, ligne_dphr, '--', color = 'red')
# plt.plot(alpha, ligne_dphrtrim, '--', color = 'red')
# plt.annotate("dphr = "+str(dphrtrim_1), xy = (alpha[0], valeur_trim[1]), xytext = (-0.1, -0.1),
#              arrowprops = {"facecolor" : "black", "width" : 0, "headwidth" : 7}, color = 'red')
# plt.annotate("dphre = "+str(dphr[idx]), xy = (alpha[0], dphr[idx]), xytext = (-0.1, -0.08),
#              arrowprops = {"facecolor" : "black", "width" : 0, "headwidth" : 7}, color = 'red')
# plt.annotate("alphae = "+str(alpha[idx]), xy = (alpha[idx], dphr[0]), xytext = (0.2, -0.12),
#              arrowprops = {"facecolor" : "black", "width" : 0, "headwidth" : 7}, color = 'red')
# plt.xlim(alpha[0], alpha[nb-1])
# plt.ylim(dphr[0], dphr[nb-1])
# plt.legend()
# plt.show()
#
# #mvgammapoint = L-mgcos(gamma) Fsin(alpha) est négligé
# #L=mgcos(gamma)=mg
# #1/2rhoSv**2CL=mg

''''''''''''''''''

#simule la trajectoire de l'avion sur 100 secondes
avion.set_mass_and_static_margin(km[0],ms[0])
va_trim = dynamic.va_of_mach(Ma_trim,h_trim)
X_trim, U_trim = dynamic.trim(avion, {'va':va_trim, 'gamm':0, 'h':h_trim})
time=np.array([i for i in range(1,100)])
traj=odeint(dynamic.dyn, X_trim,time,args=(U_trim,avion))
fig=plt.figure()
dynamic.plot(time,traj,figure=fig)
plt.show()

''''''''''''''''''''

valeurspropres=np.zeros((len(h),len(Ma), len(ms),len(km),4),dtype=complex)

for lb,l in enumerate(h):
    for kb, k in enumerate(Ma):
        for ib,i in enumerate(ms):
            for jb,j in enumerate(km):
                avion.set_mass_and_static_margin(j, i)
                va_h0= dynamic.va_of_mach(k, l, k=1.4, Rs=287.05)
                X,U = dynamic.trim(avion, {'va':va_h0,'gamm': 0.,'h': l})
                A_6, B_6 = utils.num_jacobian(X, U, avion, dynamic.dyn)
                A_4,B_4 = A_6[2:,2:],B_6[2:,2:]
                lvaleurp, lvecteurp = np.linalg.eig(A_4)
                for m, n in enumerate(lvaleurp):
                    valeurspropres[lb, kb, ib, jb, m] = n



for kb, k in enumerate(Ma):
    for ib,i in enumerate(ms):
        for jb,j in enumerate(km):
            for lb, l in enumerate(h):
                fig,ax = plt.subplots()
                ax.scatter(valeurspropres[lb,kb,ib,jb].real,valeurspropres[lb,kb,ib,jb].imag)
                ax.grid()
                plt.show()
                fig, ax = plt.subplots()
                ax.scatter(valeurspropres[lb, kb, ib, jb].real, valeurspropres[lb, kb, ib, jb].imag)
                ax.grid()
                plt.show()

