#-*- coding: utf-8 -*-
import dynamic, utils
import numpy as np
import matplotlib.pyplot as plt
import math


nb = 100
h = [3000.0,10000.0] #m
Ma = [0.4,0.9]
ms = [-0.2,0.5]
km = [0.1,1]

avion = dynamic.Param_A321()
rho0 = utils.isa(0)[1]

alphatrim=np.zeros((len(h),len(Ma), len(ms),len(km)))
dphrtrim= np.zeros((len(h),len(Ma), len(ms),len(km)))
dthtrim = np.zeros((len(h),len(Ma), len(ms),len(km)))
prop =  np.zeros((len(h),len(Ma), len(ms),len(km)))

for ib,i in enumerate(ms):
    for jb,j in enumerate(km):
        avion.set_mass_and_static_margin(i, j)
        for kb,k in enumerate(Ma):
            for lb,l in enumerate(h):
                va_h0= dynamic.va_of_mach(k, l, k=1.4, Rs=287.05)
                X,U = dynamic.trim(avion, {'va':va_h0,'gamm': 0.,'h': l})
                alphatrim[ib,jb,kb,lb] = X[3]
                dphrtrim[ib,jb,kb,lb] = U[0]
                dthtrim[ib,jb,kb,lb] = U[1]
                prop[ib,jb,kb,lb] = dynamic.propulsion_model(X, U, avion)

def affiche(valeur,titre,xlabel,ylabel):
    for i in range(len(ms)):
        for j in range(len(km)):
            for k in range(len(Ma)):
                    plt.title(titre)
                    plt.plot(h,valeur[:,k,i,j],label = 'ms='+str(ms[i])+', km='+str((km[j])))
                    plt.ylabel(xlabel)
                    plt.xlabel(ylabel)
                    plt.legend()
    plt.show()

affiche(alphatrim,'alpha en fonction de l\'altitude','angle en radian','altitude')
affiche(dphrtrim,'dphr en fonction de l\'altitude','dphr','altitude')
affiche(dthtrim,'dthr en fonction de l\'altitude','dth en radian','altitude')

for i in range(len(ms)):
    for j in range(len(km)):
        for k in range(len(h)):
            plt.title('prop')
            plt.plot(Ma, prop[k,:, i, j], label='ms='+str(ms[i])+', km='+str((km[j])))
            plt.ylabel('prop')
            plt.xlabel('Mach')
            plt.legend()
plt.show()

valeur_trim = alphatrim[0,0,0,0], dphrtrim[0,0,0,0], dthtrim[0,0,0,0]
point_trim = [3000.0, 0.4, -0.2, 0.1]

h, Ma, ms, km = point_trim
p, rho, T = utils.isa(h)
avion.set_mass_and_static_margin(km, ms)
va = dynamic.va_of_mach(Ma, h, k = 1.4, Rs = 287.05)
Cl = (2)*avion.m*avion.g/(rho*avion.S*va**2)
print (Cl)


interv_alpha = [-10*math.pi/180,20*math.pi/180]

def CLE_alpha(va,q,P,alpha,nb,ms):
    angle_alpha = np.linspace(interv_alpha[0],interv_alpha[1],nb)
    dphr = np.array([-(P.Cm0 - ms * P.CLa * (k - P.a0) + P.Cmq * P.lt / va * q) / P.Cmd for k in angle_alpha])
    CLE = np.array([dynamic.get_aero_coefs(va, alpha, q, dphr[i], P)[0] for i,alpha in enumerate(angle_alpha)])
    return angle_alpha, CLE

alpha, CLE = CLE_alpha(va,0,avion,interv_alpha,100,ms)



plt.plot(alpha,CLE,label = "ms = "+str(ms))
plt.xlabel("Alpha (rad)")
plt.ylabel("CLE")
plt.title("CLE en fonction \n d'alpha ")
ligne_cl = [Cl]*nb
idx = np.argwhere(np.diff(np.sign(CLE - ligne_cl)) != 0).reshape(-1)
print(idx)
x = alpha[:idx]
colonne_alpha = [alpha[idx]]*2
y = [CLE[0], Cl]
plt.plot(x, ligne_cl[:idx], '--', color = 'red')
plt.plot(colonne_alpha, y, '--', color = 'red')
idx = idx[0]
plt.annotate("alpha = "+str(alpha[idx]), xy = (alpha[idx], -1.0), xytext = (0.175, -0.5),
             arrowprops = {"facecolor" : "black", "width" : 0, "headwidth" : 7}, color = 'red')
plt.annotate("cl = "+str(Cl), xy = (alpha[0], Cl), xytext = (-0.1, 1.0),
             arrowprops = {"facecolor" : "black", "width" : 0, "headwidth" : 7}, color = 'red')
plt.legend()
plt.show()
print(alphatrim[0,0,0,0])

