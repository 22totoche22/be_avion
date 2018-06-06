#-*- coding: utf-8 -*-
import dynamic, utils
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint

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

point_trim = [3000.0, 0.4, -0.2, 0.1]

h, Ma, ms, km = point_trim
p, rho, T = utils.isa(h)
avion.set_mass_and_static_margin(km, ms)
va = dynamic.va_of_mach(Ma, h, k = 1.4, Rs = 287.05)
Cl = (2)*avion.m*avion.g/(rho*avion.S*va**2)
print (Cl)





interv_alpha = [-10*math.pi/180,20*math.pi/180]
alpha, cl = séance1.get_CL_alpha(va, 0, dphrtrim[0,0,0,0], avion, interv_alpha, 100)
ligne = [Cl]*100
colonne = [-1.0, 2.5]

alphatrim = [alphatrim[0,0,0,0]]*2



plt.plot(alpha, cl, label = "dphr="+str(dphrtrim[0,0,0,0])+" rad")
plt.title("Cl en fonction de alpha")
plt.xlabel("alpha (rad)")
plt.ylabel("Cl")
plt.plot(alpha, ligne, '--')
plt.plot(alphatrim, colonne, '--')

plt.legend()
plt.show()

print(alphatrim[0,0,0,0])


point_trim = [alphatrim[0,0,0,0],dphrtrim[0,0,0,0],dthtrim[0,0,0,0]]

#mvgammapoint = L-mgcos(gamma) Fsin(alpha) est négligé
#L=mgcos(gamma)=mg
#1/2rhoSv**2CL=mg


time=np.array([i for i in range(1,100)])
traj=odeint(dynamic.dyn, X,time,args=(U,avion))
fig=plt.figure()
dynamic.plot(time,traj,figure=fig)
plt.show()
#valeurspropres,vecteurpropres=numpy.linalg.eig(A)

vector_select = np.array([0,0,1,1,1,1])
A_6, B_6 = utils.num_jacobian(X, U, avion, dynamic.dyn)

A_4, B_4 = np.dot(vector_select,A_6),np.dot(vector_select,B_6)
print(A_6)
print(A_4)

