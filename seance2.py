#-*- coding: utf-8 -*-
import dynamic, utils
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint

h = [3000,10000] #altitude en m
Ma = [0.4,0.9] #intervalel nombre de mach
ms = [-0.2,0.5] #marge statique
km = [0.1,1] #coefficient de réglage de la masse
nb = 100 #nb de points pour les courbes

avion = dynamic.Param_A321()
rho0 = utils.isa(0)[1]

h1 = np.linspace(h[0],h[1],nb)
Ma1 = np.linspace(Ma[0], Ma[1], nb)

''''''''''''

#calcul des points de Trim

alphatrim=np.zeros((len(ms),len(km),len(Ma),len(h1)))
dphrtrim= np.zeros((len(ms),len(km),len(Ma),len(h1)))
dthtrim = np.zeros((len(ms),len(km),len(Ma),len(h1)))
prop = np.zeros((len(ms),len(km),len(Ma1),len(h)))

for ib,i in enumerate(ms):
    for jb,j in enumerate(km):
        avion.set_mass_and_static_margin(j, i)
        for kb,k in enumerate(Ma):
            for lb,l in enumerate(h1):
                va_h0= dynamic.va_of_mach(k, l, k=1.4, Rs=287.05)
                X,U = dynamic.trim(avion, {'va':va_h0,'gamm': 0.,'h': l})
                alphatrim[ib,jb,kb,lb] = X[3]
                dphrtrim[ib,jb,kb,lb] = U[0]
                dthtrim[ib,jb,kb,lb] = U[1]
        for kbb, k in enumerate(Ma1):
            for lbb, l in enumerate(h):
                va_h0b = dynamic.va_of_mach(k, l, k=1.4, Rs=287.05)
                Xb, Ub = dynamic.trim(avion, {'va': va_h0b, 'gamm': 0., 'h': l})
                prop[ib,jb,kbb,lbb] = dynamic.propulsion_model(Xb, Ub, avion)


def affiche():
    for i in range(len(ms)):
        for j in range(len(km)):
            for k in range(len(Ma)):
                        plt.suptitle('ms=' + str(ms[i]) + ', km=' + str((km[j])))
                        plt.subplot(1,4,1)
                        plt.title("alpha")
                        plt.plot(h1,utils.deg_of_rad(alphatrim[i,j,k,:]),label = "Mach = "+str(Ma[k])+"")
                        plt.ylabel('alpha en degre')
                        plt.xlabel('altitude')
                        plt.legend()
                        plt.subplot(1, 4, 2)
                        plt.title("dphr")
                        plt.plot(h1, utils.deg_of_rad(dphrtrim[i, j, k, :]), label = "Mach = "+str(Ma[k])+"")
                        plt.ylabel('dphr en degre')
                        plt.xlabel('altitude')
                        plt.legend()
                        plt.subplot(1, 4, 3)
                        plt.title("dth")
                        plt.plot(h1, utils.deg_of_rad(dthtrim[i, j, k, :]), label = "Mach = "+str(Ma[k])+"")
                        plt.ylabel('dth en degre')
                        plt.xlabel('altitude')
                        plt.legend()
            for l in range(len(h)):
                        plt.subplot(1, 4, 4)
                        plt.title("poussee")
                        plt.plot(Ma1, prop[i, j,: ,l ], label="altitude = " + str(h[l]) + "")
                        plt.ylabel('poussee')
                        plt.xlabel('Mach')
                        plt.legend()
            plt.show()

affiche()

''''''''''''''''''
#
valeur_trim = alphatrim[0,0,0,0], dphrtrim[0,0,0,0], dthtrim[0,0,0,0]
point_trim = [3000.0, 0.4, -0.2, 0.1]

h_trim, Ma_trim, ms_trim, km_trim = point_trim
p_trim, rho_trim, T_trim = utils.isa(h_trim)


avion.set_mass_and_static_margin(km_trim, ms_trim)
va_trim = dynamic.va_of_mach(Ma_trim, h_trim, k = 1.4, Rs = 287.05)
Cl = 2*avion.m*avion.g/(rho_trim*avion.S*va_trim**2)

interv_alpha = [-10*math.pi/180,20*math.pi/180]

def CLE_alpha(va,q,P,alpha,nb,ms):
    angle_alpha = np.linspace(interv_alpha[0],interv_alpha[1],nb)
    dphr = np.array([-(P.Cm0 - ms * P.CLa * (k - P.a0) + P.Cmq * P.lt / va * q) / P.Cmd for k in angle_alpha])
    CLE = np.array([dynamic.get_aero_coefs(va, alpha, q, dphr[i], P)[0] for i,alpha in enumerate(angle_alpha)])
    return angle_alpha, CLE

alpha, CLE = CLE_alpha(va_trim,0,avion,interv_alpha,100,ms_trim)


alphatrim_0 = valeur_trim[0]
alphatrim_0 = float(int(alphatrim_0*1000 + 0.5))/1000
Cl = float(int(Cl*1000 +0.5))/1000

plt.plot(alpha,CLE,label = "ms = "+str(ms_trim))
plt.xlabel("Alpha (rad)")
plt.ylabel("CLE")
plt.title("CLE en fonction \n d'alpha ")
ligne_cl = [Cl]*nb
idx = np.argwhere(np.diff(np.sign(CLE - ligne_cl)) != 0).reshape(-1)
idx = idx[0]
x = alpha[:idx]
alpha[idx] = float(int(alpha[idx]*1000 + 0.5))/1000
colonne_alpha = [alpha[idx]]*2
y = [CLE[0], Cl]
colonne_alphatrim = [alphatrim_0]*2
plt.plot(x, ligne_cl[:idx], '--', color = 'red')
plt.plot(colonne_alpha, y, '--', color = 'red')
plt.plot(colonne_alphatrim, [CLE[0], CLE[nb-1]], '--', color = 'red')
plt.annotate("alphae = "+str(alpha[idx]), xy = (alpha[idx], CLE[0]), xytext = (0.175, -0.5),
             arrowprops = {"facecolor" : "black", "width" : 0, "headwidth" : 7}, color = 'red')
plt.annotate("cl = "+str(Cl), xy = (alpha[0], Cl), xytext = (-0.1, 1.0),
             arrowprops = {"facecolor" : "black", "width" : 0, "headwidth" : 7}, color = 'red')
plt.annotate("alpha = "+str(alphatrim_0), xy = (valeur_trim[0], CLE[0]), xytext = (-0.1, -0.5),
             arrowprops = {"facecolor" : "black", "width" : 0, "headwidth" : 7}, color = 'red')
plt.xlim(alpha[0], alpha[nb-1])
plt.ylim(CLE[0], CLE[nb-1])
plt.legend()
plt.show()


def dpha_alpha(va,q,P,alpha,nb,ms):
    angle_alpha = np.linspace(alpha[0], alpha[1], nb)
    P.set_mass_and_static_margin(P.m_k, ms)
    dphr = np.array([-(P.Cm0 - ms * P.CLa * (k - P.a0) + P.Cmq * P.lt / va * q )/ P.Cmd  for k in angle_alpha])
    return angle_alpha, dphr

alpha, dphr = dpha_alpha(va_trim, 0, avion, interv_alpha, 100, ms_trim)


dphrtrim_1 = valeur_trim[1]
dphrtrim_1 = float(int(dphrtrim_1*1000 + 0.5))/1000


plt.plot(alpha, dphr, label = "ms = "+str(ms))
plt.xlabel("alpha (rad)")
plt.ylabel("dphr (rad)")
plt.title("dphr en fonction \n d'alpha")
ligne_dphr = [dphr[idx]]*idx
dphr[idx] = float(int(dphr[idx]*1000 + 0.5))/1000
ligne_dphrtrim = [dphrtrim_1]*nb
alpha[idx] = float(int(alpha[idx]*1000 + 0.5))/1000
plt.plot(colonne_alphatrim, [dphr[0], valeur_trim[1]], '--', color = 'red')
plt.plot(x, ligne_dphr, '--', color = 'red')
plt.plot(alpha, ligne_dphrtrim, '--', color = 'red')
plt.annotate("dphr = "+str(dphrtrim_1), xy = (alpha[0], valeur_trim[1]), xytext = (-0.1, -0.1),
             arrowprops = {"facecolor" : "black", "width" : 0, "headwidth" : 7}, color = 'red')
plt.annotate("dphre = "+str(dphr[idx]), xy = (alpha[0], dphr[idx]), xytext = (-0.1, -0.08),
             arrowprops = {"facecolor" : "black", "width" : 0, "headwidth" : 7}, color = 'red')
plt.annotate("alphae = "+str(alpha[idx]), xy = (alpha[idx], dphr[0]), xytext = (0.2, -0.12),
             arrowprops = {"facecolor" : "black", "width" : 0, "headwidth" : 7}, color = 'red')
plt.xlim(alpha[0], alpha[nb-1])
plt.ylim(dphr[0], dphr[nb-1])
plt.legend()
plt.show()
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
#calcul des valeurs propres


nb2 = 10
h2 = np.linspace(h[0],h[1],nb2)
Ma2 = np.linspace(Ma[0], Ma[1], nb2)
ms2 = np.linspace(ms[0],ms[1],nb2)
km2 = np.linspace(km[0],km[1],nb2)

valeurspropres=np.zeros((len(h2),len(Ma2), len(ms2),len(km2),4),dtype=complex)


for lb,l in enumerate(h2):
    for kb, k in enumerate(Ma2):
        for ib,i in enumerate(ms2):
            for jb,j in enumerate(km2):
                avion.set_mass_and_static_margin(j, i)
                va_h0= dynamic.va_of_mach(k, l, k=1.4, Rs=287.05)
                X,U = dynamic.trim(avion, {'va':va_h0,'gamm': 0.,'h': l})
                A_6, B_6 = utils.num_jacobian(X, U, avion, dynamic.dyn)
                A_4 = A_6[2:,2:]
                lvaleurp, lvecteurp = np.linalg.eig(A_4)
                for m, n in enumerate(lvaleurp):
                    valeurspropres[lb, kb, ib, jb, m] = n

def decorate(vp_Va, vp_a, vp_theta, vp_q, ss,nss,title,inter, val):
    plt.suptitle('${} \in [{},{}]$'.format(val, inter[0],inter[1]))
    plt.subplot(1, nss, ss)
    plt.title(title)
    plt.plot(vp_Va.real, vp_Va.imag, label="V_a")
    plt.plot(vp_a.real, vp_a.imag, label="alpha")
    plt.plot(vp_theta.real, vp_theta.imag, label="theta")
    plt.plot(vp_q.real, vp_q.imag, label="q")
    plt.plot([vp_Va[0].real], [vp_Va[0].imag], marker='x')
    plt.plot([vp_a[0].real], [vp_a[0].imag], marker='x')
    plt.plot([vp_theta[0].real], [vp_theta[0].imag], marker='x')
    plt.plot([vp_q[0].real], [vp_q[0].imag], marker='x')
    plt.xlabel("Re")
    plt.ylabel("Im")
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.grid(linestyle='--')
    plt.legend()

def affiche_vp():
    ss = 1
    nomb = len(Ma)*len(ms)*len(km)
    for kb, k in enumerate(Ma):
        for ib,i in enumerate(ms):
            for jb,j in enumerate(km):
                    vp_Va = valeurspropres[:,nb2 * kb - 1, nb2 * ib - 1, nb2 * jb - 1, 0]
                    vp_a = valeurspropres[:, nb2 * kb - 1, nb2 * ib - 1, nb2 * jb - 1, 1]
                    vp_theta = valeurspropres[:, nb2 * kb - 1, nb2 * ib - 1, nb2 * jb - 1, 2]
                    vp_q = valeurspropres[:, nb2 * kb - 1, nb2 * ib - 1, nb2 * jb - 1, 3]
                    decorate(vp_Va, vp_a, vp_theta, vp_q,ss, nomb,'Ma={}, ms={}, km={}'.format(kb,ib,jb),[h[0], h[1]],'h')
                    ss+=1
    plt.show()
    ss=1
    for kb, k in enumerate(Ma):
        for ib,i in enumerate(ms):
            for lb, l in enumerate(h):
                vp_Va = valeurspropres[nb2 * lb -1, nb2 * kb - 1, nb2 * ib - 1, :, 0]
                vp_a = valeurspropres[nb2 * lb -1, nb2 * kb - 1, nb2 * ib - 1, :, 1]
                vp_theta = valeurspropres[nb2 * lb -1, nb2 * kb - 1, nb2 * ib - 1, :, 2]
                vp_q = valeurspropres[nb2 * lb -1, nb2 * kb - 1, nb2 * ib - 1, :, 3]
                decorate(vp_Va, vp_a, vp_theta, vp_q, ss, nomb, 'Ma={}, ms={}, h={}'.format(kb, ib, lb), [km[0], km[1]],'km')
                ss += 1
    plt.show()
    ss=1
    for kb, k in enumerate(Ma):
        for lb, l in enumerate(h):
            for jb,j in enumerate(km):
                vp_Va = valeurspropres[nb2 * lb -1, nb2 * kb - 1, :, nb2 * jb - 1, 0]
                vp_a = valeurspropres[nb2 * lb -1, nb2 * kb - 1, :, nb2 * jb - 1, 1]
                vp_theta = valeurspropres[nb2 * lb -1, nb2 * kb - 1, :, nb2 * jb - 1, 2]
                vp_q = valeurspropres[nb2 * lb -1, nb2 * kb - 1, :, nb2 * jb - 1, 3]
                decorate(vp_Va, vp_a, vp_theta, vp_q, ss, nomb, 'Ma={}, h={}, km={}'.format(kb, lb, jb), [ms[0], ms[1]],'ms')
                ss += 1
    plt.show()
    ss=1
    for lb, l in enumerate(h):
        for ib,i in enumerate(ms):
            for jb,j in enumerate(km):
                vp_Va = valeurspropres[nb2 * lb -1, :, nb2 * ib - 1, nb2 * jb - 1, 0]
                vp_a = valeurspropres[nb2 * lb -1, :, nb2 * ib - 1, nb2 * jb - 1, 1]
                vp_theta = valeurspropres[nb2 * lb -1, :, nb2 * ib - 1, nb2 * jb - 1, 2]
                vp_q = valeurspropres[nb2 * lb -1, :, nb2 * ib - 1, nb2 * jb - 1, 3]
                decorate(vp_Va, vp_a, vp_theta, vp_q, ss, nomb, 'h={}, ms={}, km={}'.format(lb, ib, jb), [Ma[0], Ma[1]],'Ma')
                ss += 1
    plt.show()


affiche_vp()


