# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 09:59:56 2024

@author: Flavien

L'objectif est de résoudre numériquement sur un maillage en 2D le problème inverse associé à l'EDP :
    d/dt A- d^2/dx^2 A - d^2/dy^2 A = mu * f1(A,B)
    d/dt B- d * (d^2/dx^2 B - d^2/dy^2 B) = mu * f2(A,B)
    avec des conditions initiales en temps et des conditions de Neumann sur les bords
    
    
La résolution du problème inverse revient à résoudre :
    (PF1)
        d/dt A- d^2/dx^2 A - d^2/dy^2 A = mu * f1(A,B) + (KqqL)1
        d/dt B- d * (d^2/dx^2 B - d^2/dy^2 B) = mu * f2(A,B) + (KqqL)2
        A(t=0,x,y) = A0(x,y) + (KaaL^0)1
        B(t=0,x,y) = B0(x,y) + (KaaL^0)2
        Neumann homogène aux bords
    (PF2)
        d/dt L1- d^2/dx^2 L1 - d^2/dy^2 L1 = -mu * d/dA f1(A,B) * L1 -mu * d/dB f1(A,B) * L2  + C(A,B)1
        d/dt L2- d * (d^2/dx^2 L2 - d^2/dy^2 L2) = -mu * d/dA f2(A,B) * L1 -mu * d/dB f2(A,B) * L2  + C(A,B)2
        L1(t=tf,x,y) = L2(t=tf,x,y) = 0
        Neumann homogène aux bords
    (PF3)
        integrale_I integrale_Omega [0, f1(A,B);(d^2/dx^2 B - d^2/dy^2 B), f2(A,B)]^T . [L1 ; L2] = [0 ; 0]
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.sparse import diags
from small_functions import *
from big_functions import *
from calc_PF1 import *
from calc_PF2 import *
from calc_PF3 import *
from calc_PF1_and_PF2_for_eta_known import *
from Newton_Raphson_proj import *
from save_animation import *
from math import *
from graphique import *

# On travaille sur le carré [0,1]*[0,1] entre t0 et tf:
##Espace : 
N = 21;
x = np.linspace(0,1,N); deltax = 1/(N-1);
y = np.linspace(0,1,N);
xx, yy = np.meshgrid(x,y)
Vol_Omega = 1
##Temps :
deltat = 0.001
tf = 1.0
Ntime = int(tf/deltat)+1

#Coefficients d et mu :

Coef_d = 23.23;    Coef_d_real = 20;
Coef_mu = 26.75;   Coef_mu_real = 30;

#Hyperparamètres : 
sigma_eps = 10.0
sigma_a = 0.0001
sigma_q = 0.0001
lb = np.array([15,25])
ub = np.array([25,35])

#Conditions initiales : 

A0 = 9.9338+0.01*(np.cos(np.pi*xx)*np.cos(np.pi*yy))
B0 = 9.2892+0.01*(np.cos(np.pi*xx)*np.cos(np.pi*yy))

A0_vec = A0.reshape(-1)
B0_vec = B0.reshape(-1)

fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 2, 1, projection='3d')
surf = ax.plot_surface(xx, yy, A0)
ax = fig.add_subplot(1, 2, 2, projection='3d')
surf = ax.plot_surface(xx, yy, B0)
plt.show()

#Création de la matrice discrète du laplacien en 2D
##On met la structure avec des conditions aux bords cycliques ce qui marche avec du Neumann homogène
Mat_Lap = diags([-1,-1,4,-1,-1],offsets=[-N,-1,0,1,N],shape=(N*N,N*N))
Mat_Lap = Mat_Lap + diags([-1,-1],offsets=[-N*(N-1),N*(N-1)],shape=(N*N,N*N))
Mat_Lap = Mat_Lap + diags([-1,-1],offsets=[-N*N+1,N*N-1],shape=(N*N,N*N))
Mat_Lap_for_A = sp.sparse.eye(N*N) + deltat/(deltax**2)*Mat_Lap
Mat_Lap_for_B = sp.sparse.eye(N*N) + Coef_d*deltat/(deltax**2)*Mat_Lap
Mat_Lap_for_A_csc = Mat_Lap_for_A.tocsc()
Mat_Lap_for_B_csc = Mat_Lap_for_B.tocsc()
Mat_Lap_for_A_array = Mat_Lap_for_A.toarray()
Mat_Lap_for_B_array = Mat_Lap_for_B.toarray()

##Pour aider le temps de calcul de KqqL
vec_time = np.array([np.linspace(0, tf,num=Ntime)])
Mat_smart_for_kqql = np.minimum(vec_time,np.transpose(vec_time))

# Fabrication des données :
KqqL1 = np.zeros((N*N,Ntime))#.reshape(-1)
KqqL2 = np.zeros((N*N,Ntime))#.reshape(-1)
N_mea = 20
Mat_Lap_for_B_corrected = sp.sparse.eye(N*N) + Coef_d_real*deltat/(deltax**2)*Mat_Lap
(A_all_time,B_all_time) = calc_PF1(A0_vec,B0_vec,deltax,N,deltat,tf,Coef_mu_real,Coef_d_real,Mat_Lap,KqqL1,KqqL2)
print('Données done')

A_exact = A_all_time
B_exact = B_all_time

t_mea = np.linspace(0,tf,N_mea+1,endpoint=False)[1:]
index_mea = np.round(t_mea/deltat).astype(int)
D_mea_A = A_all_time[:,index_mea]
D_mea_B = B_all_time[:,index_mea]


# #Première résolution de (PF1) : 
# KqqL1 = np.zeros((N*N,Ntime))#.reshape(-1)
# KqqL2 = np.zeros((N*N,Ntime))#.reshape(-1)

# (A_all_time,B_all_time) = calc_PF1(A0_vec,B0_vec,N,deltat,tf,Coef_mu,Mat_Lap_for_A_csc,KqqL1,Mat_Lap_for_B_csc,KqqL2)
# print('PF1 done')
# A_tf = A_next_vec.reshape((N,N))
# B_tf = B_next_vec.reshape((N,N))
A_tf = A_all_time[:,-1].reshape((N,N))
B_tf = B_all_time[:,-1].reshape((N,N))

fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 2, 1, projection='3d')
surf = ax.plot_surface(xx, yy, A_tf)
ax = fig.add_subplot(1, 2, 2, projection='3d')
surf = ax.plot_surface(xx, yy, B_tf)
plt.show()

# # Première résolution de (PF2) :
    
# (L1_all_time,L2_all_time) = calc_PF2(A_all_time,B_all_time,index_mea,D_mea_A,D_mea_B,sigma_eps,N,deltax,deltat,tf,Coef_mu,Coef_d,Mat_Lap)
# print('PF2 done')

eta = (Coef_d,Coef_mu)
(A_eta,B_eta,L1_eta,L2_eta) = calc_PF1_PF2(eta,N,deltax,deltat,tf,Vol_Omega,Mat_smart_for_kqql,index_mea,D_mea_A,D_mea_B,sigma_eps,sigma_q,sigma_a,Mat_Lap,A0_vec,B0_vec)

# print('Sauvegarde de A_eta')
# save_animation(A_eta,'A_eta')

# print('Sauvegarde de B_eta')
# save_animation(B_eta,'B_eta')

# print('Sauvegarde de L1_eta')
# save_animation(L1_eta,'L1_eta')

# print('Sauvegarde de L2_eta')
# save_animation(L2_eta,'L2_eta')

#%%
eta_0 = np.array([Coef_d,Coef_mu])
all_eta,A_final,B_final,L1_final,L2_final = Methode_NP_proj(lb,ub,eta_0,N,deltax,deltat,tf,Vol_Omega,Mat_Lap_for_A,Mat_Lap_for_B,Mat_smart_for_kqql,index_mea,D_mea_A,D_mea_B,sigma_eps,sigma_q,sigma_a,Mat_Lap,A0_vec,B0_vec)
# print('Sauvegarde de A_eta')
# save_animation(A_final,'A_final')

# print('Sauvegarde de B_eta')
# save_animation(B_final,'B_final')

# print('Sauvegarde de L1_eta')
# save_animation(L1_final,'L1_final')

# print('Sauvegarde de L2_eta')
# save_animation(L2_final,'L2_final')
#%% 
graphique(A_exact,B_exact,A_eta,B_eta,L1_eta,L2_eta,all_eta,A_final,B_final,L1_final,L2_final,tf,N)