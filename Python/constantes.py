#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 13:51:40 2022

@author: Flavien
"""
import numpy as np
import module_libraries
import param_simu
import scipy.sparse 
import scipy.sparse.linalg
import fonctions

global simu 

simu = module_libraries.Simulation(**param_simu.param)

global h, P, x, Cq, C_eps, W_eps, C_a, W_a, solution

h = 1/(simu.N-1)
P = int((simu.tf-simu.t0)/simu.dt+1)
x = np.linspace(0,1,simu.N)
u = np.array(range(1,P+1))
v=u.reshape(-1,1)
Cq = pow(simu.sigma_q,2)*simu.dt*np.minimum(u,v)
#C_eps = scipy.sparse.eye(2*simu.N)*(pow(simu.sigma_eps,2))
C_eps = scipy.sparse.csc_matrix(pow(simu.sigma_eps,2)*scipy.sparse.identity(2*simu.N))
W_eps = scipy.sparse.linalg.inv(C_eps)
C_a = scipy.sparse.csc_matrix(pow(simu.sigma_a,2)*scipy.sparse.identity(2*simu.N))
W_a = scipy.sparse.linalg.inv(C_a)

solution = module_libraries.Solution_exacte(**param_simu.param,**param_simu.sol)

estimation = module_libraries.Estimation(**param_simu.estimation_initiale)

def calc_Psi_0():
    a = 1/2
    b = 7/10
    c = (a+b)/2
    u_0 = np.zeros(simu.N)
    for i in range(0,simu.N):
        if np.abs(x[i]-c)<(b-a)/2:
            y = (x[i]-c)/((b-a)/2)
            u_0[i] = np.exp(1-1/(1-pow(y,2)))
    c_0 = 0.5*u_0
    solution.Psi_0 = np.stack([u_0,c_0])
        

def calc_Psi_exact():
    u_prev = solution.Psi_0[0,:]
    c_prev = solution.Psi_0[1,:]
    U = [u_prev] 
    C = [c_prev]
    for k in range(1,P):
        (u_inter,c_inter) = fonctions.calc_Runge_Kutta_4(simu,estimation,u_prev,c_prev)
        c_next = fonctions.calc_c_next(simu,h,c_inter)
        u_next = fonctions.calc_u_next(simu,h,u_inter,c_next)
        U = np.append(U,[u_next],axis=0)
        C = np.append(C,[c_next],axis=0)
        u_prev = u_next
        c_prev = c_next
    solution.Psi_exact = np.array([U,C])

calc_Psi_0()
calc_Psi_exact()

def give_measurement():
    times = simu.dt*np.array(range(simu.dmea,P,simu.dmea))
    U = solution.Psi_exact[0]
    C = solution.Psi_exact[1]
    simu.Measures = (times,U[range(simu.dmea,P,simu.dmea),:],C[range(simu.dmea,P,simu.dmea),:])
    
give_measurement()