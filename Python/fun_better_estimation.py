#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 14:41:52 2022

@author: Flavien
"""
import fun_resolution
import numpy as np
#from scipy import optimize
import fonctions
import time

def better_estimation(simulation,estimation,solution):
    n = 0
    itermax = 10
    Coefs = [estimation.coefficients]
    while n<itermax:
        print(n)
        #estimation.Lambda = np.zeros((2,simulation.P,simulation.N))
        #estimation.Psi[:,0,:] = solution.Psi_0
        Coefs_ini = estimation.coefficients
        start_time =time.time()
        fun_resolution.resolution(simulation,estimation)
        Int_approx = fonctions.G(simulation,estimation)
        print(Int_approx)
        print("--- %s seconds ---" % (time.time() - start_time))
        Grad_approx = np.zeros((5,5))
        Vec_zeros = np.zeros(5)
        for j in range(0,5):
            vec_j = Vec_zeros 
            vec_j[j] = 1.0
            eps = 1.49e-8
            x_j = Coefs_ini + eps*vec_j
            start_time = time.time()
            #estimation.Lambda = np.zeros((2,simulation.P,simulation.N))
            #estimation.Psi[:,0,:] = solution.Psi_0
            F_j = fonctions.raccourci_G(simulation,estimation,x_j)
            #estimation.coefficients = Coefs_ini
            print("--- %s seconds ---" % (time.time() - start_time))
            Grad_approx[:,j] = (F_j-Int_approx)/eps
        #def G_coef(x,i): return fonctions.raccourci_G(simulation,estimation,x)[i]
        #G_coef_i = lambda x : G_coef(x,i)
        #for i in range(0,5):
        #    Grad_approx[i,:] = optimize.approx_fprime(estimation.coefficients,G_coef_i,1.49e-8)
        #Coefs_linsolve = np.linalg.solve(Grad_approx, -Int_approx+np.matmul(Grad_approx,Coefs_ini))
        print(Grad_approx)
        Coefs_linsolve = np.linalg.solve(Grad_approx, -Int_approx)+Coefs_ini
        start_time =time.time()
        Coefs_new = fonctions.proj(simulation,Coefs_linsolve)
        print("--- %s seconds ---" % (time.time() - start_time))
        estimation.coefficients = Coefs_new
        Coefs = np.append(Coefs,[Coefs_new],axis=0)
        n+=1
        if np.linalg.norm(Coefs_new-Coefs_ini)/np.linalg.norm(Coefs_ini)<1e-3:
            print(n)
            break
    return Coefs