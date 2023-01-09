#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 15:08:29 2022

@author: Flavien
"""
import fonctions
import numpy as np
import scipy.sparse
import scipy.sparse.linalg

def calc_Psi(simulation,estimation):
    Err_cond_ini = simulation.C_a.dot(np.concatenate((estimation.Lambda[0,0,:],estimation.Lambda[1,0,:]),axis=0))
    u_prev = estimation.Psi[0,0,:]+Err_cond_ini[0:simulation.N]
    c_prev = estimation.Psi[1,0,:]+Err_cond_ini[simulation.N:2*simulation.N]
    U = [u_prev] 
    C = [c_prev]
    Cross_product = fonctions.e1(simulation,estimation)
    for k in range(1,simulation.P):
        (u_inter,c_inter) = fonctions.calc_Runge_Kutta_4(simulation,estimation,u_prev,c_prev)
        #calc_c_next
        Mat = simulation.diff_c*simulation.dt*(2*np.eye(simulation.N)-np.diag(np.ones(simulation.N-1),-1)-np.diag(np.ones(simulation.N-1),1))/pow(simulation.h,2)
        Mat2 = scipy.sparse.csc_matrix(Mat+np.eye(simulation.N))
        c_next = scipy.sparse.linalg.spsolve(Mat2,c_prev+simulation.dt*Cross_product[1,k,:])
        #calc_u_next
        Mata = simulation.diff_u*simulation.dt*(2*np.eye(simulation.N)-np.diag(np.ones(simulation.N-1),-1)-np.diag(np.ones(simulation.N-1),1))/pow(simulation.h,2)
        Vec_plus = scipy.sparse.csc_matrix(np.eye(simulation.N)-np.diag(np.ones(simulation.N-1),1)).dot(c_next)
        Vec_plus[simulation.N-1]=0.0
        Vec_moins = scipy.sparse.csc_matrix(np.eye(simulation.N)-np.diag(np.ones(simulation.N-1),-1)).dot(c_next)
        Vec_moins[0] = 0.0
        Vec_zeros = np.zeros(simulation.N)
        C_plus_pos = np.maximum(Vec_plus,Vec_zeros)
        C_plus_neg = np.maximum(-Vec_plus,Vec_zeros)
        C_moins_pos = np.maximum(Vec_moins,Vec_zeros)
        C_moins_neg = np.maximum(-Vec_moins,Vec_zeros)
        Matb = simulation.diff_u*simulation.chi*simulation.dt*(-np.diag(C_plus_neg)-np.diag(C_moins_neg)+np.diag(C_plus_pos[0:simulation.N-1],1)+np.diag(C_moins_pos[1:simulation.N],-1))/pow(simulation.h,2)
        Matc = scipy.sparse.csc_matrix(Mata-Matb+np.eye(simulation.N))
        u_next =  scipy.sparse.linalg.spsolve(Matc,u_prev+simulation.dt*Cross_product[0,k,:])
        U = np.append(U,[u_next],axis=0)
        C = np.append(C,[c_next],axis=0)
        u_prev = u_next
        c_prev = c_next
    estimation.Psi = np.array([U,C])
