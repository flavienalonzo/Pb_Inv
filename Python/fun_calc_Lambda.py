#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 16:51:15 2022

@author: Flavien
"""
import fonctions
import numpy as np
import scipy.sparse
import scipy.sparse.linalg

def calc_lambda(simulation,estimation):
    L1_next = np.zeros(simulation.N)
    L2_next = np.zeros(simulation.N)
    L1 = [L1_next]
    L2 = [L2_next]
    Mea = fonctions.measurement(simulation,estimation)
    for k in range(simulation.P-1,0,-1):
        (L1_inter,L2_inter) = fonctions.calc_Runge_Kutta_4_inv(simulation,estimation,k,L1_next,L2_next)
        #calc_L1
        Mat1 = simulation.diff_u*simulation.dt/pow(simulation.h,2)*scipy.sparse.csc_matrix(2*np.eye(simulation.N)-np.diag(np.ones(simulation.N-1),-1)-np.diag(np.ones(simulation.N-1),1))
        M_vit = np.diag(np.ones(simulation.N-1),1) - np.diag(np.ones(simulation.N-1),-1)
        M_vit[0,0]=-1
        M_vit[simulation.N-1,simulation.N-1] = 1
        M_vit = scipy.sparse.csc_matrix(M_vit)
        Vit = 0.5/simulation.h*M_vit.dot(estimation.Psi[1,k,:])
        Mat2 = simulation.diff_u*simulation.chi*simulation.dt*0.5/simulation.h*scipy.sparse.diags(Vit).multiply(M_vit)
        Mat3 = scipy.sparse.csc_matrix(-Mat1+Mat2-scipy.sparse.eye(simulation.N))
        if k in simulation.Measures[1]:
            Corr = simulation.W_eps.dot(np.append(simulation.Measures[2][int(np.where(simulation.Measures[1]==k)[0]),:],simulation.Measures[3][int(np.where(simulation.Measures[1]==k)[0]),:])-np.append(Mea[0,int(np.where(simulation.Measures[1]==k)[0]),:],Mea[1,int(np.where(simulation.Measures[1]==k)[0]),:]))
            L1_prev = scipy.sparse.linalg.spsolve(Mat3,-L1_inter-simulation.dt*Corr[0:simulation.N])
        else:
            L1_prev = scipy.sparse.linalg.spsolve(Mat3,-L1_inter)                
        #calc_L2
        Mata = simulation.diff_c*simulation.dt/pow(simulation.h,2)*scipy.sparse.csc_matrix(2*np.eye(simulation.N)-np.diag(np.ones(simulation.N-1),-1)-np.diag(np.ones(simulation.N-1),1))
        Vec_plus = scipy.sparse.csc_matrix(np.eye(simulation.N)-np.diag(np.ones(simulation.N-1),1)).dot(L1_prev)
        Vec_plus[simulation.N-1]=0.0
        Vec_moins = scipy.sparse.csc_matrix(np.eye(simulation.N)-np.diag(np.ones(simulation.N-1),-1)).dot(L1_prev)
        Vec_moins[0] = 0.0
        Vec_zeros = np.zeros(simulation.N)
        L1_plus_pos = np.maximum(Vec_plus,Vec_zeros)
        L1_plus_neg = np.maximum(-Vec_plus,Vec_zeros)
        L1_moins_pos = np.maximum(Vec_moins,Vec_zeros)
        L1_moins_neg = np.maximum(-Vec_moins,Vec_zeros)
        Matb = simulation.diff_u*simulation.chi*simulation.dt*scipy.sparse.csc_matrix(-np.diag(L1_plus_neg)-np.diag(L1_moins_neg)+np.diag(L1_plus_pos[0:simulation.N-1],1)+np.diag(L1_moins_pos[1:simulation.N],-1))/pow(simulation.h,2)
        Conv = Matb.dot(estimation.Psi[0,k,:])
        Matc = scipy.sparse.csc_matrix(-Mata-scipy.sparse.eye(simulation.N))
        if k in simulation.Measures[1]:
            L2_prev = scipy.sparse.linalg.spsolve(Matc,-L2_inter-Conv-simulation.dt*Corr[simulation.N:2*simulation.N])
        else:
            L2_prev = scipy.sparse.linalg.spsolve(Matc,-L2_inter-Conv) 
        L1 = np.append(L1,[L1_prev],axis=0)
        L2 = np.append(L2,[L2_prev],axis=0)
        L1_next = L1_prev
        L2_next = L2_prev
    estimation.Lambda = np.flip(np.array([L1,L2]),axis=1)