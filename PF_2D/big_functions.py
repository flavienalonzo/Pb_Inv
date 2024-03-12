# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:01:18 2024

@author: Maxence
"""
from small_functions import * 
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve,eigs
import scipy as sp
from scipy.sparse import diags,hstack,vstack
from scipy.optimize import minimize,LinearConstraint
import numpy as np

def calc_A_next(A_prev,B_prev,dt,mu,Lap,KqqL):
    second_membre = A_prev+dt*mu*f1(A_prev,B_prev)+dt*KqqL
    A_next = spsolve(Lap,second_membre)
    return A_next

def calc_B_next(A_prev,B_prev,dt,mu,Lap,KqqL):
    second_membre = B_prev+dt*mu*f2(A_prev,B_prev)+dt*KqqL
    A_next = spsolve(Lap,second_membre)
    return A_next

def calc_L1_prev(L1_next,L2_next,A_next,B_next,dt,mu,Lapp,CAB1):
    #second_membre = L1_next + dt*mu*diff_A_f1(A_next,B_next)*L1_next + dt*mu*diff_B_f1(A_next,B_next)*L2_next - dt*CAB1
    second_membre = -L1_next - dt*mu*diff_A_f1(A_next,B_next)*L1_next - dt*mu*diff_B_f1(A_next,B_next)*L2_next - dt*CAB1
    #L1_prev = spsolve(Lapp,second_membre)
    L1_prev = spsolve(Lapp,second_membre)
    return L1_prev

def calc_L2_prev(L1_next,L2_next,A_next,B_next,dt,mu,Lapp,CAB2):
    # second_membre = L2_next + dt*mu*diff_A_f2(A_next,B_next)*L1_next + dt*mu*diff_B_f2(A_next,B_next)*L2_next - dt*CAB2
    # L1_prev = spsolve(Lapp,second_membre)
    second_membre = -L2_next  - dt*mu*diff_A_f2(A_next,B_next)*L1_next - dt*mu*diff_B_f2(A_next,B_next)*L2_next - dt*CAB2
    L1_prev = spsolve(Lapp,second_membre)
    return L1_prev

def calc_L1_L2_prev(L1_next,L2_next,A_next,B_next,dt,mu,Lapp_L1,Lapp_L2,CAB1,CAB2):
    second_membre_L1 = L1_next - dt*CAB1
    second_membre_L2 = L2_next - dt*CAB2
    M11 = -Lapp_L1+dt*mu*diags(diff_A_f1(A_next,B_next))
    M12 = dt*mu*diags(diff_B_f1(A_next,B_next))
    M21 = dt*mu*diags(diff_A_f2(A_next,B_next))
    M22 = -Lapp_L2+dt*mu*diags(diff_B_f2(A_next,B_next))
    Giga_matrix_L1 = vstack([M11,M21])
    Giga_matrix_L2 = vstack([M12,M22])
    Giga_matrix = hstack([Giga_matrix_L1,Giga_matrix_L2])
    Giga_matrix_csc = Giga_matrix.tocsc()
    L1_L2_prev = spsolve(Giga_matrix_csc,np.concatenate((second_membre_L1,second_membre_L2),axis=0))
    return (L1_L2_prev[0:len(L1_next)],L1_L2_prev[len(L1_next):])

def calc_CAB1(A_all,B_all,index_mea,D_mea_A,D_mea_B,sigma_eps):
    CAB1 = 0*A_all
    count_mea = 0
    for i in index_mea:
        CAB1[:,i] = -1/sigma_eps**2*(D_mea_A[:,count_mea]-A_all[:,i])
        count_mea = count_mea + 1
    return CAB1

def calc_CAB2(A_all,B_all,index_mea,D_mea_A,D_mea_B,sigma_eps):
    CAB2 = 0*B_all
    count_mea = 0
    for i in index_mea:
        CAB2[:,i] = -1/sigma_eps**2*(D_mea_B[:,count_mea]-B_all[:,i])
        count_mea = count_mea + 1
    return CAB2

def calc_KqqL_1(L1_all_time,L2_all_time,sigma_q,Vol_Omega,tf,Mat_smart_for_kqql):
    # n,m = L1_all_time.shape
    # KqqL1 = np.zeros((n,m))
    # f = lambda t,s,k,L1 : min(t,s)*L1[:,k]
    # for j in range(0,m):
    #     tj = j*tf/(m-1)
    #     int_j = tf/(m-1)*(0.5*(f(tj,0,j,L1_all_time)+f(tj,tf,j,L1_all_time))+sum(f(tj,s*tf/(m-1),j,L1_all_time) for s in range(1,m-1)))
    #     KqqL1[:,j] = Vol_Omega*sigma_q**2*int_j
    # return KqqL1
    return sigma_q**2*np.matmul(L1_all_time,Mat_smart_for_kqql)

def calc_KqqL_2(L1_all_time,L2_all_time,sigma_q,Vol_Omega,tf,Mat_smart_for_kqql):
    # n,m = L2_all_time.shape
    # KqqL2 = np.zeros((n,m))
    # f = lambda t,s,k,L2 : min(t,s)*L2[:,k]
    # for j in range(0,m):
    #     tj = j*tf/(m-1)
    #     int_j = tf/(m-1)*(0.5*(f(tj,0,j,L2_all_time)+f(tj,tf,j,L2_all_time))+sum(f(tj,s*tf/(m-1),j,L2_all_time) for s in range(1,m-1)))
    #     KqqL2[:,j] = Vol_Omega*sigma_q**2*int_j
    # return KqqL2
    return sigma_q**2*np.matmul(L2_all_time,Mat_smart_for_kqql)

def calc_KaaL_1(L1_t0_vec,L2_t0_vec,sigma_a):
    return sigma_a**2*L1_t0_vec

def calc_KaaL_2(L1_t0_vec,L2_t0_vec,sigma_a):
    return sigma_a**2*L2_t0_vec

def calc_proj(lb,ub,eta):
    f = lambda x : 0.5*np.linalg.norm(x-eta,ord = 2)**2
    f_der = lambda x : x-eta
    eta_proj = minimize(f,eta,jac = f_der,bounds = ((lb[0],ub[0]),(lb[1],ub[1])))
    return eta_proj.x

