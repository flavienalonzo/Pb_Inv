# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 14:20:24 2024

@author: Flavien
"""
from big_functions import *
import numpy as np
import scipy as sp

def calc_PF1(A0_vec,B0_vec,deltax,N,deltat,tf,Coef_mu,Coef_d,Mat_Lap,KqqL1,KqqL2):
    
    Mat_Lap_for_A = sp.sparse.eye(N*N) + deltat/(deltax**2)*Mat_Lap
    Mat_Lap_for_B = sp.sparse.eye(N*N) + Coef_d*deltat/(deltax**2)*Mat_Lap
    Mat_Lap_for_A_csc = Mat_Lap_for_A.tocsc()
    Mat_Lap_for_B_csc = Mat_Lap_for_B.tocsc()
    
    A_next_vec = calc_A_next(A0_vec,B0_vec,deltat,Coef_mu,Mat_Lap_for_A_csc,KqqL1[:,0].reshape(-1))
    B_next_vec = calc_B_next(A0_vec,B0_vec,deltat,Coef_mu,Mat_Lap_for_B_csc,KqqL2[:,0].reshape(-1))

    #time_spent = deltat
    A_all_time = A_next_vec.reshape((N*N,1))
    B_all_time = B_next_vec.reshape((N*N,1))
    
    Ntime = int(tf/deltat)+1
    
    #m=1
    #while (time_spent<tf):
    for m in range(1,Ntime):
        A_prev_vec = A_next_vec
        B_prev_vec = B_next_vec
        A_next_vec = calc_A_next(A_prev_vec,B_prev_vec,deltat,Coef_mu,Mat_Lap_for_A_csc,KqqL1[:,m].reshape(-1))
        B_next_vec = calc_B_next(A_prev_vec,B_prev_vec,deltat,Coef_mu,Mat_Lap_for_B_csc,KqqL2[:,m].reshape(-1))
        A_all_time=np.append(A_all_time,A_next_vec.reshape((N*N,1)),axis=1)
        B_all_time=np.append(B_all_time,B_next_vec.reshape((N*N,1)),axis=1)
        #time_spent = time_spent + deltat
       # m=m+1
        
    return (A_all_time,B_all_time)

