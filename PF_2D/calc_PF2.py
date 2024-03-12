# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 16:30:25 2024

@author: Flavien
"""
import numpy as np
import scipy as sp
from scipy.sparse import diags
from big_functions import *

def calc_PF2(A_all_time,B_all_time,index_mea,D_mea_A,D_mea_B,sigma_eps,N,deltax,deltat,tf,eta,Mat_Lap):    
    L1_tf = np.zeros((N,N))
    L2_tf = np.zeros((N,N))
    Coef_d,Coef_mu = eta

    L1_all_time = L1_tf.reshape((N*N,1))
    L2_all_time = L2_tf.reshape((N*N,1))

    L1_tf_vec = L1_tf.reshape(-1)
    L2_tf_vec = L2_tf.reshape(-1)
    
    ##La matrice du laplacien plus compliquée pour L1 et L2 :
    Mat_Lap_for_L1 = -sp.sparse.eye(N*N) - deltat/(deltax**2)*Mat_Lap 
    Mat_Lap_for_L2 = -sp.sparse.eye(N*N) - Coef_d*deltat/(deltax**2)*Mat_Lap
    # if (sp.sparse.linalg.eigs(Mat_Lap_for_L1,k=1,which='SR')[0]*sp.sparse.linalg.eigs(Mat_Lap_for_L1,k=1,which='LR')[0]):
    #     print('la matrice L1 nest pas definie positive')
    #     print(sp.sparse.linalg.eigs(Mat_Lap_for_L1,k=1,which='SR')[0],sp.sparse.linalg.eigs(Mat_Lap_for_L1,k=1,which='LR')[0])
    # if (sp.sparse.linalg.eigs(Mat_Lap_for_L2,k=1,which='SR')[0]*sp.sparse.linalg.eigs(Mat_Lap_for_L2,k=1,which='LR')[0]):
    #     print('la matrice L2 nest pas definie positive')
    #     print(sp.sparse.linalg.eigs(Mat_Lap_for_L2,k=1,which='SR')[0],sp.sparse.linalg.eigs(Mat_Lap_for_L2,k=1,which='LR')[0])
    Mat_Lap_for_L1_csc = Mat_Lap_for_L1.tocsc()
    Mat_Lap_for_L2_csc = Mat_Lap_for_L2.tocsc()
    # Mat_Lap_for_L1_array = Mat_Lap_for_L1.toarray()
    # Mat_Lap_for_L2_array = Mat_Lap_for_L2.toarray()

    CAB1 = calc_CAB1(A_all_time,B_all_time,index_mea,D_mea_A,D_mea_B,sigma_eps)
    CAB2 = calc_CAB2(A_all_time,B_all_time,index_mea,D_mea_A,D_mea_B,sigma_eps)

    Ntime = int(tf/deltat)+1

    # L1_prev = calc_L1_prev(L1_tf_vec,L2_tf_vec,A_all_time[:,-1],B_all_time[:,-1],deltat,Coef_mu,Mat_Lap_for_L1_csc,CAB1[:,-1])
    # L2_prev = calc_L2_prev(L1_tf_vec,L2_tf_vec,A_all_time[:,-1],B_all_time[:,-1],deltat,Coef_mu,Mat_Lap_for_L2_csc,CAB2[:,-1])

    # #(L1_prev,L2_prev) = calc_L1_L2_prev(L1_tf_vec,L2_tf_vec,A_all_time[:,-1],B_all_time[:,-1],deltat,Coef_mu,Mat_Lap_for_L1,Mat_Lap_for_L2,CAB1[:,-1],CAB2[:,-1])

    # L1_all_time=np.append(L1_all_time,L1_prev.reshape((N*N,1)),axis=1)
    # L2_all_time=np.append(L2_all_time,L2_prev.reshape((N*N,1)),axis=1)
    L1_next = L1_tf_vec
    L2_next = L2_tf_vec

    for m in range(Ntime-1,0,-1):
        L1_prev = calc_L1_prev(L1_next,L2_next,A_all_time[:,m],B_all_time[:,m],deltat,Coef_mu,Mat_Lap_for_L1_csc,CAB1[:,m])
        L2_prev = calc_L2_prev(L1_next,L2_next,A_all_time[:,m],B_all_time[:,m],deltat,Coef_mu,Mat_Lap_for_L2_csc,CAB2[:,m])
        #(L1_prev,L2_prev) = calc_L1_L2_prev(L1_next,L2_next,A_all_time[:,m],B_all_time[:,m],deltat,Coef_mu,Mat_Lap_for_L1,Mat_Lap_for_L2,CAB1[:,m],CAB2[:,m])
        L1_all_time=np.append(L1_all_time,L1_prev.reshape((N*N,1)),axis=1)
        L2_all_time=np.append(L2_all_time,L2_prev.reshape((N*N,1)),axis=1)
        L1_next = L1_prev
        L2_next = L2_prev
        if (np.isnan(L1_prev).any() or np.isnan(L1_prev).any()):
            print('Il y a des nan dans le calcul de (PF2)')
            break

    return (np.fliplr(L1_all_time),np.fliplr(L2_all_time))