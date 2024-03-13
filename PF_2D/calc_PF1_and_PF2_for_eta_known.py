# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 15:05:23 2024

@author: Flavien
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.sparse import diags
from small_functions import *
from big_functions import *
from calc_PF1 import *
from calc_PF2 import *

def calc_PF1_PF2(eta,N,deltax,deltat,tf,Vol_Omega,Mat_smart_for_kqql,index_mea,D_mea_A,D_mea_B,sigma_eps,sigma_q,sigma_a,Mat_Lap,A0,B0):
    Coef_d,Coef_mu = eta
    Ntime = int(tf/deltat)+1
    x = np.linspace(0,1,N); 
    y = np.linspace(0,1,N);
    xx, yy = np.meshgrid(x,y)
    taux_conv_1 = 0.05      #From the previous step we need less than 5% of changes (relative)
    taux_conv_2 = 0.001      #From the previous step we need less than 5% of changes (absolute)
    #First (PF1)
    KqqL1 = np.zeros((N*N,Ntime))
    KqqL2 = np.zeros((N*N,Ntime))
    (A_all_time,B_all_time) = calc_PF1(A0,B0,deltax,N,deltat,tf,Coef_mu,Coef_d,Mat_Lap,KqqL1,KqqL2)
    
    #Initialisation de Lambda:
    L1_all_time = np.zeros((N*N,Ntime))
    L2_all_time = np.zeros((N*N,Ntime))

    #boucle itérative 
    for iter_count in range(1,100):
        print('\t','\t','boucle numéro ',iter_count)
        #To check if convergence
        A_prev = A_all_time
        B_prev = B_all_time
        L1_prev = L1_all_time
        L2_prev = L2_all_time
        
        #(PF2)
        (L1_all_time,L2_all_time) = calc_PF2(A_all_time,B_all_time,index_mea,D_mea_A,D_mea_B,sigma_eps,N,deltax,deltat,tf,eta,Mat_Lap)
#        print('\t','\t','\t','(PF2) done')
        #Actualisation des KqqL et KaaL0
        KqqL1 = calc_KqqL_1(L1_all_time,L2_all_time,sigma_q,Vol_Omega,tf,Mat_smart_for_kqql)
        KqqL2 = calc_KqqL_2(L1_all_time,L2_all_time,sigma_q,Vol_Omega,tf,Mat_smart_for_kqql)
        KaaL0_1 = calc_KaaL_1(L1_all_time[:,0].reshape(-1),L2_all_time[:,0].reshape(-1),sigma_a)
        KaaL0_2 = calc_KaaL_2(L1_all_time[:,0].reshape(-1),L2_all_time[:,0].reshape(-1),sigma_a)
#        print('\t','\t','\t','KqqL and KaaL0 done')
        #(PF1)
        #                                  
        (A_all_time,B_all_time) = calc_PF1(A0+KaaL0_1,B0+KaaL0_2,deltax,N,deltat,tf,Coef_mu,Coef_d,Mat_Lap,KqqL1,KqqL2)
#        print('\t','\t','\t','(PF1) done')
        #Convergence ?
        conv_A_1 = np.linalg.norm(A_all_time-A_prev,ord='fro') < taux_conv_1*np.linalg.norm(A_prev,ord='fro')
        conv_A_2 = np.linalg.norm(A_all_time-A_prev,ord='fro') < taux_conv_2
        conv_B_1 = np.linalg.norm(B_all_time-B_prev,ord='fro') < taux_conv_1*np.linalg.norm(B_prev,ord='fro')
        conv_B_2 = np.linalg.norm(B_all_time-B_prev,ord='fro') < taux_conv_2
        conv_L1_1 = np.linalg.norm(L1_all_time-L1_prev,ord='fro') < taux_conv_1*np.linalg.norm(L1_prev,ord='fro')
        conv_L1_2 = np.linalg.norm(L1_all_time-L1_prev,ord='fro') < taux_conv_2
        conv_L2_1 = np.linalg.norm(L2_all_time-L2_prev,ord='fro') < taux_conv_1*np.linalg.norm(L2_prev,ord='fro')
        conv_L2_2 = np.linalg.norm(L2_all_time-L2_prev,ord='fro') < taux_conv_2
        if(conv_A_1 and conv_A_2 and conv_B_1 and conv_B_2 and conv_L1_1 and conv_L1_2 and conv_L2_1 and conv_L2_2):
            # print('\t','\t','\t','\t','conv A',np.linalg.norm(A_all_time-A_prev,ord='fro'),np.linalg.norm(A_prev,ord='fro'))
            # print('\t','\t','\t','\t','conv B',np.linalg.norm(B_all_time-B_prev,ord='fro'),np.linalg.norm(B_prev,ord='fro'))
            # print('\t','\t','\t','\t','conv L1',np.linalg.norm(L1_all_time-L1_prev,ord='fro'),np.linalg.norm(L1_prev,ord='fro'))
            # print('\t','\t','\t','\t','conv L2',np.linalg.norm(L2_all_time-L2_prev,ord='fro'),np.linalg.norm(L2_prev,ord='fro'))
            # fig = plt.figure(figsize=plt.figaspect(0.5))
            # ax = fig.add_subplot(2, 2, 1, projection='3d')
            # surf = ax.plot_surface(xx, yy, A_all_time[:,-1].reshape((N,N)))
            # ax = fig.add_subplot(2, 2, 2, projection='3d')
            # surf = ax.plot_surface(xx, yy, B_all_time[:,-1].reshape((N,N)))
            # ax = fig.add_subplot(2, 2, 3, projection='3d')
            # surf = ax.plot_surface(xx, yy, L1_all_time[:,0].reshape((N,N)))
            # ax = fig.add_subplot(2, 2, 4, projection='3d')
            # surf = ax.plot_surface(xx, yy, L2_all_time[:,0].reshape((N,N)))
            # plt.show()
            return (A_all_time,B_all_time,L1_all_time,L2_all_time)
        # print('\t','\t','\t','\t','conv A',np.linalg.norm(A_all_time-A_prev,ord='fro'),np.linalg.norm(A_prev,ord='fro'))
        # print('\t','\t','\t','\t','conv B',np.linalg.norm(B_all_time-B_prev,ord='fro'),np.linalg.norm(B_prev,ord='fro'))
        # print('\t','\t','\t','\t','conv L1',np.linalg.norm(L1_all_time-L1_prev,ord='fro'),np.linalg.norm(L1_prev,ord='fro'))
        # print('\t','\t','\t','\t','conv L2',np.linalg.norm(L2_all_time-L2_prev,ord='fro'),np.linalg.norm(L2_prev,ord='fro'))
    print('\t','\t','\t','Convergence has not been met: ',conv_A_1,conv_A_2,conv_B_1,conv_B_2,conv_L1_1,conv_L1_2,conv_L2_1,conv_L2_2)
    return (A_all_time,B_all_time,L1_all_time,L2_all_time)