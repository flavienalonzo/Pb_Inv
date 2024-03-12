# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 17:32:20 2024

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
from calc_PF3 import *
from calc_PF1_and_PF2_for_eta_known import *

def Methode_NP_proj(lb,ub,eta_0,N,deltax,deltat,tf,Vol_Omega,Mat_Lap_for_A,Mat_Lap_for_B,Mat_smart_for_kqql,index_mea,D_mea_A,D_mea_B,sigma_eps,sigma_q,sigma_a,Mat_Lap,A0_vec,B0_vec):
    #First A,B,L1,L2
    f = lambda eta : calc_PF1_PF2(eta,N,deltax,deltat,tf,Vol_Omega,Mat_smart_for_kqql,index_mea,D_mea_A,D_mea_B,sigma_eps,sigma_q,sigma_a,Mat_Lap,A0_vec,B0_vec)
    h = 0.0001
    n_eta = len(eta_0)
    eta = eta_0
    all_eta = eta_0.reshape((n_eta,1))
    for iter_meth in range(1,50):
        print('Boucle algo Newton numéro : ',iter_meth)
        print('\t','calcul de F(eta)')
        (A_eta,B_eta,L1_eta,L2_eta) = f(eta)
        err_eta = calc_PF3(A_eta,B_eta,L1_eta,L2_eta,Vol_Omega,tf,Mat_Lap)

        print('\t','approximation de la jacobienne')
        vec_zeros = np.zeros(n_eta)
        Jac = np.zeros((n_eta,n_eta))
        for i in range(0,n_eta):
            e_i = np.zeros(n_eta)
            e_i[i]= 1
            eta_i = eta +h*e_i
            (A_eta_i,B_eta_i,L1_eta_i,L2_eta_i) = f(eta_i)
            err_eta_i = calc_PF3(A_eta_i,B_eta_i,L1_eta_i,L2_eta_i,Vol_Omega,tf,Mat_Lap)
            for j in range(0,n_eta):
                e_j = np.zeros(n_eta)
                e_j[j] = 1 
                Jac[j,i] = e_j.dot(err_eta_i-err_eta)/h
                
        print('\t','\t','Jacobienne : ',Jac)
        second_membre = -err_eta
        
        print('\t','\t','Second membre : ',second_membre)
        new_eta = eta + np.linalg.solve(Jac,second_membre)
        
        print('\t','\t','Nouveaux coefficients sans projection : ',new_eta)
        new_eta_proj = calc_proj(lb,ub,new_eta)
        
        all_eta=np.append(all_eta,new_eta_proj.reshape((n_eta,1)),axis=1)
        if (np.linalg.norm(new_eta_proj-eta,ord=2)<0.05*np.linalg.norm(eta,ord=2) and np.linalg.norm(new_eta_proj-eta,ord=2)<0.01):
            print('conv',np.linalg.norm(new_eta_proj-eta,ord=2),np.linalg.norm(eta,ord=2))
            return all_eta,A_eta,B_eta,L1_eta,L2_eta
        print('\t','\t','Après la projection : ',new_eta_proj)
        eta = new_eta_proj
    print('\t','La méthode de Newton-Raphson n a pas convergé')
    return eta,A_eta,B_eta,L1_eta,L2_eta