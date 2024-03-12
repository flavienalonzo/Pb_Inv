# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 19:31:22 2024

@author: Flavien
"""

import numpy as np
import math
import matplotlib.pyplot as plt

def graphique(A_exact,B_exact,A_eta,B_eta,L1_eta,L2_eta,all_eta,A_final,B_final,L1_final,L2_final,tf):

    n,m = B_eta.shape
    ## 1er graphique 4 sous-figures :
    fig, ax = plt.subplots(2,2,figsize=(10,8))
        ##Le tracé de ||A-A_ex||_L2(Omega) au cours du temps pour A le premier test et le A final après calibration
    ax[0,0].plot(np.linspace(0,tf,m),np.sqrt(1/n*np.sum((A_eta-A_exact)**2,axis=0)),c='blue',label='without calibration')
    ax[0,0].plot(np.linspace(0,tf,m),np.sqrt(1/n*np.sum((A_final-A_exact)**2,axis=0)),c='orangered',label='after calibration')
    ax[0,0].set_title('$||A-A_{ex}||_{L^2(\Omega)}$', size=20, fontweight='bold')
    ax[0,0].set_yscale('log')
    ax[0,0].set_xlim([0, tf])
        ##Le tracé de ||B-B_ex||_L2(Omega) au cours du temps pour B le premier test et le B final après calibration
    ax[0,1].plot(np.linspace(0,tf,m),np.sqrt(1/n*np.sum((B_eta-B_exact)**2,axis=0)),c='blue',label='without calibration')
    ax[0,1].plot(np.linspace(0,tf,m),np.sqrt(1/n*np.sum((B_final-B_exact)**2,axis=0)),c='orangered',label='after calibration')
    ax[0,1].set_title('$||B-B_{ex}||_{L^2(\Omega)}$', size=20, fontweight='bold')
    ax[0,1].set_yscale('log') 
    ax[0,1].set_xlim([0, tf])
    ax[0,1].legend()
        
        ##Le tracé de ||L1||_L2(Omega) au cours du temps pour L1 le premier test et le L1 final après calibration
    ax[1,0].plot(np.linspace(0,tf,m),np.sqrt(1/n*np.sum((L1_eta)**2,axis=0)),c='blue',label='without calibration')
    ax[1,0].plot(np.linspace(0,tf,m),np.sqrt(1/n*np.sum((L1_final)**2,axis=0)),c='orangered',label='after calibration')
    ax[1,0].set_title('$||\lambda_1||_{L^2(\Omega)}$', size=20, fontweight='bold')
    ax[1,0].set_yscale('log')
    ax[1,0].set_xlim([0, tf])
        
        ##Le tracé de ||L2||_L2(Omega) au cours du temps pour L2 le premier test et le L2 final après calibration
    ax[1,1].plot(np.linspace(0,tf,m),np.sqrt(1/n*np.sum((L2_eta)**2,axis=0)),c='blue',label='without calibration')
    ax[1,1].plot(np.linspace(0,tf,m),np.sqrt(1/n*np.sum((L2_final)**2,axis=0)),c='orangered',label='after calibration')
    ax[1,1].set_title('$||\lambda_2||_{L^2(\Omega)}$', size=20, fontweight='bold')
    ax[1,1].set_yscale('log')
    ax[1,1].set_xlim([0, tf])
        
    plt.savefig('analyse_2D.png', dpi=120)
    plt.show()
        
        
        
        
    # ## 2nd graphique 2 sous-graphiques :
    fig, ax = plt.subplots(1,2,figsize=(10,5))
        ## Les itérations du coefficient d durant le process de l'algo de Newton
    ax[0].scatter(range(all_eta.shape[1]),all_eta[0,:],marker='*',c='blue',label='$(d^k)_{k\geq 0}$')
    ax[0].hlines(y=20,xmin=0,xmax=all_eta.shape[1]-1,linestyles='solid',colors='blue',label='$d_{exact}$')
    ax[0].hlines(y=15,xmin=0,xmax=all_eta.shape[1]-1,linestyles='dashed',colors='blue',label='lower bound')
    ax[0].hlines(y=24.95,xmin=0,xmax=all_eta.shape[1],linestyles='dashed',colors='blue',label='upper bound')
    ax[0].set_title('Estimation of d', size=20, fontweight='bold')
    ax[0].grid(axis='y')
    ax[0].set_yticks(range(15,26))
    ax[0].set_xlim([0, all_eta.shape[1]-0.95])
    ax[0].set_ylim([15, 25])
    ax[0].legend()
        ## Les itérations du coefficient mu durant le process de l'algo de Newton
    ax[1].scatter(range(all_eta.shape[1]),all_eta[1,:],marker='*',c='orangered',label='$(\mu^k)_{k\geq 0}$')
    ax[1].hlines(y=30,xmin=0,xmax=all_eta.shape[1]-1,linestyles='solid',colors='orangered',label='$d_{exact}$')
    ax[1].hlines(y=25,xmin=0,xmax=all_eta.shape[1]-1,linestyles='dashed',colors='orangered',label='lower bound')
    ax[1].hlines(y=34.95,xmin=0,xmax=all_eta.shape[1],linestyles='dashed',colors='orangered',label='upper bound')
    ax[1].set_title('Estimation of $\mu$', size=20, fontweight='bold')
    ax[1].grid(axis='y')
    ax[1].set_yticks(range(25,36))
    ax[1].set_xlim([0, all_eta.shape[1]-0.95])
    ax[1].set_ylim([25, 35])
    ax[1].legend()
    plt.savefig('calibration_2D.png', dpi=120)
    plt.show()

    return