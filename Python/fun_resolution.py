#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 13:51:33 2022

@author: Flavien
"""
import fun_calc_Lambda
import fun_calc_Psi
import numpy as np
import math
import matplotlib.pyplot as plt

def resolution(simulation,estimation):
    #xpoints = simulation.x
    #ypoints = simulation.dt*np.array(range(0,simulation.P))
    #xpoints,ypoints = np.meshgrid(xpoints,ypoints)
    fun_calc_Psi.calc_Psi(simulation, estimation)
    fun_calc_Lambda.calc_lambda(simulation, estimation)
    Psi1_prev = estimation.Psi[0]
    Psi2_prev = estimation.Psi[1]
    Lambda1_prev = estimation.Lambda[0]
    Lambda2_prev = estimation.Lambda[1]
    norm_prev = math.sqrt(np.linalg.norm(Psi1_prev,'fro')**2+np.linalg.norm(Psi2_prev,'fro')**2+np.linalg.norm(Lambda1_prev,'fro')**2+np.linalg.norm(Lambda2_prev,'fro')**2)
    i = 1
    itermax = 100
    cond = False
    while i<itermax and not cond:
        #fig = plt.figure()
        #ax00 = fig.add_subplot(2, 2, 1, projection='3d')
        #ax01 = fig.add_subplot(2, 2, 2, projection='3d')
        #ax10 = fig.add_subplot(2, 2, 3, projection='3d')
        #ax11 = fig.add_subplot(2, 2, 4, projection='3d')
        fun_calc_Psi.calc_Psi(simulation, estimation)
        fun_calc_Lambda.calc_lambda(simulation, estimation)
        Psi1_new = estimation.Psi[0]
        Psi2_new = estimation.Psi[1]
        Lambda1_new = estimation.Lambda[0]
        Lambda2_new = estimation.Lambda[1]
        #ax00.plot_surface(xpoints, ypoints, Psi1_new)
        #ax01.plot_surface(xpoints, ypoints, Psi2_new)
        #ax10.plot_surface(xpoints, ypoints, Lambda1_new)
        #ax11.plot_surface(xpoints, ypoints, Lambda2_new)
        dif_norm = math.sqrt(np.linalg.norm(Psi1_new-Psi1_prev,'fro')**2+np.linalg.norm(Psi2_new-Psi2_prev,'fro')**2+np.linalg.norm(Lambda1_new-Lambda1_prev,'fro')**2+np.linalg.norm(Lambda2_new-Lambda2_prev,'fro')**2)
        cond = (dif_norm/norm_prev)<5e-3
        if (i>=2): 
            print(i)
            print(dif_norm)
            print(norm_prev)
            print(dif_norm/norm_prev)
        plt.show()
        i+=1
        Psi1_prev = Psi1_new
        Psi2_prev = Psi2_new
        Lambda1_prev = Lambda1_new
        Lambda2_prev = Lambda2_new
        norm_prev = math.sqrt(np.linalg.norm(Psi1_prev,'fro')**2+np.linalg.norm(Psi2_prev,'fro')**2+np.linalg.norm(Lambda1_prev,'fro')**2+np.linalg.norm(Lambda2_prev,'fro')**2)
    #fig = plt.figure()
    #xpoints = simulation.x
    #ypoints = simulation.dt*np.array(range(0,simulation.P))
    #xpoints,ypoints = np.meshgrid(xpoints,ypoints)
    #ax00 = fig.add_subplot(2, 2, 1, projection='3d')
    #ax01 = fig.add_subplot(2, 2, 2, projection='3d')
    #ax10 = fig.add_subplot(2, 2, 3, projection='3d')
    #ax11 = fig.add_subplot(2, 2, 4, projection='3d')
    #ax00.plot_surface(xpoints, ypoints, Psi1_new)
    #ax01.plot_surface(xpoints, ypoints, Psi2_new)
    #ax10.plot_surface(xpoints, ypoints, Lambda1_new)
    #ax11.plot_surface(xpoints, ypoints, Lambda2_new)