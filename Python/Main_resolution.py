#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 11:45:06 2022

@author: Flavien
"""
import numpy as np
import module_libraries
import param_simu
from initialisation import simu, solution, estimation
import fun_better_estimation

#Coefs = fun_better_estimation.better_estimation(simu,estimation,solution)

import matplotlib.pyplot as plt


#ax = fig.add_subplot(projection='3d')

#xpoints = simu.x
#ypoints = simu.dt*np.array(range(0,simu.P))
#ypoints = simu.Measures[0]
#xpoints,ypoints = np.meshgrid(xpoints,ypoints)

#ax.plot_surface(xpoints, ypoints, estimation.Psi[0])
#ax.plot_surface(xpoints, ypoints, simu.Measures[2],color='r')

#ax.set_zlim(0.0,1.0)
#plt.plot(xpoints, ypoints[0,:],color='r')
#plt.plot(xpoints, ypoints[1,:],color='b')
#axs[0,0].plot(np.exp(Coefs[:,0]),np.exp(Coefs[:,1]),'k',linestyle='dashed', marker='+')
#plt.plot(range(0,np.size(Coefs,axis=0)),np.exp(Coefs[:,0]),'k+')
#axs[0,0].axvline(np.exp(solution.Theta_exact[0]))
#axs[0,0].axvline(np.exp(simu.lb[0]))
#axs[0,0].axvline(np.exp(simu.ub[0]))
#axs[0,0].axhline(np.exp(solution.Theta_exact[1]))
#axs[0,0].axhline(np.exp(simu.lb[1]))
#axs[0,0].axhline(np.exp(simu.ub[1]))
#axs[0,0].xlim([np.min(np.exp(Coefs[:,0])),np.max(np.exp(Coefs[:,0]))])
#axs[0,0].ylim([np.min(np.exp(Coefs[:,1])),np.max(np.exp(Coefs[:,1]))])
#axs[1,0].plot(np.exp(Coefs[:,2]),np.exp(Coefs[:,3]),'r',linestyle='dashed', marker='+')
#plt.plot(range(0,np.size(Coefs,axis=0)),np.exp(Coefs[:,0]),'k+')
#axs[1,0].axvline(np.exp(solution.Theta_exact[2]))
#axs[1,0].axvline(np.exp(simu.lb[2]))
#axs[1,0].axvline(np.exp(simu.ub[2]))
#axs[1,0].axhline(np.exp(solution.Theta_exact[3]))
#axs[1,0].axhline(np.exp(simu.lb[3]))
#axs[1,0].axhline(np.exp(simu.ub[3]))
#axs[0,0].xlim([np.min(np.exp(Coefs[:,0])),np.max(np.exp(Coefs[:,0]))])
#axs[0,0].ylim([np.min(np.exp(Coefs[:,1])),np.max(np.exp(Coefs[:,1]))])

fig, axs = plt.subplots(2,3)

#axs[0,0].plot(range(0,np.size(Coefs,axis=0)),np.exp(Coefs[:,0]),'k',linestyle='', marker='+')
#plt.plot(range(0,np.size(Coefs,axis=0)),np.exp(Coefs[:,0]),'k+')
axs[0,0].axhline(np.exp(solution.Theta_exact[0]))
axs[0,0].axhline(np.exp(simu.lb[0]))
axs[0,0].axhline(np.exp(simu.ub[0]))
axs[0,0].set_yscale('log')
axs[0,0].grid(axis='y')


#axs[0,1].plot(range(0,np.size(Coefs,axis=0)),np.exp(Coefs[:,1]),'r',linestyle='', marker='+')
#plt.plot(range(0,np.size(Coefs,axis=0)),np.exp(Coefs[:,0]),'k+')
axs[0,1].axhline(np.exp(solution.Theta_exact[1]))
axs[0,1].axhline(np.exp(simu.lb[1]))
axs[0,1].axhline(np.exp(simu.ub[1]))
axs[0,1].set_yscale('log')
axs[0,1].grid(axis='y')

#axs[1,0].plot(range(0,np.size(Coefs,axis=0)),np.exp(Coefs[:,2]),'y',linestyle='', marker='+')
#plt.plot(range(0,np.size(Coefs,axis=0)),np.exp(Coefs[:,0]),'k+')
axs[1,0].axhline(np.exp(solution.Theta_exact[2]))
axs[1,0].axhline(np.exp(simu.lb[2]))
axs[1,0].axhline(np.exp(simu.ub[2]))
axs[1,0].set_yscale('log')
axs[1,0].grid(axis='y')


#axs[1,1].plot(range(0,np.size(Coefs,axis=0)),np.exp(Coefs[:,3]),'g',linestyle='', marker='+')
#plt.plot(range(0,np.size(Coefs,axis=0)),np.exp(Coefs[:,0]),'k+')
axs[1,1].axhline(np.exp(solution.Theta_exact[3]))
axs[1,1].axhline(np.exp(simu.lb[3]))
axs[1,1].axhline(np.exp(simu.ub[3]))
axs[1,1].set_yscale('log')
axs[1,1].grid(axis='y')


#axs[1,2].plot(range(0,np.size(Coefs,axis=0)),np.exp(Coefs[:,4]),'m',linestyle='', marker='+')
#plt.plot(range(0,np.size(Coefs,axis=0)),np.exp(Coefs[:,0]),'k+')
axs[1,2].axhline(np.exp(solution.Theta_exact[4]))
axs[1,2].axhline(np.exp(simu.lb[4]))
axs[1,2].axhline(np.exp(simu.ub[4]))
axs[1,2].set_yscale('log')
axs[1,2].grid(axis='y')


#axs[0,2].plot(range(0,np.size(Coefs,axis=0)),np.linalg.norm(Coefs-np.matmul(np.ones(np.size(Coefs,axis=0)).reshape((np.size(Coefs,axis=0),1)),solution.Theta_exact.reshape((1,5))),axis=1),'b',linestyle='', marker='+')
plt.show()
