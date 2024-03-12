# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 18:47:34 2024

@author: Flavien
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from random import randint
import matplotlib.animation 
from matplotlib.animation import FuncAnimation
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from mpl_toolkits import mplot3d

import pandas as pd

def save_animation(Variable_all_time,nom_var):
    zList = []
    n,p = Variable_all_time.shape
    N = int(np.sqrt(n))
    x = np.linspace(0,1,N); 
    y = np.linspace(0,1,N);
    xx, yy = np.meshgrid(x,y)
    
    for j in range(0,p):
        zList.append(Variable_all_time[:,j].reshape((N,N)))
        
    fig = plt.figure()
    ax = fig.add_subplot(111,projection = '3d')
    
    ax.set_xlim3d([0., 1.])
    ax.set_xlabel('x')
    
    
    ax.set_ylim3d([0., 1.])
    ax.set_ylabel('Y')
    
    
    ax.set_zlim3d([np.min(Variable_all_time), np.max(Variable_all_time)])
    ax.set_title(nom_var)
    
    surf = ax.plot_surface(xx, yy, zList[0],cmap=cm.coolwarm)
    
    #pts = ax.scatter(xList[0], yList[0], zList[0], s=30)
    
    def animate(i, zList, plot):
        # setting data to new scatter coordinates
        ax.clear()
        return ax.plot_surface(xx, yy, zList[i], cmap=cm.coolwarm)
        
        #surf._offsets3d = (xx, yy, zList[i])
        #print(i, zList[i])
        #return surf
    #fig.colorbar(surf, shrink=0.5, aspect=5)
    
    anim = animation.FuncAnimation(fig, animate, fargs=(zList, surf), frames = len(zList), interval=10, repeat=False, blit=False)
        
    #plt.show()
    #anim.save('test.mp4')
    FFwriter = animation.FFMpegWriter(fps=10)
    anim.save(nom_var + '.mp4', writer = FFwriter)
    
    return 