# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 16:01:13 2024

@author: Flavien
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.sparse import diags
from small_functions import *
from big_functions import *

def calc_PF3(A_eta,B_eta,L1_eta,L2_eta,Vol_Omega,tf,Mat_Lap):
    #first_f = lambda i,j,B,L2 : 
    n,m = B_eta.shape
    Mat_Lap_B = np.matmul(Mat_Lap.toarray(),B_eta)
    return Vol_Omega*1/n*tf*1/m*np.array([np.sum(Mat_Lap_B*L2_eta),np.sum(L1_eta*f1(A_eta,B_eta)+L2_eta*f2(A_eta,B_eta))])