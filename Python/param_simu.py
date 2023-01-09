#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 14:22:00 2022

@author: Flavien
"""
import math
import numpy as np

param = {
    "title": "KS simulation",
    "description": "Keller-Segel simulation in 1D",
    "t0": 0.0,
    "tf": 7.1,
    "dt": 0.01,
    "N": 100,
    "diff_u": 5e-5,
    "diff_c": 1e-3,
    "chi": 1e1,
    "dmea": 100,
    "Measures": [],
    "sigma_q": 1e-2,
    "sigma_eps": math.sqrt(1e-2),
    "sigma_a": 1e-1,
    "lb": np.array([np.log(7e-2),np.log(6e-2),np.log(7e-2),np.log(9e-3),np.log(5e-2)]),
    "ub": np.array([np.log(6e-1),np.log(5e-1),np.log(4e-1),np.log(6e-2),np.log(2e-1)]),
    "A": np.array([[-1,1,0,0,0],[0,0,-1,1,0]]),
    "b": np.array([-1e-8,-1e-8]),
    }

param_auto_complet = {
    "h": [], 
    "P": [],
    "x": [], 
    "C_q": [], 
    "C_eps": [], 
    "W_eps": [], 
    "C_a": [], 
    "W_a": [],
    }

sol = {
       "Psi_exact": [],
       "Psi_0": [],
       "Theta_exact": np.array([np.log(0.2),np.log(0.1),np.log(0.1),np.log(0.03),np.log(0.08)]),
       }

estimation_initiale = {
    "coefficients": np.array([np.log(0.2864),np.log(0.0647),np.log(0.3075),np.log(0.0529),np.log(0.1281)]),
    "Psi": [],
    "Lambda": [],
    }

#