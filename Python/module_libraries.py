#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 13:49:56 2022

@author: Flavien
"""

class Simulation:
    
    def __init__(self, title, description, t0, tf, dt, N, diff_u, diff_c, chi, dmea, Measures, sigma_q, sigma_eps, sigma_a, lb, ub, A, b, h, P, x, C_q, C_eps, W_eps, C_a, W_a):
        self.title = title
        self.description = description
        self.t0 = t0
        self.tf = tf
        self.dt = dt
        self.N = N
        self.diff_u = diff_u
        self.diff_c = diff_c
        self.chi = chi
        self.dmea = dmea
        self.Measures = Measures
        self.sigma_q = sigma_q
        self.sigma_eps = sigma_eps
        self.sigma_a = sigma_a
        self.lb = lb
        self.ub = ub
        self.A = A
        self.b = b
        self.h = h
        self.P = P
        self.x = x
        self.C_q = C_q
        self.C_eps = C_eps
        self.W_eps = W_eps
        self.C_a = C_a 
        self.W_a = W_a
    
    def get_bnds_for_minimize(self):
        bnds = []
        for i in range(0,len(self.lb)):
            bnds.append((self.lb[i],self.ub[i]))
        return tuple(bnds)

    
class Solution_exacte(Simulation):
    def __init__(self, title, description, t0, tf, dt, N, diff_u, diff_c, chi, dmea, Measures, sigma_q, sigma_eps, sigma_a, lb, ub, A, b, h, P, x, C_q, C_eps, W_eps, C_a, W_a,Psi_exact,Psi_0,Theta_exact):
        Simulation.__init__(self, title, description, t0, tf, dt, N, diff_u, diff_c, chi, dmea, Measures, sigma_q, sigma_eps, sigma_a, lb, ub, A, b, h, P, x, C_q, C_eps, W_eps, C_a, W_a)
        self.Psi_exact = Psi_exact
        self.Psi_0 = Psi_0
        self.Theta_exact= Theta_exact
        
class Estimation:
    def __init__(self,coefficients,Psi,Lambda):
        self.coefficients = coefficients
        self.Psi = Psi
        self.Lambda = Lambda
    