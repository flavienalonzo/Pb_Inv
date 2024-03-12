# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 11:52:59 2024

@author: Flavien
"""
k1 = 92; k2 = 64; k3 = 1.5; k4 = 18.5; k5 = 0.1

def f1(A,B):
    return k1 - A - h(A,B)

def f2(A,B):
    return k3 * ( k2 - B ) - h(A,B)

def h(A,B):
    return ( k4 * A * B )/( 1 + A + k5 * A**2)

def diff_A_f1(A,B):
    return -1 - diff_A_h(A,B)

def diff_B_f1(A,B):
    return -diff_B_h(A,B)
    
def diff_A_f2(A,B):
    return -diff_A_h(A,B)
        
def diff_B_f2(A,B):
    return -k3 - diff_B_h(A,B)
    
def diff_A_h(A,B):
    return ((k4 * B)*(1 + A + k5 * A**2) - (k4 * A * B )*(1 + 2*k5*A) )/(( 1 + A + k5 * A**2)**2)
        
def diff_B_h(A,B):
    return ( k4 * A )/( 1 + A + k5 * A**2)