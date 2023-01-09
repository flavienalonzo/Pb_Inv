#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:14:42 2022

@author: Flavien
"""
import numpy as np
import scipy.sparse 
import scipy.sparse.linalg
import fun_resolution
import module_libraries
import scipy.optimize

def reaction_KS(simu,estim,x,y):
    rho = np.exp(estim.coefficients[0])
    delta = np.exp(estim.coefficients[1])
    alpha = np.exp(estim.coefficients[2])
    beta = np.exp(estim.coefficients[3])
    gamma = np.exp(estim.coefficients[4])
    return np.stack((rho*x*(1-x)-delta*x,alpha*x-beta*y-gamma*x*y))

def reaction_inv_KS(simu,estim,k,x,y):
    rho = np.exp(estim.coefficients[0])
    delta = np.exp(estim.coefficients[1])
    alpha = np.exp(estim.coefficients[2])
    beta = np.exp(estim.coefficients[3])
    gamma = np.exp(estim.coefficients[4])
    MatE = simu.dt*scipy.sparse.csc_matrix(2*rho*np.diag(estim.Psi[0,k,:])+(delta-rho)*np.eye(simu.N))
    MatF = simu.dt*scipy.sparse.csc_matrix(gamma*np.diag(estim.Psi[1,k,:])-alpha*np.eye(simu.N))
    MatG = scipy.sparse.csc_matrix(np.zeros((simu.N,simu.N)))
    MatH = simu.dt*scipy.sparse.csc_matrix(gamma*np.diag(estim.Psi[0,k,:])+beta*np.eye(simu.N))
    Mat1 = scipy.sparse.bmat([[MatE,MatF]])
    Mat2 = scipy.sparse.bmat([[MatG,MatH]])
    return np.stack((Mat1.dot(np.append(x,y)),Mat2.dot(np.append(x,y))))


def calc_Runge_Kutta_4(simu,estim,x,y):
    (k1u,k1c) = reaction_KS(simu,estim,x,y)
    (k2u,k2c) = reaction_KS(simu,estim,x+simu.dt/2*k1u,y+simu.dt/2*k1c)
    (k3u,k3c) = reaction_KS(simu,estim,x+simu.dt/2*k2u,y+simu.dt/2*k2c)
    (k4u,k4c) = reaction_KS(simu,estim,x+simu.dt*k3u,y+simu.dt*k3c)
    return np.stack((x+simu.dt/6*(k1u+2*k2u+2*k3u+k4u),y+simu.dt/6*(k1c+2*k2c+2*k3c+k4c)))

def calc_Runge_Kutta_4_inv(simu,estim,k,x,y):
    (k1u,k1c) = reaction_inv_KS(simu,estim,k,x,y)
    (k2u,k2c) = reaction_inv_KS(simu,estim,k,x+simu.dt/2*k1u,y+simu.dt/2*k1c)
    (k3u,k3c) = reaction_inv_KS(simu,estim,k,x+simu.dt/2*k2u,y+simu.dt/2*k2c)
    (k4u,k4c) = reaction_inv_KS(simu,estim,k,x+simu.dt*k3u,y+simu.dt*k3c)
    return np.stack((x-simu.dt/6*(k1u+2*k2u+2*k3u+k4u),y-simu.dt/6*(k1c+2*k2c+2*k3c+k4c)))

def calc_c_next(simu,c_prev):
    Mat = simu.diff_c*simu.dt*(2*np.eye(simu.N)-np.diag(np.ones(simu.N-1),-1)-np.diag(np.ones(simu.N-1),1))/pow(simu.h,2)
    Mat2 = scipy.sparse.csc_matrix(Mat+np.eye(simu.N))
    return scipy.sparse.linalg.spsolve(Mat2,c_prev)
    
def calc_u_next(simu,u_prev,c_next):
    Mat = simu.diff_u*simu.dt*(2*np.eye(simu.N)-np.diag(np.ones(simu.N-1),-1)-np.diag(np.ones(simu.N-1),1))/pow(simu.h,2)
    Vec_plus = scipy.sparse.csc_matrix(np.eye(simu.N)-np.diag(np.ones(simu.N-1),1)).dot(c_next)
    Vec_plus[simu.N-1]=0.0
    Vec_moins = scipy.sparse.csc_matrix(np.eye(simu.N)-np.diag(np.ones(simu.N-1),-1)).dot(c_next)
    Vec_moins[0] = 0.0
    Vec_zeros = np.zeros(simu.N)
    C_plus_pos = np.maximum(Vec_plus,Vec_zeros)
    C_plus_neg = np.maximum(-Vec_plus,Vec_zeros)
    C_moins_pos = np.maximum(Vec_moins,Vec_zeros)
    C_moins_neg = np.maximum(-Vec_moins,Vec_zeros)
    Mat2 = simu.diff_u*simu.chi*simu.dt*(-np.diag(C_plus_neg)-np.diag(C_moins_neg)+np.diag(C_plus_pos[0:simu.N-1],1)+np.diag(C_moins_pos[1:simu.N],-1))/pow(simu.h,2)
    Mat3 = scipy.sparse.csc_matrix(Mat-Mat2+np.eye(simu.N))
    return scipy.sparse.linalg.spsolve(Mat3,u_prev)

def measurement(simulation,estimation):
    index_times = simulation.Measures[1]
    U = estimation.Psi[0]
    C = estimation.Psi[1]
    return np.array([U[index_times,:],C[index_times,:]])

def G(simulation,estimation):
    int_1 = np.exp(estimation.coefficients[0])*(estimation.Lambda[0]*estimation.Psi[0]*(1-estimation.Psi[0])).sum()
    int_2 = -np.exp(estimation.coefficients[1])*(estimation.Lambda[0]*estimation.Psi[0]).sum()
    int_3 = np.exp(estimation.coefficients[2])*(estimation.Lambda[1]*estimation.Psi[0]).sum()
    int_4 = -np.exp(estimation.coefficients[3])*(estimation.Lambda[1]*estimation.Psi[1]).sum()
    int_5 = -np.exp(estimation.coefficients[4])*(estimation.Lambda[1]*estimation.Psi[0]*estimation.Psi[1]).sum()
    return simulation.dt*simulation.h*np.array([int_1,int_2,int_3,int_4,int_5])

def e1(simulation,estimation):
    return simulation.dt*np.stack((np.matmul(simulation.C_q,estimation.Lambda[0]),np.matmul(simulation.C_q,estimation.Lambda[1])))

def proj(simulation,coefs):
    Mat_C = np.concatenate((simulation.A,-np.eye(np.size(simulation.lb)),np.eye(np.size(simulation.ub))),axis=0)
    Vec_d = np.concatenate((simulation.b,-simulation.lb,simulation.ub),axis=0)
    pas_besoin_de_proj = all(Mat_C.dot(coefs)-Vec_d<=0)
    if (pas_besoin_de_proj):
        return coefs
    new_x_proj = coefs
    conv = False
    while (not(conv)):
        prior_x_proj = new_x_proj 
        x_proj = prior_x_proj
        for k in range(0,np.size(simulation.b)):
            if ((simulation.A.dot(x_proj)-simulation.b)[k]>0):
                x_proj = x_proj - (simulation.A.dot(x_proj)-simulation.b)[k]/pow(np.linalg.norm(simulation.A[k,:]),2)*simulation.A[k,:]
        new_x_proj = np.maximum(simulation.lb,np.minimum(simulation.ub,x_proj))
        conv = (np.linalg.norm(new_x_proj-prior_x_proj)/np.linalg.norm(prior_x_proj)<1e-8)
    return new_x_proj

#    mu = np.zeros(np.size(simulation.b)+2*np.size(simulation.ub))
#    Mat_C = np.concatenate((simulation.A,-np.eye(np.size(simulation.lb)),np.eye(np.size(simulation.ub))),axis=0)
#    Vec_d = np.concatenate((simulation.b,-simulation.lb,simulation.ub),axis=0)
#    if (all(Mat_C.dot(coefs)-Vec_d<=0)):
#        return coefs
#    prior_mu = mu+1e-1
#    prior_x = coefs
#    rho = 1/np.sum(np.linalg.norm(Mat_C,axis=1))
#    conv = False
#    while (not(conv)):
#        new_x = coefs - prior_mu.dot(Mat_C)
#        new_mu = np.maximum(np.zeros(np.size(mu)),prior_mu+rho*(Mat_C.dot(new_x)-Vec_d))
#        conv = (np.linalg.norm(new_x-prior_x)/np.linalg.norm(prior_x)<1e-8)
#        prior_mu = new_mu
#        prior_x = new_x
#    return new_x
    
    
    
    #fun = lambda x: 0.5*pow(np.linalg.norm(x-coefs),2)
    #der_fun = lambda x : x-coefs
    #bnds = simulation.get_bnds_for_minimize()
    #cons = scipy.optimize.LinearConstraint(-simulation.A, lb =np.array([-np.Inf,-np.Inf]), ub = simulation.b, keep_feasible=True)
    #return scipy.optimize.minimize(fun, coefs, method='SLSQP', jac = der_fun, bounds=bnds,constraints=cons)

def raccourci_G(simulation,estimation,coefficients):
    estimation.coefficients = coefficients
    fun_resolution.resolution(simulation,estimation)
    return G(simulation,estimation)


