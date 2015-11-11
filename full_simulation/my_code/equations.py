# -*- coding: utf-8 -*-
"""
Created on Wed Nov 04 17:37:37 2015

@author: Kevin Multani (kvmu)

This module serves to define the RK4 solver to solve the system of
ordinary differential equations to find the closed orbit in an
FFAG accelerator. This code has been translated/interpreted 
from C++ to Python by me.

Original code by: Yuya Horita  
"""

def RK4(dt,
        u, v, w, pu, pv, pw,
        dudt, dvdt, dwdt, dpudt, dpvdt, dpwdt,
        init_momentum_magnitude):
    
    ku_1 = dt * dudt(u, v, w, pu, pv, pw)
    kv_1 = dt * dvdt(u , v, pu, pv, pw)
    kw_1 = dt * dwdt(u , v, pu, pv, pw)
    kpu_1 = dt * dpudt(u , v, pu, pv, pw)
    kpv_1 = dt * dpvdt(u , v, pu, pv, pw)
    kpw_1 = dt * dpwdt(u , v, pu, pv, pw)
    
    u_temp = u + ku_1 / 2
    v_temp = v + kv_1 / 2
    w_temp = w + kw_1 / 2
    pu_temp = pu + kpu_1 / 2
    pv_temp = pv + kpv_1 / 2
    pw_temp = pw + kpw_1 / 2
 
    ku_2 = dt * dudt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    kv_2 = dt * dvdt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    kw_2 = dt * dwdt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    kpu_2 = dt * dpudt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    kpv_2 = dt * dpvdt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    kpw_2 = dt * dpwdt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
      
    u_temp = u + ku_2 / 2
    v_temp = v + kv_2 / 2
    w_temp = w + kw_2 / 2
    pu_temp = pu + kpu_2 / 2
    pv_temp = pv + kpv_2 / 2
    pw_temp = pw + kpw_2 / 2

    ku_3 = dt * dudt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    kv_3 = dt * dvdt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    kw_3 = dt * dwdt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    kpu_3 = dt * dpudt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    kpv_3 = dt * dpvdt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    kpw_3 = dt * dpwdt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
 
    u_temp = u + ku_3 / 2
    v_temp = v + kv_3 / 2
    w_temp = w + kw_3 / 2
    pu_temp = pu + kpu_3 / 2
    pv_temp = pv + kpv_3 / 2
    pw_temp = pw + kpw_3 / 2
   
    ku_4 = dt * dudt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    kv_4 = dt * dvdt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    kw_4 = dt * dwdt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    kpu_4 = dt * dpudt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    kpv_4 = dt * dpvdt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    kpw_4 = dt * dpwdt(u_temp, v_temp, w_temp, pu_temp, pv_temp, pw_temp)
    
    du = (ku_1 + 2 * ku_2 + 2 * ku_3 + ku_4) / 6
    dv = (kv_1 + 2 * kv_2 + 2 * kv_3 + kv_4) / 6
    dw = (kw_1 + 2 * kw_2 + 2 * kw_3 + kw_4) / 6
    dpu = (kpu_1 + 2 * kpu_2 + 2 * kpu_3 + kpu_4) / 6
    dpv = (kpv_1 + 2 * kpv_2 + 2 * kpv_3 + kpv_4) / 6
    dpw = (kpw_1 + 2 * kpw_2 + 2 * kpw_3 + kpw_4) / 6
    
    u += du
    v += dv
    w += dw
    pu += dpu
    pv += dpv
    pw += dpw

    check_momentum = (pu**2 + pv**2 + pw**2)**0.5   
    
    pu *= init_momentum_magnitude / check_momentum
    pv *= init_momentum_magnitude / check_momentum
    pw *= init_momentum_magnitude / check_momentum
    
    return u, v, w, pu, pv, pw

    
def lorentz_u(pu, pv, pw, Bu, Bv, Bw):
    return pv * Bw - pw * Bv


def lorentz_v(pu, pv, pw, Bu, Bv, Bw):
    return pw * Bu - pu * Bw
    

def lorentz_w(pu, pv, pw, Bu, Bv, Bw):
    return pu * Bv - pv * Bu
    

def gamma(kinetic_energy, natural_mass = 938.272046):
    return (kinetic_energy + natural_mass) / natural_mass


def drdth(r, th, z, pr, pth, pz, kinetic_energy, natural_mass = 938.272046):
    return r * pr / pth
    

def dzdth(r, th, z, pr, pth, pz, kinetic_energy, natural_mass = 938.272046):
    return r * pz / pth


def dtdth(r, th, z, pr, pth, pz, kinetic_energy, natural_mass = 938.272046,
          mass = 1.6726217817181187e-27):
    return r * gamma(kinetic_energy, natural_mass) * mass / pth
    

def dprdth(r, th, z, pr, pth, pz, Br, Bth, Bz, Q = 1.60217657e-19):
    return Q * r * lorentz_u(pr, pth, pz, Br, Bth, Bz) / pth + pth
  
  
def dpthdth(r, th, z, pr, pth, pz, Br, Bth, Bz, Q = 1.60217657e-19):
    return Q * r * lorentz_v(pr, pth, pz, Br, Bth, Bz) / pth - pr
    
    
def dpzdth(r, th, z, pr, pth, pz, Br, Bth, Bz, Q = 1.60217657e-19):
    return Q * r * lorentz_w(pr, pth, pz, Br, Bth, Bz) / pth 
  
  
  