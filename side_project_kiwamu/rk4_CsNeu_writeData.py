# -*- coding: utf-8 -*-
"""
Created on Tue Sep 01 00:22:11 2015

@author: Kevin Multani 

This python script is to solve a system of ordinary differential equations.
These set of differential equations looks at a Cs-135 muon capture -> 
neutron production experiment. This script uses a generic fourth order 
Runge-Kutta solver to find the solution of the system and writes the values
a file.

Last updated: 2/9/2015-13:49
"""

from math import log

## Initialize variables

# gamma: muon beam deliver rate [mol/sec]
gamma = (3600*24*365)**(-1)

# conversion: conversion factor years -> seconds
conversion = gamma**(-1)

# lambda_cs135: mean life-time of Cesium-135 [1/sec]
lambda_cs135 = log(2)/(2.3*10**6*conversion) 

# lambda_cs133: mean life-time of Xenon-133 [1/sec]
lambda_xe133 = log(2)/(5.2475/365*conversion)

# lambda_cs135: mean life-time of Xenon-135 [1/sec]
lambda_xe135 = log(2)/(9.14/(24*365)*conversion)


def rk4_solve(CS135, CS133, XE135, XE134, XE133, XE132, BA135, fCS135, fCS133, fXE135, fXE134, fXE133, fXE132, fBA135, dt):

    """ 
    
    This function solves the system of equations describing the dynamics
    of neutron-production via muon transmutation. It takes in the
    initial values of the problem and returns the value after 1 step
    of the Runge-Kutta algorithm.
    
    """
      
    CS135_1 = fCS135(CS135, CS133, XE135, XE134, XE133, XE132, BA135)*dt
    CS133_1 = fCS133(CS135, CS133, XE135, XE134, XE133, XE132, BA135)*dt
    XE135_1 = fXE135(CS135, CS133, XE135, XE134, XE133, XE132, BA135)*dt
    XE134_1 = fXE134(CS135, CS133, XE135, XE134, XE133, XE132, BA135)*dt
    XE133_1 = fXE133(CS135, CS133, XE135, XE134, XE133, XE132, BA135)*dt
    XE132_1 = fXE132(CS135, CS133, XE135, XE134, XE133, XE132, BA135)*dt
    BA135_1 = fBA135(CS135, CS133, XE135, XE134, XE133, XE132, BA135)*dt

    CS135_k = CS135 + CS135_1*0.5   
    CS133_k = CS133 + CS133_1*0.5
    XE135_k = XE135 + XE135_1*0.5
    XE134_k = XE134 + XE134_1*0.5
    XE133_k = XE133 + XE133_1*0.5
    XE132_k = XE132 + XE132_1*0.5
    BA135_k = BA135 + BA135_1*0.5
 
    CS135_2 = fCS135(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    CS133_2 = fCS133(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    XE135_2 = fXE135(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    XE134_2 = fXE134(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    XE133_2 = fXE133(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    XE132_2 = fXE132(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    BA135_2 = fBA135(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    
    CS135_k = CS135 + CS135_2*0.5
    CS133_k = CS133 + CS133_2*0.5
    XE135_k = XE135 + XE135_2*0.5
    XE134_k = XE134 + XE134_2*0.5
    XE133_k = XE133 + XE133_2*0.5
    XE132_k = XE132 + XE132_2*0.5
    BA135_k = BA135 + BA135_2*0.5
     
    CS135_3 = fCS135(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    CS133_3 = fCS133(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    XE135_3 = fXE135(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    XE134_3 = fXE134(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    XE133_3 = fXE133(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    XE132_3 = fXE132(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    BA135_3 = fBA135(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
     
    CS135_k = CS135 + CS135_3
    CS133_k = CS133 + CS133_3
    XE135_k = XE135 + XE135_3
    XE134_k = XE134 + XE134_3
    XE133_k = XE133 + XE133_3
    XE132_k = XE132 + XE132_3
    BA135_k = BA135 + BA135_3

    CS135_4 = fCS135(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    CS133_4 = fCS133(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    XE135_4 = fXE135(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    XE134_4 = fXE134(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    XE133_4 = fXE133(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    XE132_4 = fXE132(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    BA135_4 = fBA135(CS135_k, CS133_k, XE135_k, XE134_k, XE133_k, XE132_k, BA135_k)*dt
    
    CS135 = CS135 + (CS135_1 + 2*(CS135_2 + CS135_3) + CS135_4)/6
    CS133 = CS133 + (CS133_1 + 2*(CS133_2 + CS133_3) + CS133_4)/6
    XE135 = XE135 + (XE135_1 + 2*(XE135_2 + XE135_3) + XE135_4)/6
    XE134 = XE134 + (XE134_1 + 2*(XE134_2 + XE134_3) + XE134_4)/6
    XE133 = XE133 + (XE133_1 + 2*(XE133_2 + XE133_3) + XE133_4)/6
    XE132 = XE132 + (XE132_1 + 2*(XE132_2 + XE132_3) + XE132_4)/6
    BA135 = BA135 + (BA135_1 + 2*(BA135_2 + BA135_3) + BA135_4)/6
     
    return CS135, CS133, XE135, XE134, XE133, XE132, BA135


# Functions that define each row of the ODE matrix    
def cs135(CS135, CS133, XE135, XE134, XE133, XE132, BA135):
    return -(gamma+lambda_cs135)*CS135 + lambda_xe135*XE135

def cs133(CS135, CS133, XE135, XE134, XE133, XE132, BA135):
    return lambda_xe133*XE133

def xe135(CS135, CS133, XE135, XE134, XE133, XE132, BA135):
    return gamma*0.1*CS135 - lambda_xe135*XE135

def xe134(CS135, CS133, XE135, XE134, XE133, XE132, BA135):
    return gamma*0.5*CS135

def xe133(CS135, CS133, XE135, XE134, XE133, XE132, BA135):
    return gamma*0.3*CS135 - lambda_xe133*XE133
    
def xe132(CS135, CS133, XE135, XE134, XE133, XE132, BA135):
    return gamma*0.1*CS135

def ba135(CS135, CS133, XE135, XE134, XE133, XE132, BA135):
    return lambda_cs135*CS135

# Initial conditions     
CS135, CS133, XE135, XE134, XE133, XE132, BA135, dt = 1, 0, 0, 0, 0, 0, 0, 60

day = 0
neutron_expected = 0*0.1 + 1*0.5 + 2*0.3 + 3*0.1

fh = open('rk4solution.txt', 'w+')
neu_fh = open('neutron_production_rate.txt', 'w+')

# Apply Runge-Kutta and take 144000000 steps, which equals to 10,000 days since
# each time step is 60 seconds (1 minute).
for i in xrange(14400000):
    CS135, CS133, XE135, XE134, XE133, XE132, BA135 = rk4_solve(CS135, CS133, XE135, XE134, XE133, XE132, BA135, cs135, cs133, xe135, xe134, xe133, xe132, ba135, dt)
    
    if i % 1440 == 0:
        day+=1
        fh.write(str(day)+'\t'+str(CS135)+'\t'+str(CS133)+'\t'+str(XE135)
                 +'\t'+str(XE134)+'\t'+str(XE133)+'\t'+str(XE132)
                 +'\t'+str(BA135)+'\n')
        neu_fh.write(str(day) +'\t'+ str(neutron_expected*gamma*CS135)+'\n')

neu_fh.close()
fh.close()

     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     