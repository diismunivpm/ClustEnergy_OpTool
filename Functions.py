# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 13:32:04 2020

@author: ZEN
"""
import numpy as np
import math
from SingleFamilyHouses import SFH_06, SFH_05, SFH_90, SFH_60
"""
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 

Funzione che permette di scrivere il sistema per i vincoli lineari 
(1) Per edifici (optForm_build())
richiede: 
    samples -> numero timestep simulazione
    Adopt   -> matrice A del sistema discretizzato scritto nella forma SSM
    Xo      -> condizione iniziale del sistema (variabili di stato)
    Edopt   -> matrice E del sistema discretizzato che isola la variabile decisionale 
    Uopt    -> vettore di input 
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
"""

def optForm_build(samples, Adopt, Xo, Bdopt, Uopt, Edopt, p_vd):
    """The optForm_build() function allows to write the building space state system for the linear constraints of the optimization problem. 
        For the calculation the following are required:
            samples -> simulation timestep number
            Adopt -> matrix A of the discretized system written in SSM form
            Xo -> initial condition of the system (state variables)
            Edopt -> matrix E of the discretized system isolating the decision variable 
            Uopt -> input vector """
    
    Nvs       = Adopt.shape[0]
    Pvv       = p_vd   # variable position (system state) to be bound - Tair
    duration  = range(0,samples)
    
    M_opt_max = np.empty([samples,1])
    M         = np.empty([Nvs,samples])
    Z_opt     = np.zeros([samples,samples])
    
    for k in duration:
        if k == 0:
            m          = np.dot(Adopt, Xo) + np.dot(Bdopt, Uopt[:,k]).reshape(Nvs,1)
            Z_opt[0,0] = Edopt[Pvv,0]
        else: 
            Adis  = np.dot(np.linalg.matrix_power(Adopt,(k+1)), Xo)
            b_dis = np.empty([Nvs,k])
            
            for i in range(0,k):
                b_dis[:,i] = np.dot(np.linalg.matrix_power(Adopt,i),np.dot(Bdopt,Uopt[:,(k-i)]).reshape(Nvs,1)).reshape(Nvs,1)[:,0]
           
            for i in range(0,k+1):
                Z_opt[k,(k-i)] = np.dot(np.linalg.matrix_power(Adopt,i),Edopt)[Pvv,0]  
                        
            Bdis = np.sum(b_dis, axis = 1, keepdims = True)
            m = (Adis+Bdis).reshape(Nvs,1)
        M[:,k] = m[:,0]
    
    M_opt_max = M[Pvv,:].reshape(samples,1)
    return Z_opt, M_opt_max, M

def DHW_draw(floor_area,samples):
    Nu = floor_area*2 # m2
    if Nu > 30:
        a = (62*math.log(Nu) - 160)/Nu
    else:
        a = 2
    Vwater_day = a*Nu # water volume per day - m3/day
    Vwater_step = Vwater_day/samples
    
    return(Vwater_step)


def floor_area(SFH,ACH, wall_absorption, radiative_part, convective_part):
    if SFH == 'SFH 2006-today':
        Archetype = SFH_06(ACH, wall_absorption, radiative_part, convective_part)
    if SFH == 'SFH 1991-2005':
        Archetype = SFH_05(ACH, wall_absorption, radiative_part, convective_part)
    if SFH == 'SFH 1976-1990':
        Archetype = SFH_90(ACH, wall_absorption, radiative_part, convective_part)
    if SFH == 'SFH 1946-1960':
        Archetype = SFH_60(ACH, wall_absorption, radiative_part, convective_part)
        
    return(Archetype.floor_area)

def fixed_intGains(floor_area):
    Af = floor_area*2   
    if Af<=120: 
       gains = 7.987 * Af - 0.0353 * Af^2  # Internal gains (UNI 11300-1 pr. 13) 
    else: 
       gains = 450
       
    return(gains)