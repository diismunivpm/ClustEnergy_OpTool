# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 10:50:41 2020

@author: ZEN
"""

# Occupancy pattern 


import os
import numpy as np
import pandas as pd


def Occupancy_bl(duration,timestep):
    samples  = int(duration/timestep)
    upper_th = []
    lower_th = []
    for t in range(0,samples):
        if 25 <= t <= 26: 
            upper_th.append(26)
            lower_th.append(20)
        else:  
            upper_th.append(26)
            lower_th.append(20)
            
    return(upper_th,lower_th)
    
def Occupancy_bl_floor(duration,timestep):
    samples  = int(duration/timestep)
    upper_th = []
    lower_th = []
    for t in range(0,samples):
        if 25 <= t <= 26: 
            upper_th.append(22)
            lower_th.append(19)
        else:  
            upper_th.append(22)
            lower_th.append(19)
            
    return(upper_th,lower_th)    
    
#def DayNight_DR(timestep, t_in, t_end, variation_up, variation_lo):
#    
#    samples  = int(25/timestep)
#    setpoint = []
#    upper_th = []
#    lower_th = []
#    
#    for t in range(0,samples): 
#        if t < int(6/timestep) or t > int(22/timestep): 
#            setpoint.append(18)
#            upper_th.append(20.5)
#            lower_th.append(18)
#        else: 
#            setpoint.append(20)
#            upper_th.append(20.5)
#            lower_th.append(20)
#    
#    setpoint_DR = []
#    upper_th_DR = []
#    lower_th_DR = []        
#    
#    for t in range(0,samples):        
#        if t_in <= t <=t_end:
#            setpoint_DR.append(setpoint[t])
#            upper_th_DR.append(upper_th[t]+ variation_up)
#            lower_th_DR.append(lower_th[t]+ variation_lo)
#            
#        else: 
#            setpoint_DR.append(setpoint[t])
#            upper_th_DR.append(upper_th[t])
#            lower_th_DR.append(lower_th[t])
#            
#    return(setpoint_DR,upper_th_DR,lower_th_DR) 
#    
#    
#    
#def MoringEvening(timestep):
#    samples  = int(25/timestep)
#    setpoint = []
#    upper_th = []
#    lower_th = []
#    for t in range(0,samples): 
#        if int(6/timestep) <= t <= int(8/timestep) or int(18/timestep) <= t <= int(22/timestep): 
#            setpoint.append(20)
#            upper_th.append(20.5)
#            lower_th.append(20)
#        else: 
#            setpoint.append(18)
#            upper_th.append(20.5)
#            lower_th.append(18)
#    return(setpoint,upper_th,lower_th)    
#    
#    
#def MoringEvening_DR(timestep, t_in, t_end, variation_up, variation_lo):
#    samples  = int(25/timestep)
#    setpoint = []
#    upper_th = []
#    lower_th = []
#    for t in range(0,samples): 
#        if int(6/timestep) <= t <= int(8/timestep) or int(18/timestep) <= t <= int(22/timestep): 
#            setpoint.append(20)
#            upper_th.append(20.5)
#            lower_th.append(20)
#        else: 
#            setpoint.append(18)
#            upper_th.append(20.5)
#            lower_th.append(18)
#    
#    setpoint_DR = []
#    upper_th_DR = []
#    lower_th_DR = []        
#    
#    for t in range(0,samples):        
#        if t_in <= t <=t_end:
#            setpoint_DR.append(setpoint[t])
#            upper_th_DR.append(upper_th[t]+ variation_up)
#            lower_th_DR.append(lower_th[t]+ variation_lo)
#            
#        else: 
#            setpoint_DR.append(setpoint[t])
#            upper_th_DR.append(upper_th[t])
#            lower_th_DR.append(lower_th[t])
#            
#    return(setpoint_DR,upper_th_DR,lower_th_DR)     