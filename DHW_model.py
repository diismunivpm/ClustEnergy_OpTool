# -*- coding: utf-8 -*-
"""
ClustEnergy OpTool
Released on January, 2024
@author: DIISM UNIVPM
"""
import numpy as np
from scipy import signal

class DHW_tank(object):
    """
    Insulated DHW tank.
    """
    def __init__(self, water_volume, surf_area, wall_U):
        """Definition of building properties. """
    
        self.density = 997          # kg/m3
        self.specific_heat = 4186   # J/(kg/K)
        self.volume = water_volume  # water volume of the tank - m3
        self.A_tank = surf_area
        self.U_tank = wall_U
    """RC properties."""
    
    @property
    def Cwater(self):
        return (self.density*self.specific_heat*self.volume/3600) # Water thermal capacitance (Wh/K)
    
    @property
    def Rwa(self):       
        return (1/(self.U_tank*self.A_tank))

    @property
    def Kwa(self):       
        return (self.U_tank*self.A_tank)            # Thermal conductance water-ambient (W/K)
    


    def stateSpace_DHW(self) :
        """State space model formulation for air heating system. """
        # Thermal node of the water (Twater)
        A11 = -(self.Kwa)/self.Cwater
         
        B11 = self.Kwa/self.Cwater
        B12 = 1/self.Cwater
        B13 = 1/self.Cwater
       
        A = [[A11]]
        B = [[B11,B12,B13]]
        C = [[1.0]]
        D = [[0.0,0.0,0.0]]

        Ac =  np.array([[A11]])
        Bc =  np.array([[B11,B12,B13]])
        Cc =  np.array([[1.0]])
        Dc =  np.array([[0.0,0.0,0.0]])

        self.sys = signal.StateSpace(A, B, C, D)    
        self.sys_tuple = (Ac, Bc, Cc, Dc)
        
    def sys_dis(self, sys, dt):             
      (A, B, C, D, dt) = signal.cont2discrete(sys, dt, method='zoh', alpha=None)
      self.sys_disc_A = A
      self.sys_disc_B = B
      
    def sys_dis_forOp(self, sys, dt):            
      (A, B, C, D, dt) = signal.cont2discrete(sys, dt, method='zoh', alpha=None)      
      A_opt = A
      B_opt = np.delete(B,0,1)
      E_opt = B[:,0].reshape(A.shape[1],1)  
      self.sys_disc_A = A_opt
      self.sys_disc_B = B_opt   
      self.sys_disc_E = E_opt      
    
    def simulation_sys(self,sys,U,T,Xo): 
        t,y,x = signal.lsim(sys, U, T, Xo, interp=True)
        self.output = y