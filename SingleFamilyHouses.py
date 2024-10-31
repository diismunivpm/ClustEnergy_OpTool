# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 10:14:42 2024

@author: Utente
"""

"""
ClustEnergy OpTool
Released on January, 2024
@author: DIISM UNIVPM
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
from scipy import signal
import numpy as np
import math
from MeteoFile import SF_double_glazing,SF_single_glazing

"""
The classes of Single-Family Houses (SFHs) are modeled below. 
    In particular, buildings with different age classes are considered, each characterized by a different degree of insulation. 
    With reference to the building types contained in the Tabula Project (Corrado et al., 2014), the following classes were considered:
                Class year of construction    	before-1970           	1971-1990	1991-2005    	2006-today
                Single-Family House archetype 	1946-1960             	1976-1990	1991-2005    	2006-today
                Thermal insulation level      	no thermal insulation 	low level	medium level 	high level
                
    To represent the energy demand of individual buildings, lumped-parameter modeling based on the thermoelectric analogy (RC networks) was adopted. 
    From the U values, RC values useful for modeling buildings with networks were calculated (see excel spreadsheet attached in the documentation). 
    In the tool, there are two different model architectures to represent dwelling archetypes. 
    One is for thermally insulated buildings (i.e., SFH 2006-present, SFH 1991-2005, and SFH 1976-1990). For each building envelope element, 
    the front side (layers from the insulation inward) and the back side (layers from the insulation outward) of the insulation are distinguished. 
    In buildings without thermal insulation (SFH 1946-1960), on the other hand, each component of the building envelope is represented by a single thermal node.
    
    Solar gains and internal gains are applied to nodes representing the building envelope. 
    Specifically, following the methodology described in ANSI/ASHRAE Standard 140 (ANSI/ASHRAE, 2020) (i.e., BESTEST) incident solar radiation is applied 
    to the nodes representing the exterior walls and roof. In contrast, the transmitted solar radiation and internal gains are assigned to the nodes 
    representing the internal surfaces. 
    The contribution of air systems is applied directly to the indoor air thermal node. 
    While, the contribution of the underfloor heating system is applied to the node representing the layers from pipes facing inward (Tfi,ap). 
    As a result, the energy dynamics of the building can be represented through the formulation of a state-space model.
    
    REFERENCES:
    Corrado, V., Ballarini, I., Corgnati, S. P., 2014. Building Typology Brochure – Italy. EPISCOPE. Accessed February 2023.
    ANSI/ASHRAE, 2020. ANSI/ASHRAE Standard 140-2017 - Standard Method of Test for the Evaluation of Building Energy Analysis Computer Programs. ASHRAE. Accessed February 2023.
    
    (external calculation of RC parameters allows a generalized form of modeling to be maintained, 
     which is useful for those without knowledge of the building stratigraphy)
    """


class SFH_06(object):
    """
    Single-Family House 2006-today class definition (high level of thermal insulation). 
        Two emission systems are available for defining the archetype:
            - air emission system for space cooling or heating (contribution applied directly on the thermal air node - Tair);
            - floor heating system (contribution applied on the node representing the internal floor surface - Tfi,ap).
    """
 
    def __init__(self, ACH, wall_absorption, radiative_part, convective_part):
        """Definition of building properties. """
        self.ACH         = ACH                  # Air change per hours (1/hr)
        self.wall_absorption = wall_absorption  # Wall absorption
        self.floor_area  = 96.4    # Floor Area (m2)
        self.roof_area   = 96.4    # Roof Area (m2)
        self.volume      = 607.0   # Total volume (m3)
        
        self.wall_wide         = math.sqrt(self.floor_area)            # Building width (m)
        self.wall_long         = math.sqrt(self.floor_area)            # Building length (m)
        self.wall_high         = self.volume/self.floor_area           # Building height (m)
        
        self.windows_area = 21.7 # Total window area (m2)
        self.S_windows_area      = self.windows_area/4        # Total south-facing window area (m2)
        self.E_windows_area      = self.windows_area/4        # Total east-facing window area (m2)
        self.W_windows_area      = self.windows_area/4        # Total west-facing window area (m2)
        self.N_windows_area      = self.windows_area/4        # Total north-facing window area (m2)
        
        self.south_area    = self.wall_wide*self.wall_high - self.S_windows_area # Total south-facing wall area (m2)
        self.north_area    = self.wall_wide*self.wall_high - self.N_windows_area # Total north-facing wall area (m2)
        self.east_area     = self.wall_long*self.wall_high - self.E_windows_area # Total east-facing wall area (m2)
        self.west_area     = self.wall_long*self.wall_high - self.W_windows_area # Total west-facing wall area (m2)
        
        self.total_opaque_area = (self.wall_wide*self.wall_high)*2+ (self.wall_long*self.wall_high)*2 - (self.windows_area) # Total opaque Area (m2)
        self.total_surf_area   = self.floor_area + self.roof_area + self.total_opaque_area                                  # Total surface area (m2)
        
        """Window features of SFH 2006-today (double glazing). """
        self.U_window   = 2.2   # Window thermal transmittance (W/(m2K))
        self.g_window   = 0.724 # Window glass solar factor (%)
        self.emissivity  = 0.2  # Window emissivity
        self.int_window_absorption = 0.156 # Internal window absorption
        self.ext_window_absorption = 0.104 # External window absorption
        self.int_heat_trans_coef = 3.6 + 4.4*(self.emissivity)/0.837 # Internal heat transfer coefficient (W/Km2) - according to UNI EN ISO 10077/1
        self.ext_heat_trans_coef = 23 # External heat transfer coefficient (W/Km2) - according to UNI EN ISO 10077/1
        self.window_reflectance = 0.206 # Window reflectance
        
        """Defining internal distribution fractions of transmitted solar radiation according to the methodology described in ANSI/ASHRAE Standard 140 (ANSI/ASHRAE, 2020). 
            Through these factors, the internal distribution of the solar gains transmitted through the glass surfaces is defined. """

        self.SF_south   = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[0]    # South wall interior solar distribution fraction
        self.SF_east    = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[1]     # East wall interior solar distribution fraction
        self.SF_west    = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[2]     # West wall interior solar distribution fraction
        self.SF_north   = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[3]    # North wall interior solar distribution fraction
        self.SF_ceiling = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[4]  # Ceiling interior solar distribution fraction
        self.SF_floor   = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[5]    # Floor wall interior solar distribution fraction
        
        self.rad_part  = radiative_part  # Radiative part of internal gains
        self.conv_part = convective_part # Convective part of internal gains
    
    """RC properties (calculated externally via excel file - see documentation)."""
    @property
    def R_Window(self):         
        Rwind  = 1/(self.U_window*self.windows_area)         # Window thermal resistance (K/W)
        return (Rwind)
    
    @property
    def NatVent(self):       
        Rinf  = 1/(self.volume*(self.ACH/3600)*1.204*1012)   # Natural venilation thermal resistance (K/W)
        return (Rinf)
    @property
    def Kinfwind(self):
        return (1/self.NatVent)+(1/self.R_Window)            # Total thermal conductance due to ventilation (W/K)
    
    # Walls
    @property
    def Kwi(self):
        return 285    # Thermal conductance of the wall layers from the thermal insulation inward (W/K)
    @property
    def Kwin(self):
        return 106    # Thermal conductance of the thermal insulation layer (W/K)
    @property
    def Kwe(self):
        return 3589   # Thermal conductance of the wall layers from the thermal insulation outward (W/K)
    @property
    def Cwi(self):
        return (62077400.00/3600)    # Thermal capacitance of inner layers (Wh/K)
    @property
    def Cwe(self):
        return (8944616.45/3600)     # Thermal capacitance of the outer layers (Wh/K)
    
    # Roof
    @property
    def Kri(self):
        return 178    # Thermal conductance of the roof layers from the thermal insulation inward (W/K)
    @property
    def Krin(self):
        return 34     # Thermal conductance of the thermal insulation layer (W/K)
    @property
    def Kre(self):
        return 2410   # Thermal conductance of the roof layers from the thermal insulation outward (W/K)
    @property
    def Cri(self):
        return (50513600.00/3600)    # Thermal capacitance of inner layers (Wh/K)
    @property
    def Cre(self):
        return (518651.28/3600)      # Thermal capacitance of the outer layers (Wh/K)
    
    # Floor
    @property
    def Kfi(self):
        return 465      # Thermal conductance of the floor layers from the thermal insulation inward (W/K)
    
    @property
    def Kfiap(self):
        return 596      # Thermal conductance of the floor layers after pipes facing inward (W/K) - underfloor heating system
    @property
    def Kfibp(self):
        return 2129     # Thermal conductance of the floor layers before pipes facing outward (W/K) - underfloor heating system
    
    @property
    def Kfin(self):
        return 35       # Thermal conductance of the thermal insulation layer (W/K)
    @property
    def Kfe(self):
        return 1841     # Thermal conductance of the floor layers from the thermal insulation outward (W/K)
    @property
    def Cfi(self):
        return (23666200.00/3600)    # Thermal capacitance of inner layers (Wh/K)

    @property
    def Cfiap(self):
        return (14411800.00/3600)    # Thermal capacitance of inner layers after pipes facing inward (W/K) - underfloor heating system
    @property
    def Cfibp(self):
        return (9254400.00/3600)     # Thermal capacitance of inner layers before pipes facing outward (W/K) - underfloor heating system
    
    @property
    def Cfe(self):
        return(23650535.00/3600)     # Thermal capacitance of the outer layers (Wh/K)
    
    # Air
    @property
    def Cair(self):
        return self.volume*1.204*1012/3600  # Thermal capacitance of the air (Wh/K)

    def Uvalues(self):
        """Thermal resistance calculation. """
        R_wall   = ((1/self.Kwe)+(1/self.Kwi)+(1/self.Kwin))*self.total_opaque_area
        R_roof   = ((1/self.Kre)+(1/self.Kri)+(1/self.Krin))*self.roof_area
        R_floor  = ((1/self.Kfe)+(1/self.Kfi)+(1/self.Kfin))*self.floor_area
        """Thermal transmittance calculation. """
        self.Uwall  = 1/R_wall
        self.Uroof  = 1/R_roof
        self.Ufloor = 1/R_floor
        """Thermal loss. """
        self.Loss   = ((self.Uwall*self.total_opaque_area)+(self.Uroof*self.roof_area)+(self.Ufloor*self.floor_area)+(self.ACH*self.volume*1.204*1012/3600))

    def stateSpace_air_heating(self) :
        """State space model formulation for air heating system. """
        # Thermal node of the wall layers from the thermal insulation outward (Twe)
        A11 = -(self.Kwe+self.Kwin)/self.Cwe
        A12 = self.Kwin/self.Cwe
        A13 = 0.0
        A14 = 0.0
        A15 = 0.0
        A16 = 0.0
        A17 = 0.0 
               
        B11 = self.Kwe/self.Cwe
        B12 = 0.0
        B13 = 0.0
        B14 = 1/self.Cwe
        B15 = 0.0
        B16 = 0.0
        B17 = 0.0
       
        # Thermal node of the wall layers from the thermal insulation inward (Twi)
        A21 = self.Kwin/self.Cwi
        A22 = -(self.Kwi+self.Kwin)/self.Cwi
        A23 = 0.0
        A24 = 0.0
        A25 = 0.0
        A26 = 0.0
        A27 = self.Kwi/self.Cwi       

        B21 = 0.0
        B22 = 0.0
        B23 = (self.SF_south + self.SF_east + self.SF_west + self.SF_north)/self.Cwi   # solar gains
        B24 = 0.0
        B25 = 0.0
        B26 = (self.rad_part*(self.total_opaque_area/self.total_surf_area))/self.Cwi   # radiative part of internal gains
        B27 = 0.0

        # Thermal node of the roof layers from the thermal insulation outward (Tre)
        A31 = 0.0
        A32 = 0.0
        A33 = -(self.Kre+self.Krin)/self.Cre
        A34 = self.Krin/self.Cre
        A35 = 0.0
        A36 = 0.0
        A37 = 0.0
     
        B31 = self.Kre/self.Cre
        B32 = 0.0
        B33 = 0.0
        B34 = 0.0
        B35 = 1/self.Cre
        B36 = 0.0
        B37 = 0.0
        
        # Thermal node of the roof layers from the thermal insulation inward (Tri)
        A41 = 0.0
        A42 = 0.0
        A43 = self.Krin/self.Cri
        A44 = -(self.Kri+self.Krin)/self.Cri
        A45 = 0.0
        A46 = 0.0
        A47 = self.Kri/self.Cri
       
        B41 = 0.0
        B42 = 0.0
        B43 = self.SF_ceiling/self.Cri # solar gains
        B44 = 0.0
        B45 = 0.0
        B46 = (self.rad_part*(self.roof_area/self.total_surf_area))/self.Cri  # radiative part of internal gains
        B47 = 0.0

        # Thermal node of the floor layers from the thermal insulation outward (Tfe)
        A51 = 0.0
        A52 = 0.0
        A53 = 0.0
        A54 = 0.0
        A55 = -(self.Kfe+self.Kfin)/self.Cfe
        A56 = self.Kfin/self.Cfe
        A57 = 0.0
                
        B51 = 0.0
        B52 = self.Kfe/self.Cfe
        B53 = 0.0
        B54 = 0.0
        B55 = 0.0
        B56 = 0.0 
        B57 = 0.0
        
        # Thermal node of the floor layers from the thermal insulation inward (Tfi)
        A61 = 0.0
        A62 = 0.0
        A63 = 0.0
        A64 = 0.0
        A65 = self.Kfin/self.Cfi
        A66 = -(self.Kfi+self.Kfin)/self.Cfi
        A67 = self.Kfi/self.Cfi        
        
        B61 = 0.0
        B62 = 0.0
        B63 = self.SF_floor/self.Cfi # solar gains
        B64 = 0.0
        B65 = 0.0
        B66 = (self.rad_part*(self.floor_area/self.total_surf_area))/self.Cfi  # radiative part of internal gains
        B67 = 0.0
      
        # Thermal node of the air (Tair)
        A71 = 0.0
        A72 = self.Kwi/self.Cair
        A73 = 0.0
        A74 = self.Kri/self.Cair
        A75 = 0.0
        A76 = self.Kfi/self.Cair
        A77 = -(self.Kinfwind+self.Kwi+self.Kri+self.Kfi)/self.Cair
      
        B71 = (self.Kinfwind)/self.Cair
        B72 = 0.0
        B73 = 0.0
        B74 = 0.0
        B75 = 0.0
        B76 = self.conv_part/self.Cair   # convective part of internal gains
        B77 = 1/self.Cair                # thermal contribution (positive)
                    
        A = [[A11,A12,A13,A14,A15,A16,A17],[A21,A22,A23,A24,A25,A26,A27],[A31,A32,A33,A34,A35,A36,A37],[A41,A42,A43,A44,A45,A46,A47],[A51,A52,A53,A54,A55,A56,A57],[A61,A62,A63,A64,A65,A66,A67],[A71,A72,A73,A74,A75,A76,A77]]
        B = [[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47],[B51,B52,B53,B54,B55,B56,B57],[B61,B62,B63,B64,B65,B66,B67],[B71,B72,B73,B74,B75,B76,B77]]
        C = [[1.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,1.0]]
        D = [[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]]

        Ac =  np.array([[A11,A12,A13,A14,A15,A16,A17],[A21,A22,A23,A24,A25,A26,A27],[A31,A32,A33,A34,A35,A36,A37],[A41,A42,A43,A44,A45,A46,A47],[A51,A52,A53,A54,A55,A56,A57],[A61,A62,A63,A64,A65,A66,A67],[A71,A72,A73,A74,A75,A76,A77]])
        Bc =  np.array([[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47],[B51,B52,B53,B54,B55,B56,B57],[B61,B62,B63,B64,B65,B66,B67],[B71,B72,B73,B74,B75,B76,B77]])
        Cc =  np.array([[1.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,1.0]])
        Dc =  np.array([[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]])

        self.sys = signal.StateSpace(A, B, C, D)    
        self.sys_tuple = (Ac, Bc, Cc, Dc)

    def stateSpace_air_cooling(self) :
        """State space model formulation for air cooling system. """
        # Thermal node of the wall layers from the thermal insulation outward (Twe)
        A11 = -(self.Kwe+self.Kwin)/self.Cwe
        A12 = self.Kwin/self.Cwe
        A13 = 0.0
        A14 = 0.0
        A15 = 0.0
        A16 = 0.0
        A17 = 0.0 
               
        B11 = self.Kwe/self.Cwe
        B12 = 0.0
        B13 = 0.0
        B14 = 1/self.Cwe
        B15 = 0.0
        B16 = 0.0
        B17 = 0.0
       
        # Thermal node of the wall layers from the thermal insulation inward (Twi)
        A21 = self.Kwin/self.Cwi
        A22 = -(self.Kwi+self.Kwin)/self.Cwi
        A23 = 0.0
        A24 = 0.0
        A25 = 0.0
        A26 = 0.0
        A27 = self.Kwi/self.Cwi       
    
        B21 = 0.0
        B22 = 0.0
        B23 = (self.SF_south + self.SF_east + self.SF_west + self.SF_north)/self.Cwi   # solar gains
        B24 = 0.0
        B25 = 0.0
        B26 = (self.rad_part*(self.total_opaque_area/self.total_surf_area))/self.Cwi   # radiative part of internal gains
        B27 = 0.0
    
        # Thermal node of the roof layers from the thermal insulation outward (Tre)
        A31 = 0.0
        A32 = 0.0
        A33 = -(self.Kre+self.Krin)/self.Cre
        A34 = self.Krin/self.Cre
        A35 = 0.0
        A36 = 0.0
        A37 = 0.0
     
        B31 = self.Kre/self.Cre
        B32 = 0.0
        B33 = 0.0
        B34 = 0.0
        B35 = 1/self.Cre
        B36 = 0.0
        B37 = 0.0
        
        # Thermal node of the roof layers from the thermal insulation inward (Tri)
        A41 = 0.0
        A42 = 0.0
        A43 = self.Krin/self.Cri
        A44 = -(self.Kri+self.Krin)/self.Cri
        A45 = 0.0
        A46 = 0.0
        A47 = self.Kri/self.Cri
       
        B41 = 0.0
        B42 = 0.0
        B43 = self.SF_ceiling/self.Cri # solar gains
        B44 = 0.0
        B45 = 0.0
        B46 = (self.rad_part*(self.roof_area/self.total_surf_area))/self.Cri  # radiative part of internal gains
        B47 = 0.0
    
        # Thermal node of the floor layers from the thermal insulation outward (Tfe)
        A51 = 0.0
        A52 = 0.0
        A53 = 0.0
        A54 = 0.0
        A55 = -(self.Kfe+self.Kfin)/self.Cfe
        A56 = self.Kfin/self.Cfe
        A57 = 0.0
                
        B51 = 0.0
        B52 = self.Kfe/self.Cfe
        B53 = 0.0
        B54 = 0.0
        B55 = 0.0
        B56 = 0.0 
        B57 = 0.0
        
        # Thermal node of the floor layers from the thermal insulation inward (Tfi)
        A61 = 0.0
        A62 = 0.0
        A63 = 0.0
        A64 = 0.0
        A65 = self.Kfin/self.Cfi
        A66 = -(self.Kfi+self.Kfin)/self.Cfi
        A67 = self.Kfi/self.Cfi        
        
        B61 = 0.0
        B62 = 0.0
        B63 = self.SF_floor/self.Cfi # solar gains
        B64 = 0.0
        B65 = 0.0
        B66 = (self.rad_part*(self.floor_area/self.total_surf_area))/self.Cfi  # radiative part of internal gains  
        B67 = 0.0
      
        # Thermal node of the air (Tair)
        A71 = 0.0
        A72 = self.Kwi/self.Cair
        A73 = 0.0
        A74 = self.Kri/self.Cair
        A75 = 0.0
        A76 = self.Kfi/self.Cair
        A77 = -(self.Kinfwind+self.Kwi+self.Kri+self.Kfi)/self.Cair
      
        B71 = (self.Kinfwind)/self.Cair
        B72 = 0.0
        B73 = 0.0
        B74 = 0.0
        B75 = 0.0
        B76 = self.conv_part/self.Cair   # convective part of internal gains
        B77 = (- 1)/self.Cair            # thermal contribution (negative)
                    
        A = [[A11,A12,A13,A14,A15,A16,A17],[A21,A22,A23,A24,A25,A26,A27],[A31,A32,A33,A34,A35,A36,A37],[A41,A42,A43,A44,A45,A46,A47],[A51,A52,A53,A54,A55,A56,A57],[A61,A62,A63,A64,A65,A66,A67],[A71,A72,A73,A74,A75,A76,A77]]
        B = [[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47],[B51,B52,B53,B54,B55,B56,B57],[B61,B62,B63,B64,B65,B66,B67],[B71,B72,B73,B74,B75,B76,B77]]
        C = [[1.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,1.0]]
        D = [[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]]
    
        Ac =  np.array([[A11,A12,A13,A14,A15,A16,A17],[A21,A22,A23,A24,A25,A26,A27],[A31,A32,A33,A34,A35,A36,A37],[A41,A42,A43,A44,A45,A46,A47],[A51,A52,A53,A54,A55,A56,A57],[A61,A62,A63,A64,A65,A66,A67],[A71,A72,A73,A74,A75,A76,A77]])
        Bc =  np.array([[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47],[B51,B52,B53,B54,B55,B56,B57],[B61,B62,B63,B64,B65,B66,B67],[B71,B72,B73,B74,B75,B76,B77]])
        Cc =  np.array([[1.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,1.0]])
        Dc =  np.array([[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]])
    
        self.sys = signal.StateSpace(A, B, C, D)    
        self.sys_tuple = (Ac, Bc, Cc, Dc)

    def stateSpace_floor_heating(self) :
        """State space model formulation for underfloor heating system. """
        # Thermal node of the wall layers from the thermal insulation outward (Twe)
        A11 = -(self.Kwe+self.Kwin)/self.Cwe
        A12 = self.Kwin/self.Cwe
        A13 = 0.0
        A14 = 0.0
        A15 = 0.0
        A16 = 0.0
        A17 = 0.0 
        A18 = 0.0 
        
        B11 = self.Kwe/self.Cwe
        B12 = 0.0
        B13 = 0.0
        B14 = 1/self.Cwe
        B15 = 0.0
        B16 = 0.0
        B17 = 0.0
        
        # Thermal node of the wall layers from the thermal insulation inward (Twi)
        A21 = self.Kwin/self.Cwi
        A22 = -(self.Kwi+self.Kwin)/self.Cwi
        A23 = 0.0
        A24 = 0.0
        A25 = 0.0
        A26 = 0.0
        A27 = 0  
        A28 = self.Kwi/self.Cwi  

        B21 = 0.0
        B22 = 0.0
        B23 = (self.SF_south + self.SF_east + self.SF_west + self.SF_north)/self.Cwi   # solar gains
        B24 = 0.0
        B25 = 0.0
        B26 = (self.rad_part*(self.total_opaque_area/self.total_surf_area))/self.Cwi   # radiative part of internal gains
        B27 = 0.0

        # Thermal node of the roof layers from the thermal insulation outward (Tre)
        A31 = 0.0
        A32 = 0.0
        A33 = -(self.Kre+self.Krin)/self.Cre
        A34 = self.Krin/self.Cre
        A35 = 0.0
        A36 = 0.0
        A37 = 0.0
        A38 = 0.0
        
        B31 = self.Kre/self.Cre
        B32 = 0.0
        B33 = 0.0
        B34 = 0.0
        B35 = 1/self.Cre
        B36 = 0.0
        B37 = 0.0
        
        # Thermal node of the roof layers from the thermal insulation inward (Tri)
        A41 = 0.0
        A42 = 0.0
        A43 = self.Krin/self.Cri
        A44 = -(self.Kri+self.Krin)/self.Cri
        A45 = 0.0
        A46 = 0.0
        A47 = 0.0
        A48 = self.Kri/self.Cri
        
        B41 = 0.0
        B42 = 0.0
        B43 = self.SF_ceiling/self.Cri # solar gains
        B44 = 0.0
        B45 = 0.0
        B46 = (self.rad_part*(self.roof_area/self.total_surf_area))/self.Cri  # radiative part of internal gains
        B47 = 0.0

        # Thermal node of the floor layers from the thermal insulation outward (Tfe)
        A51 = 0.0
        A52 = 0.0
        A53 = 0.0
        A54 = 0.0
        A55 = -(self.Kfe+self.Kfin)/self.Cfe
        A56 = self.Kfin/self.Cfe
        A57 = 0.0
        A58 = 0.0 
        
        B51 = 0.0
        B52 = self.Kfe/self.Cfe
        B53 = 0.0
        B54 = 0.0
        B55 = 0.0
        B56 = 0.0 
        B57 = 0.0
        
        # Thermal node of inner layers before pipes facing outward (Tfi,bp)
        A61 = 0.0
        A62 = 0.0
        A63 = 0.0
        A64 = 0.0
        A65 = self.Kfin/self.Cfibp
        A66 = -(self.Kfibp+self.Kfin)/self.Cfibp
        A67 = self.Kfibp/self.Cfibp        
        A68 = 0.0
        
        B61 = 0.0
        B62 = 0.0
        B63 = 0.0
        B64 = 0.0
        B65 = 0.0
        B66 = 0.0 
        B67 = 0.0

        # Thermal node of inner layers after pipes facing inward (Tfi,ap)
        A71 = 0.0
        A72 = 0.0
        A73 = 0.0
        A74 = 0.0
        A75 = 0.0
        A76 = self.Kfibp/self.Cfiap
        A77 = -(self.Kfibp+self.Kfiap)/self.Cfiap
        A78 = self.Kfiap/self.Cfiap      
        
        B71 = 0.0
        B72 = 0.0
        B73 = self.SF_floor/self.Cfi # solar gains
        B74 = 0.0
        B75 = 0.0
        B76 = (self.rad_part*(self.floor_area/self.total_surf_area))/self.Cfi  # radiative part of internal gains  
        B77 = 1/self.Cfiap # thermal contribution (positive)
        
        # Thermal node of the air (Tair)
        A81 = 0.0
        A82 = self.Kwi/self.Cair
        A83 = 0.0
        A84 = self.Kri/self.Cair
        A85 = 0.0
        A86 = 0.0
        A87 = self.Kfiap/self.Cair
        A88 = -(self.Kinfwind+self.Kwi+self.Kri+self.Kfiap)/self.Cair
        
        B81 = (self.Kinfwind)/self.Cair
        B82 = 0.0
        B83 = 0.0
        B84 = 0.0
        B85 = 0.0
        B86 = self.conv_part/self.Cair # convective part of internal gains
        B87 = 0.0
        
        A = [[A11,A12,A13,A14,A15,A16,A17,A18],[A21,A22,A23,A24,A25,A26,A27,A28],[A31,A32,A33,A34,A35,A36,A37,A38],[A41,A42,A43,A44,A45,A46,A47,A48],[A51,A52,A53,A54,A55,A56,A57,A58],[A61,A62,A63,A64,A65,A66,A67,A68],[A71,A72,A73,A74,A75,A76,A77,A78],[A81,A82,A83,A84,A85,A86,A87,A88]]
        B = [[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47],[B51,B52,B53,B54,B55,B56,B57],[B61,B62,B63,B64,B65,B66,B67],[B71,B72,B73,B74,B75,B76,B77],[B81,B82,B83,B84,B85,B86,B87]]
        C = [[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0]]
        D = [[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]]

        Ac =  np.array([[A11,A12,A13,A14,A15,A16,A17,A18],[A21,A22,A23,A24,A25,A26,A27,A28],[A31,A32,A33,A34,A35,A36,A37,A38],[A41,A42,A43,A44,A45,A46,A47,A48],[A51,A52,A53,A54,A55,A56,A57,A58],[A61,A62,A63,A64,A65,A66,A67,A68],[A71,A72,A73,A74,A75,A76,A77,A78],[A81,A82,A83,A84,A85,A86,A87,A88]])
        Bc =  np.array([[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47],[B51,B52,B53,B54,B55,B56,B57],[B61,B62,B63,B64,B65,B66,B67],[B71,B72,B73,B74,B75,B76,B77],[B81,B82,B83,B84,B85,B86,B87]])
        Cc =  np.array([[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0]])
        Dc =  np.array([[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]])

        self.sys = signal.StateSpace(A, B, C, D)    
        self.sys_tuple = (Ac, Bc, Cc, Dc)
        
    def sys_dis(self, sys, dt):             
      (A, B, C, D, dt) = signal.cont2discrete(sys, dt, method='zoh', alpha=None)
      self.sys_disc_A = A
      self.sys_disc_B = B
      
    def sys_dis_forOp(self, sys, dt):            
      (A, B, C, D, dt) = signal.cont2discrete(sys, dt, method='zoh', alpha=None)      
      A_opt = A
      B_opt = np.delete(B,6,1)
      E_opt = B[:,6].reshape(A.shape[1],1)  
      self.sys_disc_A = A_opt
      self.sys_disc_B = B_opt   
      self.sys_disc_E = E_opt      
    
    def simulation_sys(self,sys,U,T,Xo): 
        t,y,x = signal.lsim(sys, U, T, Xo, interp=True)
        self.output = y

class SFH_05(object):
    """
    Single-Family House 1991-2005 class definition (medium level of thermal insulation). 
        Only air emission system is available for the archetype (contribution applied directly on the thermal air node - Tair).
    """
    def __init__(self, ACH, wall_absorption, radiative_part, convective_part):
        """Definition of building properties. """
        self.ACH         = ACH                  # Air change per hours (1/hr)
        self.wall_absorption = wall_absorption  # Wall absorption
        self.floor_area  = 96.0    # Floor Area (m2)
        self.roof_area   = 96.0    # Roof Area (m2)
        self.volume      = 605     # Total volume (m3)
        
        self.wall_wide         = math.sqrt(self.floor_area)            # Building width (m)
        self.wall_long         = math.sqrt(self.floor_area)            # Building length (m)
        self.wall_high         = self.volume/self.floor_area           # Building height (m)
        
        self.windows_area = 21.5 # Total window area (m2)
        self.S_windows_area      = self.windows_area/4        # Total south-facing window area (m2)
        self.E_windows_area      = self.windows_area/4        # Total east-facing window area (m2)
        self.W_windows_area      = self.windows_area/4        # Total west-facing window area (m2)
        self.N_windows_area      = self.windows_area/4        # Total north-facing window area (m2)
        
        self.south_area    = self.wall_wide*self.wall_high - self.S_windows_area # Total south-facing wall area (m2)
        self.north_area    = self.wall_wide*self.wall_high - self.N_windows_area # Total north-facing wall area (m2)
        self.east_area     = self.wall_long*self.wall_high - self.E_windows_area # Total east-facing wall area (m2)
        self.west_area     = self.wall_long*self.wall_high - self.W_windows_area # Total west-facing wall area (m2)
        
        self.total_opaque_area = (self.wall_wide*self.wall_high)*2+ (self.wall_long*self.wall_high)*2 - (self.windows_area) # Total opaque Area (m2)
        self.total_surf_area   = self.floor_area + self.roof_area + self.total_opaque_area                                  # Total surface area (m2)
        
        """Window features of SFH 1991-2005 (double glazing). """
        self.U_window   = 2.8   # Window thermal transmittance (W/(m2K))
        self.g_window   = 0.722 # Window glass solar factor (%)
        self.emissivity  = 0.9  # Window emissivity
        self.int_window_absorption = 0.156 # Internal window absorption
        self.ext_window_absorption = 0.104 # External window absorption
        self.int_heat_trans_coef = 3.6 + 4.4*(self.emissivity)/0.837 # Internal heat transfer coefficient (W/Km2) - according to UNI EN ISO 10077/1
        self.ext_heat_trans_coef = 23 # External heat transfer coefficient (W/Km2) - according to UNI EN ISO 10077/1
        self.window_reflectance = 0.206 # Window reflectance
        
        """Defining internal distribution fractions of transmitted solar radiation according to the methodology described in ANSI/ASHRAE Standard 140 (ANSI/ASHRAE, 2020). 
            Through these factors, the internal distribution of the solar gains transmitted through the glass surfaces is defined. """

        self.SF_south   = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[0]    # South wall interior solar distribution fraction
        self.SF_east    = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[1]     # East wall interior solar distribution fraction
        self.SF_west    = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[2]     # West wall interior solar distribution fraction
        self.SF_north   = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[3]    # North wall interior solar distribution fraction
        self.SF_ceiling = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[4]  # Ceiling interior solar distribution fraction
        self.SF_floor   = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[5]    # Floor wall interior solar distribution fraction
        
        self.rad_part  = radiative_part  # Radiative part of internal gains
        self.conv_part = convective_part # Convective part of internal gains
        
    """RC properties (calculated externally via excel file - see documentation)."""
    @property
    def R_Window(self):         
        Rwind  = 1/(self.U_window*self.windows_area)         # Window thermal resistance (K/W)
        return (Rwind)
    
    @property
    def NatVent(self):       
        Rinf  = 1/(self.volume*(self.ACH/3600)*1.204*1012)   # Natural venilation thermal resistance (K/W)
        return (Rinf)
    @property
    def Kinfwind(self):
        return (1/self.NatVent)+(1/self.R_Window)            # Total thermal conductance due to ventilation (W/K)
    
    # Walls
    @property
    def Kwi(self):
        return 285    # Thermal conductance of the wall layers from the thermal insulation inward (W/K)
    @property
    def Kwin(self):
        return 264    # Thermal conductance of the thermal insulation layer (W/K)
    @property
    def Kwe(self):
        return 3586   # Thermal conductance of the wall layers from the thermal insulation outward (W/K)
    @property
    def Cwi(self):
        return (62021800.00/3600)    # Thermal capacitance of inner layers (Wh/K)
    @property
    def Cwe(self):
        return  (8393602.06/3600)    # Thermal capacitance of the outer layers (Wh/K)
   
    # Roof
    @property
    def Kri(self):
        return 177    # Thermal conductance of the roof layers from the thermal insulation inward (W/K)
    @property
    def Krin(self):
        return 111     # Thermal conductance of the thermal insulation layer (W/K)
    @property
    def Kre(self):
        return 2400   # Thermal conductance of the roof layers from the thermal insulation outward (W/K)
    @property
    def Cri(self):
        return (50304000.00/3600)    # Thermal capacitance of inner layers (Wh/K)
    @property
    def Cre(self):
        return (159868.80/3600)      # Thermal capacitance of the outer layers (Wh/K)
    
    # Floor
    @property
    def Kfi(self):
        return 464      # Thermal conductance of the floor layers from the thermal insulation inward (W/K)
    @property
    def Kfin(self):
        return 92       # Thermal conductance of the thermal insulation layer (W/K)
    @property
    def Kfe(self):
        return 1834     # Thermal conductance of the floor layers from the thermal insulation outward (W/K)
    @property
    def Cfi(self):
        return (23568000.00/3600)    # Thermal capacitance of inner layers (Wh/K)
    @property
    def Cfe(self):
        return(23232662.40/3.6)     # Thermal capacitance of the outer layers (Wh/K)
    @property
    
    # Air
    def Cair(self):
        return self.volume*1.204*1012/3600  # Thermal capacitance of the air (Wh/K)

    def Uvalues(self):
        """Thermal resistance calculation. """
        R_wall   = ((1/self.Kwe)+(1/self.Kwi)+(1/self.Kwin))*self.total_opaque_area
        R_roof   = ((1/self.Kre)+(1/self.Kri)+(1/self.Krin))*self.roof_area
        R_floor  = ((1/self.Kfe)+(1/self.Kfi)+(1/self.Kfin))*self.floor_area
        """Thermal transmittance calculation. """
        self.Uwall  = 1/R_wall
        self.Uroof  = 1/R_roof
        self.Ufloor = 1/R_floor
        """Thermal loss. """
        self.Loss   = ((self.Uwall*self.total_opaque_area)+(self.Uroof*self.roof_area)+(self.Ufloor*self.floor_area)+(self.ACH*self.volume*1.204*1012/3600))

    def stateSpace_air_heating(self) :
        """State space model formulation for air heating system. """
        # Thermal node of the wall layers from the thermal insulation outward (Twe)
        A11 = -(self.Kwe+self.Kwin)/self.Cwe
        A12 = self.Kwin/self.Cwe
        A13 = 0.0
        A14 = 0.0
        A15 = 0.0
        A16 = 0.0
        A17 = 0.0 
               
        B11 = self.Kwe/self.Cwe
        B12 = 0.0
        B13 = 0.0
        B14 = 1/self.Cwe
        B15 = 0.0
        B16 = 0.0
        B17 = 0.0
       
        # Thermal node of the wall layers from the thermal insulation inward (Twi)
        A21 = self.Kwin/self.Cwi
        A22 = -(self.Kwi+self.Kwin)/self.Cwi
        A23 = 0.0
        A24 = 0.0
        A25 = 0.0
        A26 = 0.0
        A27 = self.Kwi/self.Cwi       

        B21 = 0.0
        B22 = 0.0
        B23 = (self.SF_south + self.SF_east + self.SF_west + self.SF_north)/self.Cwi   # solar gains
        B24 = 0.0
        B25 = 0.0
        B26 = (self.rad_part*(self.total_opaque_area/self.total_surf_area))/self.Cwi   # radiative part of internal gains
        B27 = 0.0

        # Thermal node of the roof layers from the thermal insulation outward (Tre)
        A31 = 0.0
        A32 = 0.0
        A33 = -(self.Kre+self.Krin)/self.Cre
        A34 = self.Krin/self.Cre
        A35 = 0.0
        A36 = 0.0
        A37 = 0.0
     
        B31 = self.Kre/self.Cre
        B32 = 0.0
        B33 = 0.0
        B34 = 0.0
        B35 = 1/self.Cre
        B36 = 0.0
        B37 = 0.0
        
        # Thermal node of the roof layers from the thermal insulation inward (Tri)
        A41 = 0.0
        A42 = 0.0
        A43 = self.Krin/self.Cri
        A44 = -(self.Kri+self.Krin)/self.Cri
        A45 = 0.0
        A46 = 0.0
        A47 = self.Kri/self.Cri
       
        B41 = 0.0
        B42 = 0.0
        B43 = self.SF_ceiling/self.Cri # solar gains
        B44 = 0.0
        B45 = 0.0
        B46 = (self.rad_part*(self.roof_area/self.total_surf_area))/self.Cri  # radiative part of internal gains
        B47 = 0.0

        # Thermal node of the floor layers from the thermal insulation outward (Tfe)
        A51 = 0.0
        A52 = 0.0
        A53 = 0.0
        A54 = 0.0
        A55 = -(self.Kfe+self.Kfin)/self.Cfe
        A56 = self.Kfin/self.Cfe
        A57 = 0.0
                
        B51 = 0.0
        B52 = self.Kfe/self.Cfe
        B53 = 0.0
        B54 = 0.0
        B55 = 0.0
        B56 = 0.0 
        B57 = 0.0
        
        # Thermal node of the floor layers from the thermal insulation inward (Tfi)
        A61 = 0.0
        A62 = 0.0
        A63 = 0.0
        A64 = 0.0
        A65 = self.Kfin/self.Cfi
        A66 = -(self.Kfi+self.Kfin)/self.Cfi
        A67 = self.Kfi/self.Cfi        
        
        B61 = 0.0
        B62 = 0.0
        B63 = self.SF_floor/self.Cfi # solar gains
        B64 = 0.0
        B65 = 0.0
        B66 = (self.rad_part*(self.floor_area/self.total_surf_area))/self.Cfi  # radiative part of internal gains
        B67 = 0.0
      
        # Thermal node of the air (Tair)
        A71 = 0.0
        A72 = self.Kwi/self.Cair
        A73 = 0.0
        A74 = self.Kri/self.Cair
        A75 = 0.0
        A76 = self.Kfi/self.Cair
        A77 = -(self.Kinfwind+self.Kwi+self.Kri+self.Kfi)/self.Cair
      
        B71 = (self.Kinfwind)/self.Cair
        B72 = 0.0
        B73 = 0.0
        B74 = 0.0
        B75 = 0.0
        B76 = self.conv_part/self.Cair   # convective part of internal gains
        B77 = 1/self.Cair                # thermal contribution (positive)
                    
        A = [[A11,A12,A13,A14,A15,A16,A17],[A21,A22,A23,A24,A25,A26,A27],[A31,A32,A33,A34,A35,A36,A37],[A41,A42,A43,A44,A45,A46,A47],[A51,A52,A53,A54,A55,A56,A57],[A61,A62,A63,A64,A65,A66,A67],[A71,A72,A73,A74,A75,A76,A77]]
        B = [[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47],[B51,B52,B53,B54,B55,B56,B57],[B61,B62,B63,B64,B65,B66,B67],[B71,B72,B73,B74,B75,B76,B77]]
        C = [[1.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,1.0]]
        D = [[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]]

        Ac =  np.array([[A11,A12,A13,A14,A15,A16,A17],[A21,A22,A23,A24,A25,A26,A27],[A31,A32,A33,A34,A35,A36,A37],[A41,A42,A43,A44,A45,A46,A47],[A51,A52,A53,A54,A55,A56,A57],[A61,A62,A63,A64,A65,A66,A67],[A71,A72,A73,A74,A75,A76,A77]])
        Bc =  np.array([[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47],[B51,B52,B53,B54,B55,B56,B57],[B61,B62,B63,B64,B65,B66,B67],[B71,B72,B73,B74,B75,B76,B77]])
        Cc =  np.array([[1.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,1.0]])
        Dc =  np.array([[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]])

        self.sys = signal.StateSpace(A, B, C, D)    
        self.sys_tuple = (Ac, Bc, Cc, Dc)

    def stateSpace_air_cooling(self) :
        """State space model formulation for air cooling system. """
        # Thermal node of the wall layers from the thermal insulation outward (Twe)
        A11 = -(self.Kwe+self.Kwin)/self.Cwe
        A12 = self.Kwin/self.Cwe
        A13 = 0.0
        A14 = 0.0
        A15 = 0.0
        A16 = 0.0
        A17 = 0.0 
               
        B11 = self.Kwe/self.Cwe
        B12 = 0.0
        B13 = 0.0
        B14 = 1/self.Cwe
        B15 = 0.0
        B16 = 0.0
        B17 = 0.0
       
        # Thermal node of the wall layers from the thermal insulation inward (Twi)
        A21 = self.Kwin/self.Cwi
        A22 = -(self.Kwi+self.Kwin)/self.Cwi
        A23 = 0.0
        A24 = 0.0
        A25 = 0.0
        A26 = 0.0
        A27 = self.Kwi/self.Cwi       
    
        B21 = 0.0
        B22 = 0.0
        B23 = (self.SF_south + self.SF_east + self.SF_west + self.SF_north)/self.Cwi   # solar gains
        B24 = 0.0
        B25 = 0.0
        B26 = (self.rad_part*(self.total_opaque_area/self.total_surf_area))/self.Cwi   # radiative part of internal gains
        B27 = 0.0
    
        # Thermal node of the roof layers from the thermal insulation outward (Tre)
        A31 = 0.0
        A32 = 0.0
        A33 = -(self.Kre+self.Krin)/self.Cre
        A34 = self.Krin/self.Cre
        A35 = 0.0
        A36 = 0.0
        A37 = 0.0
     
        B31 = self.Kre/self.Cre
        B32 = 0.0
        B33 = 0.0
        B34 = 0.0
        B35 = 1/self.Cre
        B36 = 0.0
        B37 = 0.0
        
        # Thermal node of the roof layers from the thermal insulation inward (Tri)
        A41 = 0.0
        A42 = 0.0
        A43 = self.Krin/self.Cri
        A44 = -(self.Kri+self.Krin)/self.Cri
        A45 = 0.0
        A46 = 0.0
        A47 = self.Kri/self.Cri
       
        B41 = 0.0
        B42 = 0.0
        B43 = self.SF_ceiling/self.Cri # solar gains
        B44 = 0.0
        B45 = 0.0
        B46 = (self.rad_part*(self.roof_area/self.total_surf_area))/self.Cri  # radiative part of internal gains
        B47 = 0.0
    
        # Thermal node of the floor layers from the thermal insulation outward (Tfe)
        A51 = 0.0
        A52 = 0.0
        A53 = 0.0
        A54 = 0.0
        A55 = -(self.Kfe+self.Kfin)/self.Cfe
        A56 = self.Kfin/self.Cfe
        A57 = 0.0
                
        B51 = 0.0
        B52 = self.Kfe/self.Cfe
        B53 = 0.0
        B54 = 0.0
        B55 = 0.0
        B56 = 0.0 
        B57 = 0.0
        
        # Thermal node of the floor layers from the thermal insulation inward (Tfi)
        A61 = 0.0
        A62 = 0.0
        A63 = 0.0
        A64 = 0.0
        A65 = self.Kfin/self.Cfi
        A66 = -(self.Kfi+self.Kfin)/self.Cfi
        A67 = self.Kfi/self.Cfi        
        
        B61 = 0.0
        B62 = 0.0
        B63 = self.SF_floor/self.Cfi # solar gains
        B64 = 0.0
        B65 = 0.0
        B66 = (self.rad_part*(self.floor_area/self.total_surf_area))/self.Cfi  # radiative part of internal gains  
        B67 = 0.0
      
        # Thermal node of the air (Tair)
        A71 = 0.0
        A72 = self.Kwi/self.Cair
        A73 = 0.0
        A74 = self.Kri/self.Cair
        A75 = 0.0
        A76 = self.Kfi/self.Cair
        A77 = -(self.Kinfwind+self.Kwi+self.Kri+self.Kfi)/self.Cair
      
        B71 = (self.Kinfwind)/self.Cair
        B72 = 0.0
        B73 = 0.0
        B74 = 0.0
        B75 = 0.0
        B76 = self.conv_part/self.Cair   # convective part of internal gains
        B77 = (- 1)/self.Cair            # thermal contribution (negative)
                    
        A = [[A11,A12,A13,A14,A15,A16,A17],[A21,A22,A23,A24,A25,A26,A27],[A31,A32,A33,A34,A35,A36,A37],[A41,A42,A43,A44,A45,A46,A47],[A51,A52,A53,A54,A55,A56,A57],[A61,A62,A63,A64,A65,A66,A67],[A71,A72,A73,A74,A75,A76,A77]]
        B = [[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47],[B51,B52,B53,B54,B55,B56,B57],[B61,B62,B63,B64,B65,B66,B67],[B71,B72,B73,B74,B75,B76,B77]]
        C = [[1.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,1.0]]
        D = [[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]]
    
        Ac =  np.array([[A11,A12,A13,A14,A15,A16,A17],[A21,A22,A23,A24,A25,A26,A27],[A31,A32,A33,A34,A35,A36,A37],[A41,A42,A43,A44,A45,A46,A47],[A51,A52,A53,A54,A55,A56,A57],[A61,A62,A63,A64,A65,A66,A67],[A71,A72,A73,A74,A75,A76,A77]])
        Bc =  np.array([[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47],[B51,B52,B53,B54,B55,B56,B57],[B61,B62,B63,B64,B65,B66,B67],[B71,B72,B73,B74,B75,B76,B77]])
        Cc =  np.array([[1.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,1.0]])
        Dc =  np.array([[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]])
    
        self.sys = signal.StateSpace(A, B, C, D)    
        self.sys_tuple = (Ac, Bc, Cc, Dc)
        
    def sys_dis(self, sys, dt):             
      (A, B, C, D, dt) = signal.cont2discrete(sys, dt, method='zoh', alpha=None)
      self.sys_disc_A = A
      self.sys_disc_B = B
      
    def sys_dis_forOp(self, sys, dt):            
      (A, B, C, D, dt) = signal.cont2discrete(sys, dt, method='zoh', alpha=None)      
      A_opt = A
      B_opt = np.delete(B,6,1)
      E_opt = B[:,6].reshape(A.shape[1],1)  
      self.sys_disc_A = A_opt
      self.sys_disc_B = B_opt   
      self.sys_disc_E = E_opt      
    
    def simulation_sys(self,sys,U,T,Xo): 
        t,y,x = signal.lsim(sys, U, T, Xo, interp=True)
        self.output = y

class SFH_90(object):
    """
    Single-Family House 1976-1990 class definition (low level of thermal insulation). 
        Only air emission system is available for the archetype (contribution applied directly on the thermal air node - Tair).
    """
    def __init__(self,ACH, wall_absorption, radiative_part, convective_part):
        """Definition of building properties. """
        self.ACH         = ACH                  # Air change per hours (1/hr)
        self.wall_absorption = wall_absorption  # Wall absorption
        self.floor_area  = 115.1    # Floor Area (m2)
        self.roof_area   = 132.9    # Roof Area (m2)
        self.volume      = 725.0    # Total volume (m3)
        
        self.wall_wide         = math.sqrt(self.floor_area)            # Building width (m)
        self.wall_long         = math.sqrt(self.floor_area)            # Building length (m)
        self.wall_high         = self.volume/self.floor_area           # Building height (m)
        
        self.windows_area = 20.9 # Total window area (m2)
        self.S_windows_area      = self.windows_area/4        # Total south-facing window area (m2)
        self.E_windows_area      = self.windows_area/4        # Total east-facing window area (m2)
        self.W_windows_area      = self.windows_area/4        # Total west-facing window area (m2)
        self.N_windows_area      = self.windows_area/4        # Total north-facing window area (m2)
        
        self.south_area    = self.wall_wide*self.wall_high - self.S_windows_area # Total south-facing wall area (m2)
        self.north_area    = self.wall_wide*self.wall_high - self.N_windows_area # Total north-facing wall area (m2)
        self.east_area     = self.wall_long*self.wall_high - self.E_windows_area # Total east-facing wall area (m2)
        self.west_area     = self.wall_long*self.wall_high - self.W_windows_area # Total west-facing wall area (m2)
        
        self.total_opaque_area = (self.wall_wide*self.wall_high)*2+ (self.wall_long*self.wall_high)*2 - (self.windows_area) # Total opaque Area (m2)
        self.total_surf_area   = self.floor_area + self.roof_area + self.total_opaque_area                                  # Total surface area (m2)
        
        """Window features of SFH 1976-1990 (double glazing). """
        self.U_window   = 2.8    # Window thermal transmittance (W/(m2K))
        self.g_window   = 0.722  # Window glass solar factor (%)
        self.emissivity  = 0.84  # Window emissivity
        self.int_window_absorption = 0.156 # Internal window absorption
        self.ext_window_absorption = 0.104 # External window absorption
        self.int_heat_trans_coef = 3.6 + 4.4*(self.emissivity)/0.837 # Internal heat transfer coefficient (W/Km2) - according to UNI EN ISO 10077/1
        self.ext_heat_trans_coef = 23 # External heat transfer coefficient (W/Km2) - according to UNI EN ISO 10077/1
        self.window_reflectance = 0.206 # Window reflectance
        
        """Defining internal distribution fractions of transmitted solar radiation according to the methodology described in ANSI/ASHRAE Standard 140 (ANSI/ASHRAE, 2020). 
            Through these factors, the internal distribution of the solar gains transmitted through the glass surfaces is defined. """

        self.SF_south   = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[0]    # South wall interior solar distribution fraction
        self.SF_east    = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[1]     # East wall interior solar distribution fraction
        self.SF_west    = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[2]     # West wall interior solar distribution fraction
        self.SF_north   = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[3]    # North wall interior solar distribution fraction
        self.SF_ceiling = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[4]  # Ceiling interior solar distribution fraction
        self.SF_floor   = SF_double_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high,  
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, self.N_windows_area, 
                           self.wall_absorption, self.U_window, self.ext_heat_trans_coef, self.int_heat_trans_coef, self.window_reflectance, 
                           self.ext_window_absorption, self.int_window_absorption)[5]    # Floor wall interior solar distribution fraction
        
        self.rad_part  = radiative_part  # Radiative part of internal gains
        self.conv_part = convective_part # Convective part of internal gains
        
    """RC properties (calculated externally via excel file - see documentation)."""
    @property
    def R_Window(self):         
        Rwind  = 1/(self.U_window*self.windows_area)         # Window thermal resistance (K/W)
        return (Rwind)
    
    @property
    def NatVent(self):       
        Rinf  = 1/(self.volume*(self.ACH/3600)*1.204*1012)   # Natural venilation thermal resistance (K/W)
        return (Rinf)
    @property
    def Kinfwind(self):
        return (1/self.NatVent)+(1/self.R_Window)            # Total thermal conductance due to ventilation (W/K)
    
    # Walls
    @property
    def Kwi(self):
        return 680    # Thermal conductance of the wall layers from the thermal insulation inward (W/K)
    @property
    def Kwin(self):
        return 914    # Thermal conductance of the thermal insulation layer (W/K)
    @property
    def Kwe(self):
        return 355   # Thermal conductance of the wall layers from the thermal insulation outward (W/K)
    @property
    def Cwi(self):
        return (22429600.00/3600)    # Thermal capacitance of inner layers (Wh/K)
    @property
    def Cwe(self):
        return  (69776827.76/3600)    # Thermal capacitance of the outer layers (Wh/K)
   
    # Roof
    @property
    def Kri(self):
        return 246    # Thermal conductance of the roof layers from the thermal insulation inward (W/K)
    @property
    def Krin(self):
        return 443     # Thermal conductance of the thermal insulation layer (W/K)
    @property
    def Kre(self):
        return 3323   # Thermal conductance of the roof layers from the thermal insulation outward (W/K)
    @property
    def Cri(self):
        return (69639600.00/3600)    # Thermal capacitance of inner layers (Wh/K)
    @property
    def Cre(self):
        return (76610.21/3600)      # Thermal capacitance of the outer layers (Wh/K)
    
    # Floor
    @property
    def Kfi(self):
        return 556      # Thermal conductance of the floor layers from the thermal insulation inward (W/K)
    @property
    def Kfin(self):
        return 162      # Thermal conductance of the thermal insulation layer (W/K)
    @property
    def Kfe(self):
        return 1099     # Thermal conductance of the floor layers from the thermal insulation outward (W/K)
    @property
    def Cfi(self):
        return (28257050.00/3600)    # Thermal capacitance of inner layers (Wh/K)
    @property
    def Cfe(self):
        return(55405272.64/3.6)     # Thermal capacitance of the outer layers (Wh/K)
    @property
    
    # Air
    def Cair(self):
        return self.volume*1.204*1012/3600  # Thermal capacitance of the air (Wh/K)

    def Uvalues(self):
        """Thermal resistance calculation. """
        R_wall   = ((1/self.Kwe)+(1/self.Kwi)+(1/self.Kwin))*self.total_opaque_area
        R_roof   = ((1/self.Kre)+(1/self.Kri)+(1/self.Krin))*self.roof_area
        R_floor  = ((1/self.Kfe)+(1/self.Kfi)+(1/self.Kfin))*self.floor_area
        """Thermal transmittance calculation. """
        self.Uwall  = 1/R_wall
        self.Uroof  = 1/R_roof
        self.Ufloor = 1/R_floor
        """Thermal loss. """
        self.Loss   = ((self.Uwall*self.total_opaque_area)+(self.Uroof*self.roof_area)+(self.Ufloor*self.floor_area)+(self.ACH*self.volume*1.204*1012/3600))

    def stateSpace_air_heating(self) :
        """State space model formulation for air heating system. """
        # Thermal node of the wall layers from the thermal insulation outward (Twe)
        A11 = -(self.Kwe+self.Kwin)/self.Cwe
        A12 = self.Kwin/self.Cwe
        A13 = 0.0
        A14 = 0.0
        A15 = 0.0
        A16 = 0.0
        A17 = 0.0 
               
        B11 = self.Kwe/self.Cwe
        B12 = 0.0
        B13 = 0.0
        B14 = 1/self.Cwe
        B15 = 0.0
        B16 = 0.0
        B17 = 0.0
       
        # Thermal node of the wall layers from the thermal insulation inward (Twi)
        A21 = self.Kwin/self.Cwi
        A22 = -(self.Kwi+self.Kwin)/self.Cwi
        A23 = 0.0
        A24 = 0.0
        A25 = 0.0
        A26 = 0.0
        A27 = self.Kwi/self.Cwi       

        B21 = 0.0
        B22 = 0.0
        B23 = (self.SF_south + self.SF_east + self.SF_west + self.SF_north)/self.Cwi   # solar gains
        B24 = 0.0
        B25 = 0.0
        B26 = (self.rad_part*(self.total_opaque_area/self.total_surf_area))/self.Cwi   # radiative part of internal gains
        B27 = 0.0

        # Thermal node of the roof layers from the thermal insulation outward (Tre)
        A31 = 0.0
        A32 = 0.0
        A33 = -(self.Kre+self.Krin)/self.Cre
        A34 = self.Krin/self.Cre
        A35 = 0.0
        A36 = 0.0
        A37 = 0.0
     
        B31 = self.Kre/self.Cre
        B32 = 0.0
        B33 = 0.0
        B34 = 0.0
        B35 = 1/self.Cre
        B36 = 0.0
        B37 = 0.0
        
        # Thermal node of the roof layers from the thermal insulation inward (Tri)
        A41 = 0.0
        A42 = 0.0
        A43 = self.Krin/self.Cri
        A44 = -(self.Kri+self.Krin)/self.Cri
        A45 = 0.0
        A46 = 0.0
        A47 = self.Kri/self.Cri
       
        B41 = 0.0
        B42 = 0.0
        B43 = self.SF_ceiling/self.Cri # solar gains
        B44 = 0.0
        B45 = 0.0
        B46 = (self.rad_part*(self.roof_area/self.total_surf_area))/self.Cri  # radiative part of internal gains
        B47 = 0.0

        # Thermal node of the floor layers from the thermal insulation outward (Tfe)
        A51 = 0.0
        A52 = 0.0
        A53 = 0.0
        A54 = 0.0
        A55 = -(self.Kfe+self.Kfin)/self.Cfe
        A56 = self.Kfin/self.Cfe
        A57 = 0.0
                
        B51 = 0.0
        B52 = self.Kfe/self.Cfe
        B53 = 0.0
        B54 = 0.0
        B55 = 0.0
        B56 = 0.0 
        B57 = 0.0
        
        # Thermal node of the floor layers from the thermal insulation inward (Tfi)
        A61 = 0.0
        A62 = 0.0
        A63 = 0.0
        A64 = 0.0
        A65 = self.Kfin/self.Cfi
        A66 = -(self.Kfi+self.Kfin)/self.Cfi
        A67 = self.Kfi/self.Cfi        
        
        B61 = 0.0
        B62 = 0.0
        B63 = self.SF_floor/self.Cfi # solar gains
        B64 = 0.0
        B65 = 0.0
        B66 = (self.rad_part*(self.floor_area/self.total_surf_area))/self.Cfi  # radiative part of internal gains
        B67 = 0.0
      
        # Thermal node of the air (Tair)
        A71 = 0.0
        A72 = self.Kwi/self.Cair
        A73 = 0.0
        A74 = self.Kri/self.Cair
        A75 = 0.0
        A76 = self.Kfi/self.Cair
        A77 = -(self.Kinfwind+self.Kwi+self.Kri+self.Kfi)/self.Cair
      
        B71 = (self.Kinfwind)/self.Cair
        B72 = 0.0
        B73 = 0.0
        B74 = 0.0
        B75 = 0.0
        B76 = self.conv_part/self.Cair   # convective part of internal gains
        B77 = 1/self.Cair                # thermal contribution (positive)
                    
        A = [[A11,A12,A13,A14,A15,A16,A17],[A21,A22,A23,A24,A25,A26,A27],[A31,A32,A33,A34,A35,A36,A37],[A41,A42,A43,A44,A45,A46,A47],[A51,A52,A53,A54,A55,A56,A57],[A61,A62,A63,A64,A65,A66,A67],[A71,A72,A73,A74,A75,A76,A77]]
        B = [[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47],[B51,B52,B53,B54,B55,B56,B57],[B61,B62,B63,B64,B65,B66,B67],[B71,B72,B73,B74,B75,B76,B77]]
        C = [[1.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,1.0]]
        D = [[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]]

        Ac =  np.array([[A11,A12,A13,A14,A15,A16,A17],[A21,A22,A23,A24,A25,A26,A27],[A31,A32,A33,A34,A35,A36,A37],[A41,A42,A43,A44,A45,A46,A47],[A51,A52,A53,A54,A55,A56,A57],[A61,A62,A63,A64,A65,A66,A67],[A71,A72,A73,A74,A75,A76,A77]])
        Bc =  np.array([[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47],[B51,B52,B53,B54,B55,B56,B57],[B61,B62,B63,B64,B65,B66,B67],[B71,B72,B73,B74,B75,B76,B77]])
        Cc =  np.array([[1.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,1.0]])
        Dc =  np.array([[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]])

        self.sys = signal.StateSpace(A, B, C, D)    
        self.sys_tuple = (Ac, Bc, Cc, Dc)

    def stateSpace_air_cooling(self) :
        """State space model formulation for air cooling system. """
        # Thermal node of the wall layers from the thermal insulation outward (Twe)
        A11 = -(self.Kwe+self.Kwin)/self.Cwe
        A12 = self.Kwin/self.Cwe
        A13 = 0.0
        A14 = 0.0
        A15 = 0.0
        A16 = 0.0
        A17 = 0.0 
               
        B11 = self.Kwe/self.Cwe
        B12 = 0.0
        B13 = 0.0
        B14 = 1/self.Cwe
        B15 = 0.0
        B16 = 0.0
        B17 = 0.0
       
        # Thermal node of the wall layers from the thermal insulation inward (Twi)
        A21 = self.Kwin/self.Cwi
        A22 = -(self.Kwi+self.Kwin)/self.Cwi
        A23 = 0.0
        A24 = 0.0
        A25 = 0.0
        A26 = 0.0
        A27 = self.Kwi/self.Cwi       
    
        B21 = 0.0
        B22 = 0.0
        B23 = (self.SF_south + self.SF_east + self.SF_west + self.SF_north)/self.Cwi   # solar gains
        B24 = 0.0
        B25 = 0.0
        B26 = (self.rad_part*(self.total_opaque_area/self.total_surf_area))/self.Cwi   # radiative part of internal gains
        B27 = 0.0
    
        # Thermal node of the roof layers from the thermal insulation outward (Tre)
        A31 = 0.0
        A32 = 0.0
        A33 = -(self.Kre+self.Krin)/self.Cre
        A34 = self.Krin/self.Cre
        A35 = 0.0
        A36 = 0.0
        A37 = 0.0
     
        B31 = self.Kre/self.Cre
        B32 = 0.0
        B33 = 0.0
        B34 = 0.0
        B35 = 1/self.Cre
        B36 = 0.0
        B37 = 0.0
        
        # Thermal node of the roof layers from the thermal insulation inward (Tri)
        A41 = 0.0
        A42 = 0.0
        A43 = self.Krin/self.Cri
        A44 = -(self.Kri+self.Krin)/self.Cri
        A45 = 0.0
        A46 = 0.0
        A47 = self.Kri/self.Cri
       
        B41 = 0.0
        B42 = 0.0
        B43 = self.SF_ceiling/self.Cri # solar gains
        B44 = 0.0
        B45 = 0.0
        B46 = (self.rad_part*(self.roof_area/self.total_surf_area))/self.Cri  # radiative part of internal gains
        B47 = 0.0
    
        # Thermal node of the floor layers from the thermal insulation outward (Tfe)
        A51 = 0.0
        A52 = 0.0
        A53 = 0.0
        A54 = 0.0
        A55 = -(self.Kfe+self.Kfin)/self.Cfe
        A56 = self.Kfin/self.Cfe
        A57 = 0.0
                
        B51 = 0.0
        B52 = self.Kfe/self.Cfe
        B53 = 0.0
        B54 = 0.0
        B55 = 0.0
        B56 = 0.0 
        B57 = 0.0
        
        # Thermal node of the floor layers from the thermal insulation inward (Tfi)
        A61 = 0.0
        A62 = 0.0
        A63 = 0.0
        A64 = 0.0
        A65 = self.Kfin/self.Cfi
        A66 = -(self.Kfi+self.Kfin)/self.Cfi
        A67 = self.Kfi/self.Cfi        
        
        B61 = 0.0
        B62 = 0.0
        B63 = self.SF_floor/self.Cfi # solar gains
        B64 = 0.0
        B65 = 0.0
        B66 = (self.rad_part*(self.floor_area/self.total_surf_area))/self.Cfi  # radiative part of internal gains  
        B67 = 0.0
      
        # Thermal node of the air (Tair)
        A71 = 0.0
        A72 = self.Kwi/self.Cair
        A73 = 0.0
        A74 = self.Kri/self.Cair
        A75 = 0.0
        A76 = self.Kfi/self.Cair
        A77 = -(self.Kinfwind+self.Kwi+self.Kri+self.Kfi)/self.Cair
      
        B71 = (self.Kinfwind)/self.Cair
        B72 = 0.0
        B73 = 0.0
        B74 = 0.0
        B75 = 0.0
        B76 = self.conv_part/self.Cair   # convective part of internal gains
        B77 = (- 1)/self.Cair            # thermal contribution (negative)
                    
        A = [[A11,A12,A13,A14,A15,A16,A17],[A21,A22,A23,A24,A25,A26,A27],[A31,A32,A33,A34,A35,A36,A37],[A41,A42,A43,A44,A45,A46,A47],[A51,A52,A53,A54,A55,A56,A57],[A61,A62,A63,A64,A65,A66,A67],[A71,A72,A73,A74,A75,A76,A77]]
        B = [[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47],[B51,B52,B53,B54,B55,B56,B57],[B61,B62,B63,B64,B65,B66,B67],[B71,B72,B73,B74,B75,B76,B77]]
        C = [[1.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,1.0]]
        D = [[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]]
    
        Ac =  np.array([[A11,A12,A13,A14,A15,A16,A17],[A21,A22,A23,A24,A25,A26,A27],[A31,A32,A33,A34,A35,A36,A37],[A41,A42,A43,A44,A45,A46,A47],[A51,A52,A53,A54,A55,A56,A57],[A61,A62,A63,A64,A65,A66,A67],[A71,A72,A73,A74,A75,A76,A77]])
        Bc =  np.array([[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47],[B51,B52,B53,B54,B55,B56,B57],[B61,B62,B63,B64,B65,B66,B67],[B71,B72,B73,B74,B75,B76,B77]])
        Cc =  np.array([[1.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,1.0]])
        Dc =  np.array([[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]])
    
        self.sys = signal.StateSpace(A, B, C, D)    
        self.sys_tuple = (Ac, Bc, Cc, Dc)
        
    def sys_dis(self, sys, dt):             
      (A, B, C, D, dt) = signal.cont2discrete(sys, dt, method='zoh', alpha=None)
      self.sys_disc_A = A
      self.sys_disc_B = B
      
    def sys_dis_forOp(self, sys, dt):            
      (A, B, C, D, dt) = signal.cont2discrete(sys, dt, method='zoh', alpha=None)      
      A_opt = A
      B_opt = np.delete(B,6,1)
      E_opt = B[:,6].reshape(A.shape[1],1)  
      self.sys_disc_A = A_opt
      self.sys_disc_B = B_opt   
      self.sys_disc_E = E_opt      
    
    def simulation_sys(self,sys,U,T,Xo): 
        t,y,x = signal.lsim(sys, U, T, Xo, interp=True)
        self.output = y
        
class SFH_60(object):
    """
    Single-Family House 1946-1960 class definition (without thermal insulation). 
        Only air emission system is available for the archetype (contribution applied directly on the thermal air node - Tair).
    """
    def __init__(self,ACH, wall_absorption, radiative_part, convective_part):
        """Definition of building properties. """
        self.ACH         = ACH                  # Air change per hours (1/hr)
        self.wall_absorption = wall_absorption  # Wall absorption
        self.floor_area  = 84.6    # Floor Area (m2)
        self.roof_area   = 97.6    # Roof Area (m2)
        self.volume      = 583.0   # Total volume (m3)
        
        self.wall_wide         = math.sqrt(self.floor_area)            # Building width (m)
        self.wall_long         = math.sqrt(self.floor_area)            # Building length (m)
        self.wall_high         = self.volume/self.floor_area           # Building height (m)
        
        self.windows_area = 20.3 # Total window area (m2)
        self.S_windows_area      = self.windows_area/4        # Total south-facing window area (m2)
        self.E_windows_area      = self.windows_area/4        # Total east-facing window area (m2)
        self.W_windows_area      = self.windows_area/4        # Total west-facing window area (m2)
        self.N_windows_area      = self.windows_area/4        # Total north-facing window area (m2)
        
        self.south_area    = self.wall_wide*self.wall_high - self.S_windows_area # Total south-facing wall area (m2)
        self.north_area    = self.wall_wide*self.wall_high - self.N_windows_area # Total north-facing wall area (m2)
        self.east_area     = self.wall_long*self.wall_high - self.E_windows_area # Total east-facing wall area (m2)
        self.west_area     = self.wall_long*self.wall_high - self.W_windows_area # Total west-facing wall area (m2)
        
        self.total_opaque_area = (self.wall_wide*self.wall_high)*2+ (self.wall_long*self.wall_high)*2 - (self.windows_area) # Total opaque Area (m2)
        self.total_surf_area   = self.floor_area + self.roof_area + self.total_opaque_area                                  # Total surface area (m2)
        
        """Window features of SFH 1946-1960 (single glazing). """
        self.U_window   = 4.9    # Window thermal transmittance (W/(m2K))
        self.g_window   = 0.837  # Window glass solar factor (%)
        self.emissivity  = 0.84  # Window emissivity
        self.int_window_absorption = 0.142 # Internal window absorption
        self.int_heat_trans_coef = 3.6 + 4.4*(self.emissivity)/0.837 # Internal heat transfer coefficient (W/Km2) - according to UNI EN ISO 10077/1
        self.window_reflectance = 0.206 # Window reflectance
        
        """Defining internal distribution fractions of transmitted solar radiation according to the methodology described in ANSI/ASHRAE Standard 140 (ANSI/ASHRAE, 2020). 
            Through these factors, the internal distribution of the solar gains transmitted through the glass surfaces is defined. """

        self.SF_south   = SF_single_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high, 
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, 
                           self.N_windows_area, self.wall_absorption, self.U_window, 
                           self.int_heat_trans_coef, self.window_reflectance, self.int_window_absorption)[0]    # South wall interior solar distribution fraction
        self.SF_east    = SF_single_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high, 
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, 
                           self.N_windows_area, self.wall_absorption, self.U_window, 
                           self.int_heat_trans_coef, self.window_reflectance, self.int_window_absorption)[1]     # East wall interior solar distribution fraction
        self.SF_west    = SF_single_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high, 
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, 
                           self.N_windows_area, self.wall_absorption, self.U_window, 
                           self.int_heat_trans_coef, self.window_reflectance, self.int_window_absorption)[2]     # West wall interior solar distribution fraction
        self.SF_north   = SF_single_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high, 
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, 
                           self.N_windows_area, self.wall_absorption, self.U_window, 
                           self.int_heat_trans_coef, self.window_reflectance, self.int_window_absorption)[3]    # North wall interior solar distribution fraction
        self.SF_ceiling = SF_single_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high, 
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, 
                           self.N_windows_area, self.wall_absorption, self.U_window, 
                           self.int_heat_trans_coef, self.window_reflectance, self.int_window_absorption)[4]  # Ceiling interior solar distribution fraction
        self.SF_floor   = SF_single_glazing(self.roof_area, self.wall_wide, self.wall_long, self.wall_high, 
                           self.S_windows_area, self.E_windows_area, self.W_windows_area, 
                           self.N_windows_area, self.wall_absorption, self.U_window, 
                           self.int_heat_trans_coef, self.window_reflectance, self.int_window_absorption)[5]    # Floor wall interior solar distribution fraction

        self.rad_part  = radiative_part  # Radiative part of internal gains
        self.conv_part = convective_part # Convective part of internal gains
        
    """RC properties (calculated externally via excel file - see documentation)."""
    @property
    def R_Window(self):         
        Rwind  = 1/(self.U_window*self.windows_area)         # Window thermal resistance (K/W)
        return (Rwind)
    
    @property
    def NatVent(self):       
        Rinf  = 1/(self.volume*(self.ACH/3600)*1.204*1012)   # Natural venilation thermal resistance (K/W)
        return (Rinf)
    @property
    def Kinfwind(self):
        return (1/self.NatVent)+(1/self.R_Window)            # Total thermal conductance due to ventilation (W/K)
    
    # Walls
    @property
    def Kwi(self):
        return 1787    # Thermal conductance of internal convection (W/K)
    @property
    def Kwe(self):
        return 425     # Wall thermal conductance (W/K)
    @property
    def Cw(self):
        return (151886240.00/3600)    # Wall thermal capacitance (Wh/K)
    
    # Roof
    @property
    def Kri(self):
        return 576    # Thermal conductance of internal convection (W/K)
    @property
    def Kre(self):
        return 343    # Roof thermal conductance (W/K)
    @property
    def Cr(self):
        return (27328000.00/3600)  # Wall thermal capacitance (Wh/K)   
    
    # floor
    @property
    def Kfi(self):
        return 846     # Thermal conductance of internal convection (W/K)
    @property
    def Kfe(self):
        return 216     # Floor thermal conductance (W/K)
    @property
    def Cf(self):
        return (87941700.00/3600)  # Floor thermal capacitance (Wh/K)   
    
    # Air
    @property
    def Cair(self):
        return self.volume*1.204*1012/3600  # Air thermal capacitance (Wh/K)   

    def Uvalues(self):
        """Thermal resistance calculation. """
        R_wall   = ((1/self.Kwe)+(1/self.Kwi))*self.total_opaque_area
        R_roof   = ((1/self.Kre)+(1/self.Kri))*self.roof_area
        R_floor  = ((1/self.Kfe)+(1/self.Kfi))*self.floor_area
        """Thermal transmittance calculation. """
        self.Uwall  = 1/R_wall
        self.Uroof  = 1/R_roof
        self.Ufloor = 1/R_floor
        """Thermal loss. """
        self.Loss   = ((self.Uwall*self.total_opaque_area)+(self.Uroof*self.roof_area)+(self.Ufloor*self.floor_area)+(self.ACH*self.volume*1.204*1012/3600))

    def stateSpace_air_heating(self) :
        """State space model formulation for air heating system. """
        # Thermal node of the wall (Tw)
        A11 = -(self.Kwe+self.Kwi)/self.Cw
        A12 = 0.0
        A13 = 0.0
        A14 = self.Kwi/self.Cw
         
        B11 = self.Kwe/self.Cw
        B12 = 0.0
        B13 = (self.SF_south + self.SF_east + self.SF_west + self.SF_north)/self.Cw # solar gains
        B14 = 1/self.Cw
        B15 = 0.0
        B16 = (self.rad_part*(self.total_opaque_area/self.total_surf_area))/self.Cw # radiative part of internal gains
        B17 = 0.0
       
        # Thermal node of the roof (Tr)
        A21 = 0.0
        A22 = -(self.Kre+self.Kri)/self.Cr
        A23 = 0.0
        A24 = self.Kri/self.Cr
     
        B21 = self.Kre/self.Cr
        B22 = 0.0
        B23 = self.SF_ceiling/self.Cr # solar gains
        B24 = 0.0
        B25 = 1/self.Cr
        B26 = (self.rad_part*(self.roof_area/self.total_surf_area))/self.Cr  # radiative part of internal gains
        B27 = 0.0

        # Thermal node of the floor (Tf)
        A31 = 0.0
        A32 = 0.0
        A33 = -(self.Kfe+self.Kfi)/self.Cf
        A34= self.Kfi/self.Cf
                
        B31 = 0.0
        B32 = self.Kfe/self.Cf
        B33 = self.SF_floor/self.Cf # solar gains
        B34 = 0.0
        B35 = 0.0
        B36 = (self.rad_part*(self.floor_area/self.total_surf_area))/self.Cf  # radiative part of internal gains   
        B37 = 0.0
        
        # Thermal node of the air (Tair)
        A41 = self.Kwi/self.Cair
        A42 = self.Kri/self.Cair
        A43 = self.Kfi/self.Cair
        A44 = -(self.Kinfwind+self.Kwi+self.Kri+self.Kfi)/self.Cair
      
        B41 = (self.Kinfwind)/self.Cair
        B42 = 0.0
        B43 = 0.0
        B44 = 0.0
        B45 = 0.0
        B46 = self.conv_part/self.Cair # convective part of internal gains
        B47 = 1/self.Cair              # thermal contribution (positive)
                    
        A = [[A11,A12,A13,A14],[A21,A22,A23,A24],[A31,A32,A33,A34],[A41,A42,A43,A44]]
        B = [[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47]]
        C = [[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]]
        D = [[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]]

        Ac =  np.array([[A11,A12,A13,A14],[A21,A22,A23,A24],[A31,A32,A33,A34],[A41,A42,A43,A44]])
        Bc =  np.array([[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47]])
        Cc =  np.array([[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]])
        Dc =  np.array([[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]])

        self.sys = signal.StateSpace(A, B, C, D)    
        self.sys_tuple = (Ac, Bc, Cc, Dc)
        
        
    def stateSpace_air_cooling(self) :
        """State space model formulation for air cooling system. """
        # Thermal node of the wall (Tw)
        A11 = -(self.Kwe+self.Kwi)/self.Cw
        A12 = 0.0
        A13 = 0.0
        A14 = self.Kwi/self.Cw
         
        B11 = self.Kwe/self.Cw
        B12 = 0.0
        B13 = (self.SF_south + self.SF_east + self.SF_west + self.SF_north)/self.Cw # solar gains
        B14 = 1/self.Cw
        B15 = 0.0
        B16 = (self.rad_part*(self.total_opaque_area/self.total_surf_area))/self.Cw # radiative part of internal gains
        B17 = 0.0
       
        # Thermal node of the roof (Tr)
        A21 = 0.0
        A22 = -(self.Kre+self.Kri)/self.Cr
        A23 = 0.0
        A24 = self.Kri/self.Cr
     
        B21 = self.Kre/self.Cr
        B22 = 0.0
        B23 = self.SF_ceiling/self.Cr # solar gains
        B24 = 0.0
        B25 = 1/self.Cr
        B26 = (self.rad_part*(self.roof_area/self.total_surf_area))/self.Cr  # radiative part of internal gains
        B27 = 0.0

        # Thermal node of the floor (Tf)
        A31 = 0.0
        A32 = 0.0
        A33 = -(self.Kfe+self.Kfi)/self.Cf
        A34= self.Kfi/self.Cf
                
        B31 = 0.0
        B32 = self.Kfe/self.Cf
        B33 = self.SF_floor/self.Cf # solar gains
        B34 = 0.0
        B35 = 0.0
        B36 = (self.rad_part*(self.floor_area/self.total_surf_area))/self.Cf  # radiative part of internal gains   
        B37 = 0.0
        
        # Thermal node of the air (Tair)
        A41 = self.Kwi/self.Cair
        A42 = self.Kri/self.Cair
        A43 = self.Kfi/self.Cair
        A44 = -(self.Kinfwind+self.Kwi+self.Kri+self.Kfi)/self.Cair
      
        B41 = (self.Kinfwind)/self.Cair
        B42 = 0.0
        B43 = 0.0
        B44 = 0.0
        B45 = 0.0
        B46 = self.conv_part/self.Cair   # convective part of internal gains
        B47 = (-1)/self.Cair             # thermal contribution (negative)
                    
        A = [[A11,A12,A13,A14],[A21,A22,A23,A24],[A31,A32,A33,A34],[A41,A42,A43,A44]]
        B = [[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47]]
        C = [[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]]
        D = [[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]]

        Ac =  np.array([[A11,A12,A13,A14],[A21,A22,A23,A24],[A31,A32,A33,A34],[A41,A42,A43,A44]])
        Bc =  np.array([[B11,B12,B13,B14,B15,B16,B17],[B21,B22,B23,B24,B25,B26,B27],[B31,B32,B33,B34,B35,B36,B37],[B41,B42,B43,B44,B45,B46,B47]])
        Cc =  np.array([[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]])
        Dc =  np.array([[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]])

        self.sys = signal.StateSpace(A, B, C, D)    
        self.sys_tuple = (Ac, Bc, Cc, Dc)

    def sys_dis(self, sys, dt):             
      (A, B, C, D, dt) = signal.cont2discrete(sys, dt, method='zoh', alpha=None)
      self.sys_disc_A = A
      self.sys_disc_B = B
      
    def sys_dis_forOp(self, sys, dt):            
      (A, B, C, D, dt) = signal.cont2discrete(sys, dt, method='zoh', alpha=None)      
      A_opt = A
      B_opt = np.delete(B,6,1)
      E_opt = B[:,6].reshape(A.shape[1],1)  
      self.sys_disc_A = A_opt
      self.sys_disc_B = B_opt   
      self.sys_disc_E = E_opt      
    
    def simulation_sys(self,sys,U,T,Xo): 
        t,y,x = signal.lsim(sys, U, T, Xo, interp=True)
        self.output = y