"""
ClustEnergy OpTool
Released on January, 2024
@author: DIISM UNIVPM
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
import pandas as pd
import numpy as np
from scipy.optimize import linprog
from Meteo import MeteoFile, PV_load
from Archetypes import SFH2006_air, SFH2005_air, SFH1990_air, SFH1960_air, SFH2006_floor
import HPperformance
import Archetypes
from User_pattern import intGains, occProfile
import User_pattern
import math

"""Functions"""
def base_load(Locality, day, month, duration, timestep, Thermal_load, 
                            n_SFH06_air, n_SFH05_air, n_SFH90_air, n_SFH60_air,n_SFH06_floor, 
                            HP_data_type, 
                            Supply_air, rated_cap_06_air, rated_cap_05_air, rated_cap_90_air, rated_cap_60_air, rated_COP_06_air, rated_COP_05_air, rated_COP_90_air, rated_COP_60_air,
                            Supply_floor, rated_cap_06_floor, rated_COP_06_floor,
                            Setpoint_type, setpoint, th_tolerance_BL_sup_air, th_tolerance_BL_min_air,th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor,
                            do_plot_E_bl_single, do_plot_E_bl_cluster, do_plot_Tair_bl_single, name_file):
    """base_load function defines thermal loads (cooling/heating) and electricity demand during the baseline scenario"""
    
    # number of buildings:
    n_SFH06_air = int(n_SFH06_air)
    n_SFH05_air = int(n_SFH05_air)
    n_SFH90_air = int(n_SFH90_air)
    n_SFH60_air = int(n_SFH60_air)
    
    if Thermal_load == "heating":
        n_SFH06_floor = int(n_SFH06_floor) # Radiant floor emission system for the SFH06 archetype (age class 2006-today) only during heating
    else: 
        n_SFH06_floor = 0
        
    n_Bui       = n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air + n_SFH06_floor # Total number of buildings composing the cluster
    
    samples = int(duration/timestep) # Calculation of total steps

    T_est = MeteoFile(Locality, day, month, duration, timestep)[0] # Outdoor air temperature (°C)

    """Call of the performance characteristics of heat pumps based on the definition of heat loads 
           and attribute definition for calling temperature setpoint functions"""
    if Thermal_load == 'cooling':
        """Call performance characteristics and capacity of heat pumps. """
        if n_SFH06_air > 0:
            COP_SFH06_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH06_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH05_air > 0:
            COP_SFH05_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH05_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH90_air > 0:
            COP_SFH90_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH90_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH60_air > 0:
            COP_SFH60_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH60_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[1] for i in range(0,samples)]
            
        """Thermostat mode setting in space cooling"""
        thermostat = 'th_cooling'
    elif Thermal_load == 'heating':
        if n_SFH06_air > 0:
            COP_SFH06_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH06_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH05_air > 0:
            COP_SFH05_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH05_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH90_air > 0:
            COP_SFH90_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH90_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH60_air > 0:
            COP_SFH60_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH60_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH06_floor > 0:
            COP_SFH06_floor = [HPperformance.HP_heating(T_est[i],Supply_floor,rated_cap_06_air,rated_COP_06_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH06_floor = [HPperformance.HP_heating(T_est[i],Supply_floor,rated_cap_06_air,rated_COP_06_air,HP_data_type)[1] for i in range(0,samples)]
        
        """Thermostat mode setting in space heating"""
        thermostat = 'th_heating'
    else:
        print('Thermal load not available')
    
    """Definition of occupancy profiles, relative internal gains (people, equipment and lighting) and setpoint temperature for each building 
           (useful for diversification of consumption loads within the cluster)"""
    occ   = {} # occupancy
    gains = {} # internal gains (W)
    tsp   = {} # thermostat set-points (°C)
    
    if Setpoint_type == "pre_defined":
        """Pre-defined temperature set-point definition
               For adding profiles or modifying existing ones, go to the User_pattern.py file and see the functions 
               th_heating (set-points for space heating) or th_cooling (set-points for space cooling)"""
        for i in range(n_SFH06_air):
            tsp["tspSFH06_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i]
        for i in range(n_SFH05_air):
            tsp["tspSFH05_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air]
        for i in range(n_SFH90_air):
            tsp["tspSFH90_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air]
        for i in range(n_SFH60_air):
            tsp["tspSFH60_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air]
        for i in range(n_SFH06_floor):
            tsp["tspSFH06_{0}_floor".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air]
    
    elif Setpoint_type == "user_defined":
        """User-defined temperature set-point definition"""
        t_setpoint = []
        for i in range (samples):
            t_setpoint.append(setpoint)
        for i in range(n_SFH06_air):
            tsp["tspSFH06_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH05_air):
            tsp["tspSFH05_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH90_air):
            tsp["tspSFH90_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH60_air):
            tsp["tspSFH60_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH06_floor):
            tsp["tspSFH06_{0}_floor".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
    
    """Occupancy patterns and related internal gain profiles definition"""
    for i in range(n_SFH06_air):
        occ["occupancySFH06_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i]
        gains["gainsSFH06_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i]
    for i in range(n_SFH05_air):
        occ["occupancySFH05_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air]
        gains["gainsSFH05_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air]
    for i in range(n_SFH90_air):
        occ["occupancySFH90_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air]
        gains["gainsSFH90_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air]
    for i in range(n_SFH60_air):
        occ["occupancySFH60_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air]
        gains["gainsSFH60_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air]
    for i in range(n_SFH06_floor):
        occ["occupancySFH06_{0}_floor".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air]
        gains["gainsSFH06_{0}_floor".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air]

    """In the following, thermal loads provided as input to the state-space model formulation are calculated via linear programming . 
           Specifically, the thermal load of each building is minimized to maintain the indoor temperature set-point.
           Linear programming solves problems of the following form: 
               minimize c @ x
               such that:
               A_ub @ x <= b_ub
               A_eq @ x == b_eq
               lb <= x <= ub
               
               For more information see the documentation of the python library scipy.optimize.linprog"""
           
    c      = {}     # Coefficients of the linear objective function to be minimized
    Aub    = {}     # Matrix of inequality constraints
    Bub    = {}     # Vector of inequality constraints
    bounds = {}     # Sequence of pairs (min, max) for each element in x, defining the minimum and maximum values of the decision variable
    
    Z      = {}     # Variable describing the state of the system containing the decision variable (e.g., thermal demand)
    M      = {}     # Variable describing the state of the indoor air temperture
    res_opt   = {}  # Optimization result - baseline scenario
    Qres_bl   = {}  # Values of the decision variable (thermal demand in Wh)
    Tair_bl   = {}  # Indoor air temeperature (°C) - baseline scenario
    Q_bl   = {}     # Thermal demand (kWh) - baseline scenario
    E_bl   = {}     # Electric demand (kWh) - baseline scenario
    upper   = {}    # Upper range thermostat set-point profile (°C) - baseline scenario
    lower   = {}    # Lower range thermostat set-point profile (°C) - baseline scenario
    
    """Calling archetype (SFHs) functions from the "Archetypes.py" module. """
    for i in range(n_SFH06_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH06 buildings (age class 2006-today) with air emission system."""
        c["c_SFH06_{0}".format(i+1)]             = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[0]
        Aub["Aub_SFH06_{0}".format(i+1)]         = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[1]
        Bub["Bub_SFH06_{0}".format(i+1)]         = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[2]
        bounds["bounds_SFH06_{0}".format(i+1)]   = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH06_{0}".format(i+1)]             = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[4]
        Z["Z_SFH06_{0}".format(i+1)]             = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[5]

        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH06_{0}".format(i+1)]     = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[7]
        lower["lower_SFH06_{0}".format(i+1)]     = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH06_{0}".format(i+1)] = linprog(c=c["c_SFH06_{0}".format(i+1)], A_ub=Aub["Aub_SFH06_{0}".format(i+1)], b_ub=Bub["Bub_SFH06_{0}".format(i+1)],bounds=bounds["bounds_SFH06_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 2006-today building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH06_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH06_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH06_{0}".format(i+1)] = (res_opt["res_opt_SFH06_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH06_{0}".format(i+1)] = list(M["M_SFH06_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}".format(i+1)], Qres_bl["Qres_bl_SFH06_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH06_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH06_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH06_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH06_{0}".format(i+1)][r]/COP_SFH06_air[r] for r in range(0, samples)]
    for i in range(n_SFH05_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH05 buildings (age class 1991-2005) with air emission system."""
        c["c_SFH05_{0}".format(i+1)]             = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[0]
        Aub["Aub_SFH05_{0}".format(i+1)]         = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[1]
        Bub["Bub_SFH05_{0}".format(i+1)]         = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[2]
        bounds["bounds_SFH05_{0}".format(i+1)]   = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH05_{0}".format(i+1)]             = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[4]
        Z["Z_SFH05_{0}".format(i+1)]             = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[5]

        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH05_{0}".format(i+1)] = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[7]
        lower["lower_SFH05_{0}".format(i+1)] = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH05_{0}".format(i+1)] = linprog(c=c["c_SFH05_{0}".format(i+1)], A_ub=Aub["Aub_SFH05_{0}".format(i+1)], b_ub=Bub["Bub_SFH05_{0}".format(i+1)],bounds=bounds["bounds_SFH05_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 1990-2005 building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH05_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH05_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH05_{0}".format(i+1)] = (res_opt["res_opt_SFH05_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH05_{0}".format(i+1)] = list(M["M_SFH05_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH05_{0}".format(i+1)], Qres_bl["Qres_bl_SFH05_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH05_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH05_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH05_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH05_{0}".format(i+1)][r]/COP_SFH05_air[r] for r in range(0, samples)]
    for i in range(n_SFH90_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH90 buildings (age class 1976-1990) with air emission system."""
        c["c_SFH90_{0}".format(i+1)]             = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[0]
        Aub["Aub_SFH90_{0}".format(i+1)]         = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[1]
        Bub["Bub_SFH90_{0}".format(i+1)]         = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[2]
        bounds["bounds_SFH90_{0}".format(i+1)]   = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH90_{0}".format(i+1)]             = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[4]
        Z["Z_SFH90_{0}".format(i+1)]             = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[5]

        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH90_{0}".format(i+1)]     = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[7]
        lower["lower_SFH90_{0}".format(i+1)]     = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[8]
        
        """Optimization result (baseline scenario)"""
        res_opt["res_opt_SFH90_{0}".format(i+1)] = linprog(c=c["c_SFH90_{0}".format(i+1)], A_ub=Aub["Aub_SFH90_{0}".format(i+1)], b_ub=Bub["Bub_SFH90_{0}".format(i+1)],bounds=bounds["bounds_SFH90_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 1976-1990 building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH90_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH90_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH90_{0}".format(i+1)] = (res_opt["res_opt_SFH90_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH90_{0}".format(i+1)] = list(M["M_SFH90_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH90_{0}".format(i+1)], Qres_bl["Qres_bl_SFH90_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH90_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH90_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH90_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH90_{0}".format(i+1)][r]/COP_SFH90_air[r] for r in range(0, samples)]
    for i in range(n_SFH60_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH60 buildings (age class 1946-1960) with air emission system."""
        c["c_SFH60_{0}".format(i+1)]             = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[0]
        Aub["Aub_SFH60_{0}".format(i+1)]         = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[1]
        Bub["Bub_SFH60_{0}".format(i+1)]         = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[2]
        bounds["bounds_SFH60_{0}".format(i+1)]   = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH60_{0}".format(i+1)]             = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[4]
        Z["Z_SFH60_{0}".format(i+1)]             = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[5]

        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH60_{0}".format(i+1)] = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[7]
        lower["lower_SFH60_{0}".format(i+1)] = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH60_{0}".format(i+1)] = linprog(c=c["c_SFH60_{0}".format(i+1)], A_ub=Aub["Aub_SFH60_{0}".format(i+1)], b_ub=Bub["Bub_SFH60_{0}".format(i+1)],bounds=bounds["bounds_SFH60_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 1946-1960 building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH60_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH60_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")

        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH60_{0}".format(i+1)] = (res_opt["res_opt_SFH60_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH60_{0}".format(i+1)] = list(M["M_SFH60_{0}".format(i+1)][3, :].reshape(samples, 1) + np.dot(Z["Z_SFH60_{0}".format(i+1)], Qres_bl["Qres_bl_SFH60_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH60_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH60_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH60_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH60_{0}".format(i+1)][r]/COP_SFH60_air[r] for r in range(0, samples)]
    for i in range(n_SFH06_floor):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH06 buildings (age class 2006-today) with underfloor heating system."""
        c["c_SFH06_{0}_floor".format(i+1)]             = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[0]
        Aub["Aub_SFH06_{0}_floor".format(i+1)]         = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[1]
        Bub["Bub_SFH06_{0}_floor".format(i+1)]         = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[2]
        bounds["bounds_SFH06_{0}_floor".format(i+1)]   = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[3]
        
        """Calling of variables referring to the state of indoor air temperature node. """
        M["M_SFH06_{0}_floor".format(i+1)]             = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[4]
        Z["Z_SFH06_{0}_floor".format(i+1)]             = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[5]

        """Calling of variables referring to the state of indoor air temperature. """
        upper["upper_SFH06_{0}_floor".format(i+1)]     = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[7]
        lower["lower_SFH06_{0}_floor".format(i+1)]     = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH06_{0}_floor".format(i+1)] = linprog(c=c["c_SFH06_{0}_floor".format(i+1)], A_ub=Aub["Aub_SFH06_{0}_floor".format(i+1)], b_ub=Bub["Bub_SFH06_{0}_floor".format(i+1)],bounds=bounds["bounds_SFH06_{0}_floor".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype building {0}".format(i+1) + " with underfloor heating system:")
        print(res_opt["res_opt_SFH06_{0}_floor".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH06_{0}_floor".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH06_{0}_floor".format(i+1)] = (res_opt["res_opt_SFH06_{0}_floor".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH06_{0}_floor".format(i+1)] = list(M["M_SFH06_{0}_floor".format(i+1)][7, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}_floor".format(i+1)], Qres_bl["Qres_bl_SFH06_{0}_floor".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH06_{0}_floor".format(i+1)]       = [Qres_bl["Qres_bl_SFH06_{0}_floor".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH06_{0}_floor".format(i+1)]       = [Q_bl["Q_bl_SFH06_{0}_floor".format(i+1)][r]/COP_SFH06_floor[r] for r in range(0, samples)]
    
    """Conversion of dictionary to data frame for calculation of aggregate electric power consumption. """
    E_bl = pd.DataFrame.from_dict(E_bl)
    
    """Calculating base electrical loads at the aggregate level. """
    Etot_bl = E_bl.sum(axis=1)
    
    """Do plot"""
    import matplotlib.pyplot as plt
    
    time = [i*timestep for i in range(0, samples)] # time array
    
    if do_plot_E_bl_cluster:
        """If True, plot cluster electric consumption trends. """
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Time (hours)', fontsize=14)
        ax1.set_ylabel('$\.P_\mathrm{cluster}$ (kW$_\mathrm{el}$)', fontsize=14)
        ax1.plot(time, Etot_bl, 'k', label="BL_scenario", linewidth=1.0)
        ax1.set_xlim(0, duration-1)
        ax1.grid(alpha=0.3)
        plt.legend(frameon=False, fontsize=12)
        fig.set_size_inches(w=10, h=3)
        fig.tight_layout()
        plt.savefig("Eele_BL", dpi=200)
        
        plt.show()
        plt.close('all')
    
    if do_plot_E_bl_single:
        """If True, plot single buildings electric consumption trends. """
        for i in range (n_SFH06_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Thermal power (kW)')
            ax1.plot(time,Q_bl["Q_bl_SFH06_{0}".format(i+1)],'k--', linewidth=1, label = "Q_bl_SFH06_{0}".format(i+1))
            
            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH05_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Thermal power (kW)')
            ax1.plot(time,Q_bl["Q_bl_SFH05_{0}".format(i+1)],'k--', linewidth=1, label = "Q_bl_SFH05_{0}".format(i+1))
            
            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH90_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Thermal power (kW)')
            ax1.plot(time,Q_bl["Q_bl_SFH90_{0}".format(i+1)],'k--', linewidth=1, label = "Q_bl_SFH90_{0}".format(i+1))
            
            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH60_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Thermal power (kW)')
            ax1.plot(time,Q_bl["Q_bl_SFH60_{0}".format(i+1)],'k--', linewidth=1, label = "Q_bl_SFH60_{0}".format(i+1)) 
            
            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH06_floor):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Thermal power (kW)')
            ax1.plot(time,Q_bl["Q_bl_SFH06_{0}_floor".format(i+1)],'k--', linewidth=1, label = "Q_bl_SFH06_{0}_floor".format(i+1))
            
            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
        

    if do_plot_Tair_bl_single:
        """If True, plot single buildings indoor air temperature trends. """
        for i in range (n_SFH06_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Temperature (°C)')
            ax1.plot(time,T_est,'r', linewidth=1, label = "T_outdoor")
            ax1.plot(time,Tair_bl["Tair_bl_SFH06_{0}".format(i+1)],'b', linewidth=1, label = "Tair_bl_SFH06_{0}".format(i+1))
            ax1.plot(time,upper["upper_SFH06_{0}".format(i+1)],'c', linewidth=1, label = "upper_SFH06_{0}".format(i+1))
            ax1.plot(time,lower["lower_SFH06_{0}".format(i+1)],'c--', linewidth=1, label = "lower_SFH06_{0}".format(i+1))
            
            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH05_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Temperature (°C)')
            ax1.plot(time,T_est,'r', linewidth=1, label = "T_outdoor")
            ax1.plot(time,Tair_bl["Tair_bl_SFH05_{0}".format(i+1)],'b', linewidth=1, label = "Tair_bl_SFH05_{0}".format(i+1))
            ax1.plot(time,upper["upper_SFH05_{0}".format(i+1)],'c', linewidth=1, label = "upper_SFH05_{0}".format(i+1))
            ax1.plot(time,lower["lower_SFH05_{0}".format(i+1)],'c--', linewidth=1, label = "lower_SFH05_{0}".format(i+1))
            
            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH90_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Temperature (°C)')  # we already handled the x-label with ax1
            ax1.plot(time,T_est,'r', linewidth=1, label = "T_outdoor")
            ax1.plot(time,Tair_bl["Tair_bl_SFH90_{0}".format(i+1)],'b', linewidth=1, label = "Tair_bl_SFH90_{0}".format(i+1))
            ax1.plot(time,upper["upper_SFH90_{0}".format(i+1)],'c', linewidth=1, label = "upper_SFH90_{0}".format(i+1))
            ax1.plot(time,lower["lower_SFH90_{0}".format(i+1)],'c--', linewidth=1, label = "lower_SFH90_{0}".format(i+1))
            
            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH60_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Temperature (°C)')  # we already handled the x-label with ax1
            ax1.plot(time,T_est,'r', linewidth=1, label = "T_outdoor")
            ax1.plot(time,Tair_bl["Tair_bl_SFH60_{0}".format(i+1)],'b', linewidth=1, label = "Tair_bl_SFH60_{0}".format(i+1))
            ax1.plot(time,upper["upper_SFH60_{0}".format(i+1)],'c', linewidth=1, label = "upper_SFH60_{0}".format(i+1))
            ax1.plot(time,lower["lower_SFH60_{0}".format(i+1)],'c--', linewidth=1, label = "lower_SFH60_{0}".format(i+1))
            
            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH06_floor):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Temperature (°C)')  # we already handled the x-label with ax1
            ax1.plot(time,T_est,'r', linewidth=1, label = "T_outdoor")
            ax1.plot(time,Tair_bl["Tair_bl_SFH06_{0}_floor".format(i+1)],'b', linewidth=1, label = "Tair_bl_SFH06_{0}_floor".format(i+1))
            ax1.plot(time,upper["upper_SFH06_{0}_floor".format(i+1)],'c', linewidth=1, label = "upper_SFH06_{0}_floor".format(i+1))
            ax1.plot(time,lower["lower_SFH06_{0}_floor".format(i+1)],'c--', linewidth=1, label = "lower_SFH06_{0}_floor".format(i+1))
            
            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')


    return(E_bl)

def peak_time(Locality, day, month, duration, timestep, Etot_bl):
    """The peak_time function calculates the time when the peak load occurs and its intensity. """
    
    samples = int(duration/timestep)
    n_days = int(duration/24)

    Max         = []
    time_max    = []
    
    Etot_bl = [Etot_bl[i] for i in range(0,samples)]  # (kWh)
    
    for n in range (n_days):
        Max.append(np.amax(np.array(Etot_bl[int(samples/n_days)*n:int(samples/n_days)*(n + 1)])))
        time_max.append(Etot_bl.index(Max[n]))           # time when maximum consumption occurs
    
    for n in range (n_days):
        print("\nMaximum thermal power during day",n,"of about " +
              str(round(Max[n],2))+" kWth at = "+str(time_max[n]*timestep) + " hr")
        print("\n------------------------------------------------------------------------")
        
    return(time_max)

def peak_shaving(Locality, day, month, duration, timestep, Thermal_load, 
                            n_SFH06_air, n_SFH05_air, n_SFH90_air, n_SFH60_air,n_SFH06_floor, 
                            HP_data_type, 
                            Supply_air, rated_cap_06_air, rated_cap_05_air, rated_cap_90_air, rated_cap_60_air, rated_COP_06_air, rated_COP_05_air, rated_COP_90_air, rated_COP_60_air,
                            Supply_floor, rated_cap_06_floor, rated_COP_06_floor,
                            Setpoint_type, setpoint, th_tolerance_BL_sup_air, th_tolerance_BL_min_air, th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor,
                            th_tolerance_DR_sup_air, th_tolerance_DR_min_air, th_tolerance_DR_sup_floor, th_tolerance_DR_min_floor, 
                            f_red, duration_dr, f_limit,
                            do_plot_E_dr_single, do_plot_E_dr_cluster,do_plot_Tair_dr_single, name_file):
    
    """peak_shaving is a function based on linear programming that minimizes the aggregate electrical consumption of a 
            user-defined cluster of buildings when applying a peak load reduction strategy. """
    
    duration_dr_h = duration_dr
    n_days        = int(duration/24)
    
    # number of buildings:
    n_SFH06_air = int(n_SFH06_air)
    n_SFH05_air = int(n_SFH05_air)
    n_SFH90_air = int(n_SFH90_air)
    n_SFH60_air = int(n_SFH60_air)
    
    if Thermal_load == "heating":
        n_SFH06_floor = int(n_SFH06_floor) # Radiant floor emission system for the SFH06 archetype (age class 2006-today) only during heating
    else: 
        n_SFH06_floor = 0
    
    n_Bui       = n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air + n_SFH06_floor # Total number of buildings composing the cluster
    
    samples = int(duration/timestep) # Calculation of total steps

    T_est = MeteoFile(Locality, day, month, duration, timestep)[0] # Outdoor air temperature (°C)

    """Call of the performance characteristics of heat pumps based on the definition of heat loads 
           and attribute definition for calling temperature setpoint functions"""
    if Thermal_load == 'cooling':
        """Call performance characteristics and capacity of heat pumps. """
        if n_SFH06_air > 0:
            COP_SFH06_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH06_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH05_air > 0:
            COP_SFH05_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH05_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH90_air > 0:
            COP_SFH90_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH90_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH60_air > 0:
            COP_SFH60_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH60_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[1] for i in range(0,samples)]
            
        """Thermostat mode setting in space cooling"""
        thermostat = 'th_cooling'
    elif Thermal_load == 'heating':
        if n_SFH06_air > 0:
            COP_SFH06_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH06_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH05_air > 0:
            COP_SFH05_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH05_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH90_air > 0:
            COP_SFH90_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH90_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH60_air > 0:
            COP_SFH60_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH60_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH06_floor > 0:
            COP_SFH06_floor = [HPperformance.HP_heating(T_est[i],Supply_floor,rated_cap_06_air,rated_COP_06_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH06_floor = [HPperformance.HP_heating(T_est[i],Supply_floor,rated_cap_06_air,rated_COP_06_air,HP_data_type)[1] for i in range(0,samples)]
        
        """Thermostat mode setting in space heating"""
        thermostat = 'th_heating'
    else:
        print('Thermal load not available')
    
    """Definition of occupancy profiles, relative internal gains (people, equipment and lighting) and setpoint temperature for each building 
           (useful for diversification of consumption loads within the cluster)"""
    occ   = {} # occupancy
    gains = {} # internal gains (W)
    tsp   = {} # thermostat set-points (°C)
    
    if Setpoint_type == "pre_defined":
        """Pre-defined temperature set-point definition
               For adding profiles or modifying existing ones, go to the User_pattern.py file and see the functions 
               th_heating (set-points for space heating) or th_cooling (set-points for space cooling)"""
        for i in range(n_SFH06_air):
            tsp["tspSFH06_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i]
        for i in range(n_SFH05_air):
            tsp["tspSFH05_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air]
        for i in range(n_SFH90_air):
            tsp["tspSFH90_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air]
        for i in range(n_SFH60_air):
            tsp["tspSFH60_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air]
        for i in range(n_SFH06_floor):
            tsp["tspSFH06_{0}_floor".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air]
    
    elif Setpoint_type == "user_defined":
        """User-defined temperature set-point definition"""
        t_setpoint = []
        for i in range (samples):
            t_setpoint.append(setpoint)
        for i in range(n_SFH06_air):
            tsp["tspSFH06_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH05_air):
            tsp["tspSFH05_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH90_air):
            tsp["tspSFH90_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH60_air):
            tsp["tspSFH60_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH06_floor):
            tsp["tspSFH06_{0}_floor".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
    
    """Occupancy patterns and related internal gain profiles definition"""
    for i in range(n_SFH06_air):
        occ["occupancySFH06_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i]
        gains["gainsSFH06_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i]
    for i in range(n_SFH05_air):
        occ["occupancySFH05_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air]
        gains["gainsSFH05_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air]
    for i in range(n_SFH90_air):
        occ["occupancySFH90_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air]
        gains["gainsSFH90_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air]
    for i in range(n_SFH60_air):
        occ["occupancySFH60_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air]
        gains["gainsSFH60_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air]
    for i in range(n_SFH06_floor):
        occ["occupancySFH06_{0}_floor".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air]
        gains["gainsSFH06_{0}_floor".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air]

    """In the following, thermal loads provided as input to the state-space model formulation are calculated via linear programming . 
           Specifically, the thermal load of each building is minimized to maintain the indoor temperature set-point.
           Linear programming solves problems of the following form: 
               minimize c @ x
               such that:
               A_ub @ x <= b_ub
               A_eq @ x == b_eq
               lb <= x <= ub
               
               For more information see the documentation of the python library scipy.optimize.linprog"""
           
    c      = {}     # Coefficients of the linear objective function to be minimized
    Aub    = {}     # Matrix of inequality constraints
    Bub    = {}     # Vector of inequality constraints
    bounds = {}     # Sequence of pairs (min, max) for each element in x, defining the minimum and maximum values of the decision variable
    
    Z      = {}     # Variable describing the state of the system containing the decision variable (e.g., thermal demand)
    M      = {}     # Variable describing the state of the indoor air temperture
    Mmax      = {}  # Variable describing the state of the indoor air temperture
    res_opt   = {}  # Optimization result - baseline scenario
    Qres_bl   = {}  # Values of the decision variable (thermal demand in Wh)
    Tair_bl   = {}  # Indoor air temeperature (°C) - baseline scenario
    Q_bl   = {}     # Thermal demand (kWh) - baseline scenario
    E_bl   = {}     # Electric demand (kWh) - baseline scenario
    upper   = {}    # Upper range thermostat set-point profile (°C) - baseline scenario
    lower   = {}    # Lower range thermostat set-point profile (°C) - baseline scenario
    
    Mmax_floor = {} # Variable describing the state of the floor temperture
    Z_floor    = {} # Variable describing the state of the floor system containing the decision variable (e.g., thermal demand)
    
    """Calling archetype (SFHs) functions from the "Archetypes.py" module. """
    for i in range(n_SFH06_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH06 buildings (age class 2006-today) with air emission system."""
        c["c_SFH06_{0}".format(i+1)]             = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[0]
        Aub["Aub_SFH06_{0}".format(i+1)]         = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[1]
        Bub["Bub_SFH06_{0}".format(i+1)]         = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[2]
        bounds["bounds_SFH06_{0}".format(i+1)]   = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH06_{0}".format(i+1)]             = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[4]
        Z["Z_SFH06_{0}".format(i+1)]             = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[5]
        Mmax["Mmax_SFH06_{0}".format(i+1)]       = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH06_{0}".format(i+1)]     = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[7]
        lower["lower_SFH06_{0}".format(i+1)]     = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH06_{0}".format(i+1)] = linprog(c=c["c_SFH06_{0}".format(i+1)], A_ub=Aub["Aub_SFH06_{0}".format(i+1)], b_ub=Bub["Bub_SFH06_{0}".format(i+1)],bounds=bounds["bounds_SFH06_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 2006-today building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH06_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH06_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH06_{0}".format(i+1)] = (res_opt["res_opt_SFH06_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH06_{0}".format(i+1)] = list(M["M_SFH06_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}".format(i+1)], Qres_bl["Qres_bl_SFH06_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH06_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH06_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH06_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH06_{0}".format(i+1)][r]/COP_SFH06_air[r] for r in range(0, samples)]
    for i in range(n_SFH05_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH05 buildings (age class 1991-2005) with air emission system."""
        c["c_SFH05_{0}".format(i+1)]             = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[0]
        Aub["Aub_SFH05_{0}".format(i+1)]         = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[1]
        Bub["Bub_SFH05_{0}".format(i+1)]         = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[2]
        bounds["bounds_SFH05_{0}".format(i+1)]   = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH05_{0}".format(i+1)]             = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[4]
        Z["Z_SFH05_{0}".format(i+1)]             = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[5]
        Mmax["Mmax_SFH05_{0}".format(i+1)]       = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH05_{0}".format(i+1)] = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[7]
        lower["lower_SFH05_{0}".format(i+1)] = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH05_{0}".format(i+1)] = linprog(c=c["c_SFH05_{0}".format(i+1)], A_ub=Aub["Aub_SFH05_{0}".format(i+1)], b_ub=Bub["Bub_SFH05_{0}".format(i+1)],bounds=bounds["bounds_SFH05_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 1990-2005 building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH05_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH05_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH05_{0}".format(i+1)] = (res_opt["res_opt_SFH05_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH05_{0}".format(i+1)] = list(M["M_SFH05_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH05_{0}".format(i+1)], Qres_bl["Qres_bl_SFH05_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH05_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH05_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH05_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH05_{0}".format(i+1)][r]/COP_SFH05_air[r] for r in range(0, samples)]
    for i in range(n_SFH90_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH90 buildings (age class 1976-1990) with air emission system."""
        c["c_SFH90_{0}".format(i+1)]             = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[0]
        Aub["Aub_SFH90_{0}".format(i+1)]         = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[1]
        Bub["Bub_SFH90_{0}".format(i+1)]         = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[2]
        bounds["bounds_SFH90_{0}".format(i+1)]   = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH90_{0}".format(i+1)]             = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[4]
        Z["Z_SFH90_{0}".format(i+1)]             = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[5]
        Mmax["Mmax_SFH90_{0}".format(i+1)]       = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH90_{0}".format(i+1)]     = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[7]
        lower["lower_SFH90_{0}".format(i+1)]     = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[8]
        
        """Optimization result (baseline scenario)"""
        res_opt["res_opt_SFH90_{0}".format(i+1)] = linprog(c=c["c_SFH90_{0}".format(i+1)], A_ub=Aub["Aub_SFH90_{0}".format(i+1)], b_ub=Bub["Bub_SFH90_{0}".format(i+1)],bounds=bounds["bounds_SFH90_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 1976-1990 building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH90_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH90_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH90_{0}".format(i+1)] = (res_opt["res_opt_SFH90_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH90_{0}".format(i+1)] = list(M["M_SFH90_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH90_{0}".format(i+1)], Qres_bl["Qres_bl_SFH90_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH90_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH90_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH90_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH90_{0}".format(i+1)][r]/COP_SFH90_air[r] for r in range(0, samples)]
    for i in range(n_SFH60_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH60 buildings (age class 1946-1960) with air emission system."""
        c["c_SFH60_{0}".format(i+1)]             = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[0]
        Aub["Aub_SFH60_{0}".format(i+1)]         = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[1]
        Bub["Bub_SFH60_{0}".format(i+1)]         = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[2]
        bounds["bounds_SFH60_{0}".format(i+1)]   = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH60_{0}".format(i+1)]             = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[4]
        Z["Z_SFH60_{0}".format(i+1)]             = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[5]
        Mmax["Mmax_SFH60_{0}".format(i+1)]       = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH60_{0}".format(i+1)] = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[7]
        lower["lower_SFH60_{0}".format(i+1)] = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH60_{0}".format(i+1)] = linprog(c=c["c_SFH60_{0}".format(i+1)], A_ub=Aub["Aub_SFH60_{0}".format(i+1)], b_ub=Bub["Bub_SFH60_{0}".format(i+1)],bounds=bounds["bounds_SFH60_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 1946-1960 building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH60_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH60_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")

        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH60_{0}".format(i+1)] = (res_opt["res_opt_SFH60_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH60_{0}".format(i+1)] = list(M["M_SFH60_{0}".format(i+1)][3, :].reshape(samples, 1) + np.dot(Z["Z_SFH60_{0}".format(i+1)], Qres_bl["Qres_bl_SFH60_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH60_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH60_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH60_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH60_{0}".format(i+1)][r]/COP_SFH60_air[r] for r in range(0, samples)]
    for i in range(n_SFH06_floor):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH06 buildings (age class 2006-today) with underfloor heating system."""
        c["c_SFH06_{0}_floor".format(i+1)]             = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[0]
        Aub["Aub_SFH06_{0}_floor".format(i+1)]         = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[1]
        Bub["Bub_SFH06_{0}_floor".format(i+1)]         = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[2]
        bounds["bounds_SFH06_{0}_floor".format(i+1)]   = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[3]
        
        """Calling of variables referring to the state of indoor air temperature node. """
        M["M_SFH06_{0}_floor".format(i+1)]             = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[4]
        Z["Z_SFH06_{0}_floor".format(i+1)]             = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[5]
        Mmax["Mmax_SFH06_{0}_floor".format(i+1)]       = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[6]
        
        """Calling of variables referring to the state of floor node. """
        Z_floor["Z_SFH06_{0}_f".format(i+1)]           = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[9]
        Mmax_floor["Mmax_SFH06_{0}_f".format(i+1)]     = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[10]
        
        """Calling of variables referring to the state of indoor air temperature. """
        upper["upper_SFH06_{0}_floor".format(i+1)]     = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[7]
        lower["lower_SFH06_{0}_floor".format(i+1)]     = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH06_{0}_floor".format(i+1)] = linprog(c=c["c_SFH06_{0}_floor".format(i+1)], A_ub=Aub["Aub_SFH06_{0}_floor".format(i+1)], b_ub=Bub["Bub_SFH06_{0}_floor".format(i+1)],bounds=bounds["bounds_SFH06_{0}_floor".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype building {0}".format(i+1) + " with underfloor heating system:")
        print(res_opt["res_opt_SFH06_{0}_floor".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH06_{0}_floor".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH06_{0}_floor".format(i+1)] = (res_opt["res_opt_SFH06_{0}_floor".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH06_{0}_floor".format(i+1)] = list(M["M_SFH06_{0}_floor".format(i+1)][7, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}_floor".format(i+1)], Qres_bl["Qres_bl_SFH06_{0}_floor".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH06_{0}_floor".format(i+1)]       = [Qres_bl["Qres_bl_SFH06_{0}_floor".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH06_{0}_floor".format(i+1)]       = [Q_bl["Q_bl_SFH06_{0}_floor".format(i+1)][r]/COP_SFH06_floor[r] for r in range(0, samples)]
    
    """Conversion of dictionary to data frame for calculation of aggregate electric power consumption. """
    E_bl = pd.DataFrame.from_dict(E_bl)
    
    """Calculating base electrical loads at the aggregate level. """
    Etot_bl = E_bl.sum(axis=1)
    Ebl_sum = float(Etot_bl.sum())
    Etot_bl = [Etot_bl[i] for i in range(0,samples)]  
    
    """------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    The Demand-Response (DR) event is modeled below. 
    From the definition of the PV generation profile and electricity demand during the reference scenario, 
    the electricity demand during the DR scenario is managed by shifting loads to centralized generation hours.
    To this end, an additional optimization constraint is introduced to minimize the aggregated energy demand drawn from the grid. 
    In particular, the electricity demand of the cluster of buildings (Eele_max) is constrained to remain below a certain 
    percentage reduction (1 - f_red) for the duration of the peak-shaving event (duration_dr). 
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""

    """For each simulation day, calculation of the aggregate peak load and the time at which it occurs. """
    duration_dr = int(duration_dr_h/timestep - 1)         # steps of peak shaving event duration
    Max         = []
    time_max    = []
    
    Max_dr      = []
    t_in_dr     = []
    t_end_dr    = []
    
    for n in range (n_days):
        Max.append(np.amax(np.array(Etot_bl[int(samples/n_days)*n:int(samples/n_days)*(n + 1)])))
        time_max.append(Etot_bl.index(Max[n]))            # Calculation of the critical step at which the aggregate peak load occurs
        for i in range (int(samples/n_days)*n,int(samples/n_days)*(n + 1)):
            Max_dr.append(Max[n])
            t_in_dr.append(time_max[n])                   # peak shaving event initiation step
            t_end_dr.append(time_max[n] + duration_dr)    # end step peak shaving event
        
    for n in range (n_days):
        print("\nMaximum thermal power during day",n,"of about " +
              str(round(Max[n],2))+" kWth at = "+str(time_max[n]*timestep) + " hr")
        print("\n------------------------------------------------------------------------")
    
    """Definition of the electricity demand of the cluster of buildings constrained to remain below a certain percentage reduction. """
    Eele_max = []

    for t in range(samples):
        if t_in_dr[t] <= t <= t_end_dr[t]:
            Eele_max.append(float(Etot_bl[t]*1000)*(1 - float(f_red))) # Aggregate electricity demand constrained during the peak shaving event (Pele <= (1-f_red)*Pele,max)
        else:
            Eele_max.append(float(Etot_bl[t]*1000)*(float(f_limit)))
    
    """Definition of temperature profiles during the DR event to unlock building energy flexibility. """
    if Thermal_load == 'cooling':
        for i in range(n_SFH06_air):
            for t in range(samples):
                if t_in_dr[t] <= t <= t_end_dr[t]:
                    upper["upper_SFH06_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air)  # Increase in indoor air temperature set point during the peak shaving event for cooling
                else:
                    lower["lower_SFH06_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air)  # Reduction of indoor air temperature set point outside the peak shaving event for cooling
        for i in range(n_SFH05_air):
            for t in range(samples):
                if t_in_dr[t] <= t <= t_end_dr[t]:
                    upper["upper_SFH05_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air)  # Increase in indoor air temperature set point during the peak shaving event for cooling
                else:
                    lower["lower_SFH05_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air)  # Reduction of indoor air temperature set point outside the peak shaving event for cooling
        for i in range(n_SFH90_air):
            for t in range(samples):
                if t_in_dr[t] <= t <= t_end_dr[t]:
                    upper["upper_SFH90_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air)  # Increase in indoor air temperature set point during the peak shaving event for cooling
                else:
                    lower["lower_SFH90_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air)  # Reduction of indoor air temperature set point outside the peak shaving event for cooling
        for i in range(n_SFH60_air):
            for t in range(samples):
                if t_in_dr[t] <= t <= t_end_dr[t]:
                    upper["upper_SFH60_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air)  # Increase in indoor air temperature set point during the peak shaving event for cooling
                else:
                    lower["lower_SFH60_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air)  # Reduction of indoor air temperature set point outside the peak shaving event for cooling
    elif Thermal_load == 'heating':
        for i in range(n_SFH06_air):
            for t in range(samples):
                if t_in_dr[t] <= t <= t_end_dr[t]:
                    lower["lower_SFH06_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air)  # Reduction of indoor air temperature set point during the peak shaving event for heating
                else:
                    upper["upper_SFH06_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air)  # Increase in indoor air temperature set point outside the peak shaving event for heating
        for i in range(n_SFH05_air):
            for t in range(samples):
                if t_in_dr[t] <= t <= t_end_dr[t]:
                    lower["lower_SFH05_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air)  # Reduction of indoor air temperature set point during the peak shaving event for heating
                else:
                    upper["upper_SFH05_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air)  # Increase in indoor air temperature set point outside the peak shaving event for heating
        for i in range(n_SFH90_air):
            for t in range(samples):
                if t_in_dr[t] <= t <= t_end_dr[t]:
                    lower["lower_SFH90_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air)  # Reduction of indoor air temperature set point during the peak shaving event for heating
                else:
                    upper["upper_SFH90_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air)  # Increase in indoor air temperature set point outside the peak shaving event for heating
        for i in range(n_SFH60_air):
            for t in range(samples):
                if t_in_dr[t] <= t <= t_end_dr[t]:
                    lower["lower_SFH60_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air)  # Reduction of indoor air temperature set point during the peak shaving event for heating
                else:
                    upper["upper_SFH60_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air)  # Increase in indoor air temperature set point outside the peak shaving event for heating
        for i in range(n_SFH06_floor):
            for t in range(samples):
                if t_in_dr[t] <= t <= t_end_dr[t]:
                    lower["lower_SFH06_{0}_floor".format(i+1)][t] -= float(th_tolerance_DR_min_air)  # Reduction of indoor air temperature set point during the peak shaving event for heating
                else:
                    upper["upper_SFH06_{0}_floor".format(i+1)][t] += float(th_tolerance_DR_sup_floor)  # Increase in indoor air temperature set point outside the peak shaving event for heating

    """In the following, thermal loads provided as input to the state-space model formulation are calculated via linear programming. 
           Specifically, the thermal load of each building is minimized to maintain the indoor temperature set-point, 
           and the electricity demand of the cluster of buildings is constrained to remain below a certain 
           percentage reduction (Eele_max) for the duration of the peak-shaving event (duration_dr).
           Linear programming solves problems of the following form: 
               minimize c @ x
               such that:
               A_ub @ x <= b_ub
               A_eq @ x == b_eq
               lb <= x <= ub
               
               For more information see the documentation of the python library scipy.optimize.linprog"""

    """Definition of the number of constraints at the individual building and cluster level. """
    Nc_06_air   = (2)*n_SFH06_air          # 2 for indoor air temperature range 
    Nc_05_air   = (2)*n_SFH05_air          # 2 for indoor air temperature range 
    Nc_90_air   = (2)*n_SFH90_air          # 2 for indoor air temperature range 
    Nc_60_air   = (2)*n_SFH60_air          # 2 for indoor air temperature range 
    Nc_06_floor = (4)*n_SFH06_floor        # 2 for indoor air temperature range + 2 for temperature range of underfloor heating system
    Nc_bui     = Nc_06_air + Nc_05_air + Nc_90_air + Nc_60_air + Nc_06_floor
    
    Nc_cluster = 1                         # Electric demand at cluster level
    Nc_tot     = int(Nc_bui + Nc_cluster)  # Total number of constrains
    
    c      = [] # Coefficients of the linear objective function to be minimized (electric consumption)
    bounds = [] # Sequence of pairs (min, max) for each element in x, defining the minimum and maximum values of the decision variable
    
    """Setting the bounds and coefficients of 2006-today archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH06_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH06_air[r]))            # bounds of the decision variable (thermal power, Qth)
            c.append(float(1/COP_SFH06_air[r]))          # having x = Qth and c = 1/COP minimizes the power consumption
    
    """Setting the bounds and coefficients of 1991-2005 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH05_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH05_air[r]))
            c.append(float(1/COP_SFH05_air[r]))          # having x = Qth and c = 1/COP minimizes the electric consumption 
    
    """Setting the bounds and coefficients of 1976-1990 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH90_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH90_air[r]))
            c.append(float(1/COP_SFH90_air[r]))          # having x = Qth and c = 1/COP minimizes the electric consumption 
    
    """Setting the bounds and coefficients of 1946-1960 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH60_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH60_air[r]))
            c.append(float(1/COP_SFH60_air[r]))          # having x = Qth and c = 1/COP minimizes the electric consumption 
    
    """Setting the bounds and coefficients of 2006-today archetypes with underfloor heating system. """
    for i in range(0,n_SFH06_floor):
        for r in range(0,samples):
            bounds.append((0,Q_SFH06_floor[r])) 
            c.append(float(1/COP_SFH06_floor[r]))    # having x = Qth and c = 1/COP minimizes the electric consumption 

    I       = np.identity(samples)
    
    """Setting matrix and vector size of inequality constraints. """
    Aub_opt = np.zeros([(Nc_tot)*samples,n_Bui*samples]) # Matrix of inequality constraints
    Bub_opt = np.zeros([(Nc_tot)*samples,1])             # Vector of inequality constraints
    
    """Compilation of the matrix and vector of inequality constraints. """
    for i in range (n_SFH06_air):
        for t in range (samples):
            """2006-today archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i)*2)*samples), ((i)*samples):((i + 1)*samples)] = np.array(Z["Z_SFH06_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i)*2)*samples), 0] = float(upper["upper_SFH06_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH06_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i)*2)+1)*samples), ((i)*samples):((i + 1)*samples)] = np.array(Z["Z_SFH06_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH06_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH06_{0}".format(i + 1)][t])
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot-Nc_cluster)*samples), ((i)*samples):((i + 1)*samples)] = I[t,:]*(1/COP_SFH06_air[t])
    
    for i in range (n_SFH05_air):
        for t in range (samples):
            """1991-2005 archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air)*2)*samples), ((i + n_SFH06_air)*samples):((i + n_SFH06_air + 1)*samples)] = np.array(Z["Z_SFH05_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air)*2)*samples), 0] = float(upper["upper_SFH05_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH05_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air)*2)+1)*samples), ((i + n_SFH06_air)*samples):((i + n_SFH06_air + 1)*samples)] = np.array(Z["Z_SFH05_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH05_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH05_{0}".format(i + 1)][t])
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot-Nc_cluster)*samples), ((i + n_SFH06_air)*samples):((i + n_SFH06_air + 1)*samples)] = I[t,:]*(1/COP_SFH05_air[t])

    for i in range (n_SFH90_air):
        for t in range (samples):
            """1976-1990 archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air)*2)*samples), ((i + n_SFH06_air + n_SFH05_air)*samples):((i + n_SFH06_air + n_SFH05_air + 1)*samples)] = np.array(Z["Z_SFH90_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air)*2)*samples), 0] = float(upper["upper_SFH90_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH90_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air)*2)+1)*samples), ((i + n_SFH06_air + n_SFH05_air)*samples):((i + n_SFH06_air + n_SFH05_air + 1)*samples)] = np.array(Z["Z_SFH90_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH90_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH90_{0}".format(i + 1)][t])
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot-Nc_cluster)*samples), ((i + n_SFH06_air + n_SFH05_air)*samples):((i + n_SFH06_air + n_SFH05_air + 1)*samples)] = I[t,:]*(1/COP_SFH90_air[t])

    for i in range (n_SFH60_air):
        for t in range (samples):
            """1946-1960 archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = np.array(Z["Z_SFH60_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)*samples), 0] = float(upper["upper_SFH60_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH60_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)+1)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = np.array(Z["Z_SFH60_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH60_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH60_{0}".format(i + 1)][t])
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot-Nc_cluster)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = I[t,:]*(1/COP_SFH60_air[t])

    for i in range (n_SFH06_floor):
        for t in range (samples):
            """2006-today archetype with underfloor heating system indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air+ 1)*samples)] = np.array(Z["Z_SFH06_{0}_floor".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2)*samples), 0] = float(upper["upper_SFH06_{0}_floor".format(i + 1)][t]) - float(Mmax["Mmax_SFH06_{0}_floor".format(i + 1)][t][0])
    
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+1)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air+ 1)*samples)] = np.array(Z["Z_SFH06_{0}_floor".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH06_{0}_floor".format(i + 1)][t][0]) - float(lower["lower_SFH06_{0}_floor".format(i + 1)][t])
            
            """2006-today archetype with underfloor heating system floor temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2 + 2)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air+ 1)*samples)] = np.array(Z_floor["Z_SFH06_{0}_f".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2 + 2)*samples), 0] = float(29) - float(Mmax_floor["Mmax_SFH06_{0}_f".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+ 2 + 1)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = np.array(Z_floor["Z_SFH06_{0}_f".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+ 2 + 1)*samples), 0] = float(Mmax_floor["Mmax_SFH06_{0}_f".format(i+ 1)][t][0]) - float(15)
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot-Nc_cluster)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air + 1)*samples)] = I[t,:]*(1/COP_SFH06_floor[t])
    
    """Constraint on power consumption at cluster level. """
    for t in range (samples):
        Bub_opt[t+((Nc_bui)*samples), 0] = Eele_max[t]      # constrained electric consumption at cluster level
    
    Aub = Aub_opt.reshape((Nc_tot)*samples,n_Bui*samples)
    Bub = Bub_opt.reshape((Nc_tot)*samples,1)
    
    """Optimization result (Demand-Response scenario). """
    res = linprog(c=c, A_ub=Aub, b_ub=Bub, bounds=bounds, method='highs-ipm', options={'maxiter': 10000000})
    print("Cluster level:")
    print(res.message + "(Demand-Response scenario)")
    if res.success == False:
        print("\n              ")
        print("\nWARNING")
        print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
        print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
        print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
        print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
    print("\n------------------------------------------------------------------------")
    
    """Calculation of thermal loads (decision variable x), indoor air temperatures and building electricity consumption during the DR scenario. """

    Qres_dr   = {}  # Values of the decision variable (thermal demand in Wh)
    Tair_dr   = {}  # Indoor air temeperature (°C) - Demand-Response scenario
    Q_dr   = {}     # Thermal demand (kWh) - Demand-Rsponse scenario
    E_dr   = {}     # Electric demand (kWh) - Demand-response scenario
    
    for i in range(n_SFH06_air):
        """DR load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_dr["Qres_dr_SFH06_{0}".format(i+1)] = (res.x[(i*samples):((i + 1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH06_{0}".format(i+1)] = list(M["M_SFH06_{0}".format(i + 1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}".format(i + 1)], Qres_dr["Qres_dr_SFH06_{0}".format(i + 1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH06_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH06_{0}".format(i + 1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (DR scenario). """
        E_dr["E_dr_SFH06_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH06_{0}".format(i+1)][r]/COP_SFH06_air[r] for r in range(0, samples)]
        
    for i in range(n_SFH05_air):
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        Qres_dr["Qres_dr_SFH05_{0}".format(i+1)] = (res.x[((i + n_SFH06_air)*samples):((i + n_SFH06_air + 1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH05_{0}".format(i+1)] = list(M["M_SFH05_{0}".format( i + 1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH05_{0}".format(i + 1)], Qres_dr["Qres_dr_SFH05_{0}".format(i + 1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH05_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH05_{0}".format(i + 1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (DR scenario). """
        E_dr["E_dr_SFH05_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH05_{0}".format(i + 1)][r]/COP_SFH05_air[r] for r in range(0, samples)]
        
    for i in range(n_SFH90_air):
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        Qres_dr["Qres_dr_SFH90_{0}".format(i+1)] = (res.x[((i + n_SFH06_air + n_SFH05_air)*samples):((i + n_SFH06_air + n_SFH05_air + 1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH90_{0}".format(i+1)] = list(M["M_SFH90_{0}".format(i + 1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH90_{0}".format(i + 1)], Qres_dr["Qres_dr_SFH90_{0}".format(i + 1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH90_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH90_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (DR scenario). """
        E_dr["E_dr_SFH90_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH90_{0}".format(i + 1)][r]/COP_SFH90_air[r] for r in range(0, samples)]
    
    for i in range(n_SFH60_air):
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        Qres_dr["Qres_dr_SFH60_{0}".format(i+1)] = (res.x[((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH60_{0}".format(i+1)] = list(M["M_SFH60_{0}".format(i + 1)][3, :].reshape(samples, 1) + np.dot(Z["Z_SFH60_{0}".format(i + 1)], Qres_dr["Qres_dr_SFH60_{0}".format(i + 1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH60_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH60_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (DR scenario). """
        E_dr["E_dr_SFH60_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH60_{0}".format(i + 1)][r]/COP_SFH60_air[r] for r in range(0, samples)]

    for i in range(n_SFH06_floor):
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        Qres_dr["Qres_dr_SFH06_{0}_floor".format(i+1)] = (res.x[((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air + 1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH06_{0}_floor".format(i+1)] = list(M["M_SFH06_{0}_floor".format(i + 1)][7, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}_floor".format(i + 1)], Qres_dr["Qres_dr_SFH06_{0}_floor".format(i + 1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH06_{0}_floor".format(i+1)]       = [Qres_dr["Qres_dr_SFH06_{0}_floor".format(i + 1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (DR scenario). """
        E_dr["E_dr_SFH06_{0}_floor".format(i+1)]       = [Q_dr["Q_dr_SFH06_{0}_floor".format(i + 1)][r]/COP_SFH06_floor[r] for r in range(0, samples)]
        
    """Conversion of dictionary to data frame for calculation of aggregate electric consumption. """
    E_dr = pd.DataFrame.from_dict(E_dr)
    
    """Calculating electric consumption at the aggregate level during DR scenario. """
    Etot_dr = E_dr.sum(axis=1)
    Edr_sum = float(Etot_dr.sum())
    
    Etot_dr = [Etot_dr[i] for i in range(0,samples)]  
    
    import matplotlib.pyplot as plt
    time = [i*timestep for i in range(0, samples)]  # time array
    
    """Do plot. """
    if do_plot_E_dr_cluster:
        """If True, plot cluster electric consumption trends. """
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Time (hours)', fontsize=14)
        ax1.set_ylabel('$\.P_\mathrm{cluster}$ (kW$_\mathrm{el}$)', fontsize=14)
        for i in range (n_days):
            ax1.axvspan((time_max[i]*timestep), (time_max[i]*timestep) + duration_dr, facecolor='r', alpha=0.10)
        ax1.axvspan(24*(1), 24*(n_days-1), facecolor='k', alpha=0.04)
        ax1.plot(time, Etot_bl, 'k', label="BL_scenario", linewidth=1.0)
        ax1.plot(time, Etot_dr, 'r', label="DR_scenario", linewidth=1.0)
        
        ax1.set_xlim(0, duration-1)
        ax1.grid(alpha=0.3)
        plt.legend(frameon=False, fontsize=12)
        fig.set_size_inches(w=10, h=4)
        fig.tight_layout()
        plt.savefig("Eele_cluster", dpi=200)
        
    plt.show()
    plt.close('all')    
    
    if do_plot_E_dr_single:
        """If True, plot single buildings electric consumption trends. """
        for i in range (n_SFH06_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            for n in range (n_days):
                ax1.axvspan((time_max[n]*timestep), (time_max[n]*timestep) + duration_dr, facecolor='r', alpha=0.10)
            ax1.plot(time,E_dr["E_dr_SFH06_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH06_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH06_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH06_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=3)
            fig.tight_layout()
            
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH05_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            for n in range (n_days):
                ax1.axvspan((time_max[n]*timestep), (time_max[n]*timestep) + duration_dr, facecolor='r', alpha=0.10)
            ax1.plot(time,E_dr["E_dr_SFH05_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH05_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH05_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH05_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH90_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            for n in range (n_days):
                ax1.axvspan((time_max[n]*timestep), (time_max[n]*timestep) + duration_dr, facecolor='r', alpha=0.10)
            ax1.plot(time,E_dr["E_dr_SFH90_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH90_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH90_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH90_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH60_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            for n in range (n_days):
                ax1.axvspan((time_max[n]*timestep), (time_max[n]*timestep) + duration_dr, facecolor='r', alpha=0.10)
            ax1.plot(time,E_dr["E_dr_SFH60_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH60_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH60_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH60_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
    
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH06_floor):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            for n in range (n_days):
                ax1.axvspan((time_max[n]*timestep), (time_max[n]*timestep) + duration_dr, facecolor='r', alpha=0.10)
            ax1.plot(time,E_dr["E_dr_SFH06_{0}_floor".format(i+1)],'r', linewidth=1, label = "E_dr_SFH06_{0}_floor".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH06_{0}_floor".format(i+1)],'k', linewidth=1, label = "E_bl_SFH06_{0}_floor".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
            
        plt.show()
        plt.close('all')

    if do_plot_Tair_dr_single:
        """If True, plot single buildings indoor air temperature trends. """
        for i in range (n_SFH06_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            for n in range (n_days):
                ax1.axvspan((time_max[n]*timestep), (time_max[n]*timestep) + duration_dr, facecolor='r', alpha=0.10)
            ax1.plot(time,Tair_dr["Tair_dr_SFH06_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH06_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH06_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH06_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
            
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH05_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            for n in range (n_days):
                ax1.axvspan((time_max[n]*timestep), (time_max[n]*timestep) + duration_dr, facecolor='r', alpha=0.10)
            ax1.plot(time,Tair_dr["Tair_dr_SFH05_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH05_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH05_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH05_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH90_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            for n in range (n_days):
                ax1.axvspan((time_max[n]*timestep), (time_max[n]*timestep) + duration_dr, facecolor='r', alpha=0.10)
            ax1.plot(time,Tair_dr["Tair_dr_SFH90_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH90_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH90_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH90_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH60_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            for n in range (n_days):
                ax1.axvspan((time_max[n]*timestep), (time_max[n]*timestep) + duration_dr, facecolor='r', alpha=0.10)
            ax1.plot(time,Tair_dr["Tair_dr_SFH60_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH60_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH60_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH60_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
    
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH06_floor):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            for n in range (n_days):
                ax1.axvspan((time_max[n]*timestep), (time_max[n]*timestep) + duration_dr, facecolor='r', alpha=0.10)
            ax1.plot(time,Tair_dr["Tair_dr_SFH06_{0}_floor".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH06_{0}_floor".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH06_{0}_floor".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH06_{0}_floor".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
    """Data export to excel. """
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writerEele = pd.ExcelWriter(name_file + '.xlsx', engine='xlsxwriter')

    # Write each dataframe to a different worksheet.
    E_bl.to_excel(writerEele, sheet_name='Eele_BL')
    E_dr.to_excel(writerEele, sheet_name='Eele_DR')

    # Close the Pandas Excel writer and output the Excel file.
    writerEele.close()
    
    """Electricity consumption during BL and DR scenarios. """
    print("\nCluster electricity consumption of about " +
          str(round(Ebl_sum,2)) + " kWh during baseline (BL) scenario")
    print("\n------------------------------------------------------------------------")
    
    Diff   = ((Edr_sum - Ebl_sum)/Ebl_sum)*100
    
    if Diff >= 0:
        print("\nCluster electricity consumption of about " +
              str(round(Edr_sum,2)) + " kWh during DR scenario with " + str(round(f_red*100,2)) + "% peak load reduction,"
              "\nwith a " + str(round(Diff,2)) + "% increase compared to BL scenario")
        print("\n------------------------------------------------------------------------")
    else:
        print("\nCluster electricity consumption of about " +
              str(round(Edr_sum,2)) + " kWh during DR scenario with " + str(round(f_red*100,2)) + "% peak load reduction,"
              "\nwith a " + str(round(Diff*(-1),2)) + "% reduction compared to BL scenario")
        print("\n------------------------------------------------------------------------")
        
    return(Etot_dr)

def pv_centralized(Locality, day, month, duration, timestep, Thermal_load, 
                            n_SFH06_air, n_SFH05_air, n_SFH90_air, n_SFH60_air,n_SFH06_floor, 
                            HP_data_type, 
                            Supply_air, rated_cap_06_air, rated_cap_05_air, rated_cap_90_air, rated_cap_60_air, rated_COP_06_air, rated_COP_05_air, rated_COP_90_air, rated_COP_60_air,
                            Supply_floor, rated_cap_06_floor, rated_COP_06_floor,
                            Setpoint_type, setpoint, th_tolerance_BL_sup_air, th_tolerance_BL_min_air, th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor,
                            th_tolerance_DR_sup_air, th_tolerance_DR_min_air, th_tolerance_DR_sup_floor, th_tolerance_DR_min_floor,
                            pv_size, panel_dim, pv_rp_06_air, pv_rp_05_air, pv_rp_90_air, pv_rp_60_air, pv_rp_06_floor, f_limit,
                            do_plot_E_dr_single, do_plot_E_dr_cluster,do_plot_Tair_dr_single, name_file):
    """pv_centralized is a function based on linear programming that minimizes the aggregate electricity consumption 
            of a cluster of buildings drawn from the grid during the application of a strategy of shifting loads to centralized 
            PV generation hours (i.e., shared energy resources). """
            
    # number of buildings:
    n_SFH06_air = int(n_SFH06_air)
    n_SFH05_air = int(n_SFH05_air)
    n_SFH90_air = int(n_SFH90_air)
    n_SFH60_air = int(n_SFH60_air)
    
    if Thermal_load == "heating":
        n_SFH06_floor = int(n_SFH06_floor) # Radiant floor emission system for the SFH06 archetype (age class 2006-today) only during heating
    else: 
        n_SFH06_floor = 0
    
    n_Bui       = n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air + n_SFH06_floor # Total number of buildings composing the cluster
    
    samples = int(duration/timestep) # Calculation of total steps

    T_est = MeteoFile(Locality, day, month, duration, timestep)[0] # Outdoor air temperature (°C)

    """Call of the performance characteristics of heat pumps based on the definition of heat loads 
           and attribute definition for calling temperature setpoint functions"""
    if Thermal_load == 'cooling':
        """Call performance characteristics and capacity of heat pumps. """
        if n_SFH06_air > 0:
            COP_SFH06_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH06_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH05_air > 0:
            COP_SFH05_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH05_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH90_air > 0:
            COP_SFH90_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH90_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH60_air > 0:
            COP_SFH60_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH60_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[1] for i in range(0,samples)]
            
        """Thermostat mode setting in space cooling"""
        thermostat = 'th_cooling'
    elif Thermal_load == 'heating':
        if n_SFH06_air > 0:
            COP_SFH06_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH06_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH05_air > 0:
            COP_SFH05_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH05_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH90_air > 0:
            COP_SFH90_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH90_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH60_air > 0:
            COP_SFH60_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH60_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH06_floor > 0:
            COP_SFH06_floor = [HPperformance.HP_heating(T_est[i],Supply_floor,rated_cap_06_air,rated_COP_06_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH06_floor = [HPperformance.HP_heating(T_est[i],Supply_floor,rated_cap_06_air,rated_COP_06_air,HP_data_type)[1] for i in range(0,samples)]
        
        """Thermostat mode setting in space heating"""
        thermostat = 'th_heating'
    else:
        print('Thermal load not available')
    
    """Definition of occupancy profiles, relative internal gains (people, equipment and lighting) and setpoint temperature for each building 
           (useful for diversification of consumption loads within the cluster)"""
    occ   = {} # occupancy
    gains = {} # internal gains (W)
    tsp   = {} # thermostat set-points (°C)
    
    if Setpoint_type == "pre_defined":
        """Pre-defined temperature set-point definition
               For adding profiles or modifying existing ones, go to the User_pattern.py file and see the functions 
               th_heating (set-points for space heating) or th_cooling (set-points for space cooling)"""
        for i in range(n_SFH06_air):
            tsp["tspSFH06_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i]
        for i in range(n_SFH05_air):
            tsp["tspSFH05_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air]
        for i in range(n_SFH90_air):
            tsp["tspSFH90_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air]
        for i in range(n_SFH60_air):
            tsp["tspSFH60_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air]
        for i in range(n_SFH06_floor):
            tsp["tspSFH06_{0}_floor".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air]
    
    elif Setpoint_type == "user_defined":
        """User-defined temperature set-point definition"""
        t_setpoint = []
        for i in range (samples):
            t_setpoint.append(setpoint)
        for i in range(n_SFH06_air):
            tsp["tspSFH06_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH05_air):
            tsp["tspSFH05_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH90_air):
            tsp["tspSFH90_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH60_air):
            tsp["tspSFH60_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH06_floor):
            tsp["tspSFH06_{0}_floor".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
    
    """Occupancy patterns and related internal gain profiles definition"""
    for i in range(n_SFH06_air):
        occ["occupancySFH06_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i]
        gains["gainsSFH06_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i]
    for i in range(n_SFH05_air):
        occ["occupancySFH05_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air]
        gains["gainsSFH05_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air]
    for i in range(n_SFH90_air):
        occ["occupancySFH90_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air]
        gains["gainsSFH90_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air]
    for i in range(n_SFH60_air):
        occ["occupancySFH60_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air]
        gains["gainsSFH60_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air]
    for i in range(n_SFH06_floor):
        occ["occupancySFH06_{0}_floor".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air]
        gains["gainsSFH06_{0}_floor".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air]

    """In the following, thermal loads provided as input to the state-space model formulation are calculated via linear programming . 
           Specifically, the thermal load of each building is minimized to maintain the indoor temperature set-point.
           Linear programming solves problems of the following form: 
               minimize c @ x
               such that:
               A_ub @ x <= b_ub
               A_eq @ x == b_eq
               lb <= x <= ub
               
               For more information see the documentation of the python library scipy.optimize.linprog"""
           
    c      = {}     # Coefficients of the linear objective function to be minimized
    Aub    = {}     # Matrix of inequality constraints
    Bub    = {}     # Vector of inequality constraints
    bounds = {}     # Sequence of pairs (min, max) for each element in x, defining the minimum and maximum values of the decision variable
    
    Z      = {}     # Variable describing the state of the system containing the decision variable (e.g., thermal demand)
    M      = {}     # Variable describing the state of the indoor air temperture
    Mmax      = {}  # Variable describing the state of the indoor air temperture
    res_opt   = {}  # Optimization result - baseline scenario
    Qres_bl   = {}  # Values of the decision variable (thermal demand in Wh)
    Tair_bl   = {}  # Indoor air temeperature (°C) - baseline scenario
    Q_bl   = {}     # Thermal demand (kWh) - baseline scenario
    E_bl   = {}     # Electric demand (kWh) - baseline scenario
    upper   = {}    # Upper range thermostat set-point profile (°C) - baseline scenario
    lower   = {}    # Lower range thermostat set-point profile (°C) - baseline scenario
    
    Mmax_floor = {} # Variable describing the state of the floor temperture
    Z_floor    = {} # Variable describing the state of the floor system containing the decision variable (e.g., thermal demand)

    """Calling archetype (SFHs) functions from the "Archetypes.py" module. """
    for i in range(n_SFH06_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH06 buildings (age class 2006-today) with air emission system."""
        c["c_SFH06_{0}".format(i+1)]             = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[0]
        Aub["Aub_SFH06_{0}".format(i+1)]         = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[1]
        Bub["Bub_SFH06_{0}".format(i+1)]         = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[2]
        bounds["bounds_SFH06_{0}".format(i+1)]   = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH06_{0}".format(i+1)]             = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[4]
        Z["Z_SFH06_{0}".format(i+1)]             = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[5]
        Mmax["Mmax_SFH06_{0}".format(i+1)]       = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH06_{0}".format(i+1)]     = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[7]
        lower["lower_SFH06_{0}".format(i+1)]     = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH06_{0}".format(i+1)] = linprog(c=c["c_SFH06_{0}".format(i+1)], A_ub=Aub["Aub_SFH06_{0}".format(i+1)], b_ub=Bub["Bub_SFH06_{0}".format(i+1)],bounds=bounds["bounds_SFH06_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 2006-today building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH06_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH06_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH06_{0}".format(i+1)] = (res_opt["res_opt_SFH06_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH06_{0}".format(i+1)] = list(M["M_SFH06_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}".format(i+1)], Qres_bl["Qres_bl_SFH06_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH06_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH06_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH06_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH06_{0}".format(i+1)][r]/COP_SFH06_air[r] for r in range(0, samples)]
    for i in range(n_SFH05_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH05 buildings (age class 1991-2005) with air emission system."""
        c["c_SFH05_{0}".format(i+1)]             = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[0]
        Aub["Aub_SFH05_{0}".format(i+1)]         = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[1]
        Bub["Bub_SFH05_{0}".format(i+1)]         = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[2]
        bounds["bounds_SFH05_{0}".format(i+1)]   = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH05_{0}".format(i+1)]             = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[4]
        Z["Z_SFH05_{0}".format(i+1)]             = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[5]
        Mmax["Mmax_SFH05_{0}".format(i+1)]       = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH05_{0}".format(i+1)] = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[7]
        lower["lower_SFH05_{0}".format(i+1)] = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH05_{0}".format(i+1)] = linprog(c=c["c_SFH05_{0}".format(i+1)], A_ub=Aub["Aub_SFH05_{0}".format(i+1)], b_ub=Bub["Bub_SFH05_{0}".format(i+1)],bounds=bounds["bounds_SFH05_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 1990-2005 building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH05_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH05_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH05_{0}".format(i+1)] = (res_opt["res_opt_SFH05_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH05_{0}".format(i+1)] = list(M["M_SFH05_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH05_{0}".format(i+1)], Qres_bl["Qres_bl_SFH05_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH05_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH05_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH05_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH05_{0}".format(i+1)][r]/COP_SFH05_air[r] for r in range(0, samples)]
    for i in range(n_SFH90_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH90 buildings (age class 1976-1990) with air emission system."""
        c["c_SFH90_{0}".format(i+1)]             = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[0]
        Aub["Aub_SFH90_{0}".format(i+1)]         = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[1]
        Bub["Bub_SFH90_{0}".format(i+1)]         = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[2]
        bounds["bounds_SFH90_{0}".format(i+1)]   = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH90_{0}".format(i+1)]             = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[4]
        Z["Z_SFH90_{0}".format(i+1)]             = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[5]
        Mmax["Mmax_SFH90_{0}".format(i+1)]       = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH90_{0}".format(i+1)]     = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[7]
        lower["lower_SFH90_{0}".format(i+1)]     = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[8]
        
        """Optimization result (baseline scenario)"""
        res_opt["res_opt_SFH90_{0}".format(i+1)] = linprog(c=c["c_SFH90_{0}".format(i+1)], A_ub=Aub["Aub_SFH90_{0}".format(i+1)], b_ub=Bub["Bub_SFH90_{0}".format(i+1)],bounds=bounds["bounds_SFH90_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 1976-1990 building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH90_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH90_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH90_{0}".format(i+1)] = (res_opt["res_opt_SFH90_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH90_{0}".format(i+1)] = list(M["M_SFH90_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH90_{0}".format(i+1)], Qres_bl["Qres_bl_SFH90_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH90_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH90_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH90_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH90_{0}".format(i+1)][r]/COP_SFH90_air[r] for r in range(0, samples)]
    for i in range(n_SFH60_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH60 buildings (age class 1946-1960) with air emission system."""
        c["c_SFH60_{0}".format(i+1)]             = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[0]
        Aub["Aub_SFH60_{0}".format(i+1)]         = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[1]
        Bub["Bub_SFH60_{0}".format(i+1)]         = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[2]
        bounds["bounds_SFH60_{0}".format(i+1)]   = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH60_{0}".format(i+1)]             = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[4]
        Z["Z_SFH60_{0}".format(i+1)]             = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[5]
        Mmax["Mmax_SFH60_{0}".format(i+1)]       = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH60_{0}".format(i+1)] = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[7]
        lower["lower_SFH60_{0}".format(i+1)] = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH60_{0}".format(i+1)] = linprog(c=c["c_SFH60_{0}".format(i+1)], A_ub=Aub["Aub_SFH60_{0}".format(i+1)], b_ub=Bub["Bub_SFH60_{0}".format(i+1)],bounds=bounds["bounds_SFH60_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 1946-1960 building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH60_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH60_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")

        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH60_{0}".format(i+1)] = (res_opt["res_opt_SFH60_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH60_{0}".format(i+1)] = list(M["M_SFH60_{0}".format(i+1)][3, :].reshape(samples, 1) + np.dot(Z["Z_SFH60_{0}".format(i+1)], Qres_bl["Qres_bl_SFH60_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH60_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH60_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH60_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH60_{0}".format(i+1)][r]/COP_SFH60_air[r] for r in range(0, samples)]
    for i in range(n_SFH06_floor):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH06 buildings (age class 2006-today) with underfloor heating system."""
        c["c_SFH06_{0}_floor".format(i+1)]             = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[0]
        Aub["Aub_SFH06_{0}_floor".format(i+1)]         = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[1]
        Bub["Bub_SFH06_{0}_floor".format(i+1)]         = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[2]
        bounds["bounds_SFH06_{0}_floor".format(i+1)]   = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[3]
        
        """Calling of variables referring to the state of indoor air temperature node. """
        M["M_SFH06_{0}_floor".format(i+1)]             = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[4]
        Z["Z_SFH06_{0}_floor".format(i+1)]             = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[5]
        Mmax["Mmax_SFH06_{0}_floor".format(i+1)]       = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[6]
        
        """Calling of variables referring to the state of floor node. """
        Z_floor["Z_SFH06_{0}_f".format(i+1)]           = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[9]
        Mmax_floor["Mmax_SFH06_{0}_f".format(i+1)]     = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[10]
        
        """Calling of variables referring to the state of indoor air temperature. """
        upper["upper_SFH06_{0}_floor".format(i+1)]     = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[7]
        lower["lower_SFH06_{0}_floor".format(i+1)]     = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH06_{0}_floor".format(i+1)] = linprog(c=c["c_SFH06_{0}_floor".format(i+1)], A_ub=Aub["Aub_SFH06_{0}_floor".format(i+1)], b_ub=Bub["Bub_SFH06_{0}_floor".format(i+1)],bounds=bounds["bounds_SFH06_{0}_floor".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype building {0}".format(i+1) + " with underfloor heating system:")
        print(res_opt["res_opt_SFH06_{0}_floor".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH06_{0}_floor".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH06_{0}_floor".format(i+1)] = (res_opt["res_opt_SFH06_{0}_floor".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH06_{0}_floor".format(i+1)] = list(M["M_SFH06_{0}_floor".format(i+1)][7, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}_floor".format(i+1)], Qres_bl["Qres_bl_SFH06_{0}_floor".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH06_{0}_floor".format(i+1)]       = [Qres_bl["Qres_bl_SFH06_{0}_floor".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH06_{0}_floor".format(i+1)]       = [Q_bl["Q_bl_SFH06_{0}_floor".format(i+1)][r]/COP_SFH06_floor[r] for r in range(0, samples)]
    
    """Conversion of dictionary to data frame for calculation of aggregate electric power consumption. """
    E_bl = pd.DataFrame.from_dict(E_bl)
    
    """Calculating base electrical loads at the aggregate level. """
    Etot_bl = E_bl.sum(axis=1)
    Ebl_sum = float(Etot_bl.sum())
    Etot_bl = [Etot_bl[i] for i in range(0,samples)]  
    
    """Definition of the centralized photovoltaic generation profile by application of the python pvlib library. 
            Specifically, this centralized generation scenario represents solar PV generation available to the entire cluster of buildings. 
            Consequently, during centralized generation, energy resources are shared within the building cluster. """
    
    start      = MeteoFile(Locality,int(day),int(month),duration,timestep)[10] # start time to give as input to the PV_load function
   
    if pv_size == "calculated":
        """If "calculated", the rated power of the PV systems available to archetypes are calculated according to the Italian regulation Dlgs. 28/2011. """
        panel = float(panel_dim)                                  # W
     
        PV_06_air    = (Archetypes.floor_area_SFH2006/50)*1000    # W
        PV_05_air    = (Archetypes.floor_area_SFH2005/50)*1000    # W
        PV_90_air    = (Archetypes.floor_area_SFH1990/50)*1000    # W
        PV_60_air    = (Archetypes.floor_area_SFH1960/50)*1000    # W
        PV_06_floor  = (Archetypes.floor_area_SFH2006/50)*1000    # W
        
        n_panels_06_air   = math.ceil(PV_06_air/panel)            # Number of panels defining the photovoltaic system of the 2006-today archetype with air system
        n_panels_05_air   = math.ceil(PV_05_air/panel)            # Number of panels defining the photovoltaic system of the 1991-2005 archetype with air system
        n_panels_90_air   = math.ceil(PV_90_air/panel)            # Number of panels defining the photovoltaic system of the 1976-1990 archetype with air system
        n_panels_60_air   = math.ceil(PV_60_air/panel)            # Number of panels defining the photovoltaic system of the 1946-1960 archetype with air system
        n_panels_06_floor = math.ceil(PV_06_floor/panel)          # Number of panels defining the photovoltaic system of the 2006-today archetype with underfloor heating system
        
        PV_rp_06_air   = n_panels_06_air*panel/1000               # PV array's rated power (kW)
        PV_rp_05_air   = n_panels_05_air*panel/1000               # PV array's rated power (kW)
        PV_rp_90_air   = n_panels_90_air*panel/1000               # PV array's rated power (kW)
        PV_rp_60_air   = n_panels_60_air*panel/1000               # PV array's rated power (kW)
        PV_rp_06_floor = n_panels_06_floor*panel/1000             # PV array's rated power (kW)
    elif pv_size == "user_defined":
        """If "user-defined", the rated power of the PV systems are set by the user. """
        PV_rp_06_air   = pv_rp_06_air   # PV array's rated power (kW)
        PV_rp_05_air   = pv_rp_05_air   # PV array's rated power (kW)
        PV_rp_90_air   = pv_rp_90_air   # PV array's rated power (kW)
        PV_rp_60_air   = pv_rp_60_air   # PV array's rated power (kW)
        PV_rp_06_floor = pv_rp_06_floor # PV array's rated power (kW)
    
    """Rated power of the centralized photovoltaic system. """
    rated_power = PV_rp_06_air*n_SFH06_air + PV_rp_05_air*n_SFH05_air + PV_rp_90_air*n_SFH90_air + PV_rp_60_air*n_SFH60_air + PV_rp_06_floor*n_SFH06_floor
    
    """Calculation of centralized photovoltaic generation profile. """
    P_ele_pv = PV_load(rated_power, start, duration, timestep, Locality) # (W)
    
    """------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    The Demand-Response (DR) event is modeled below. 
    From the definition of the PV generation profile and electricity demand during the baseline scenario, electricity demand 
    during the DR scenario is managed by shifting aggregate loads to centralized generation hours.
    To this end, an additional optimization constraint is introduced to minimize the energy demand of the entire cluster. 
    Specifically, aggregated electricity from the grid is minimized:
        min (Q_dr/COP - Upv * P_ele_pv), 
        with
        Q_dr = thermal load (variable 1),
        Upv  = PV utilization factor (variable 2).
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
    """For each simulation day, calculation of the aggregate peak load. """
    n_days        = int(duration/24) # number of simulation days
    
    Max         = []
    Max_dr      = []

    for n in range (n_days):
        Max.append(np.amax(np.array(Etot_bl[int(samples/n_days)*n:int(samples/n_days)*(n + 1)])))
        for i in range (int(samples/n_days)*n,int(samples/n_days)*(n + 1)):
            Max_dr.append(Max[n]) # maximum electricity consumption during the baseline scenario

    """Definition of the electricity demand of the cluster of buildings constrained to remain below a certain limit to reduce rebound effects.. """
    Eele_max = []

    for t in range(samples):
            Eele_max.append(float(Max_dr[t]*1000)*(float(f_limit)))
    
    """Definition of temperature profiles during the DR event to unlock building energy flexibility. """
    if Thermal_load == 'cooling':
        for i in range(n_SFH06_air):
            for t in range(samples):
                if P_ele_pv[t] > 0:
                    lower["lower_SFH06_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to increase PV consumption
                else:
                    upper["upper_SFH06_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for cooling to reduce electric consumption from the grid
        for i in range(n_SFH05_air):
            for t in range(samples):
                if P_ele_pv[t] > 0:
                    lower["lower_SFH05_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to increase PV consumption
                else:
                    upper["upper_SFH05_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for cooling to reduce electric consumption from the grid
        for i in range(n_SFH90_air):
            for t in range(samples):
                if P_ele_pv[t] > 0:
                    lower["lower_SFH90_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to increase PV consumption
                else:
                    upper["upper_SFH90_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for cooling to reduce electric consumption from the grid
        for i in range(n_SFH60_air):
            for t in range(samples):
                if P_ele_pv[t] > 0:
                    lower["lower_SFH60_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to increase PV consumption
                else:
                    upper["upper_SFH60_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for cooling to reduce electric consumption from the grid

    elif Thermal_load == 'heating':
        for i in range(n_SFH06_air):
            for t in range(samples):
                if P_ele_pv[t] > 0:
                    upper["upper_SFH06_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for heating to increase PV consumption
                else:
                    lower["lower_SFH06_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to reduce electric consumption from the grid
        for i in range(n_SFH05_air):
            for t in range(samples):
                if P_ele_pv[t] > 0:
                    upper["upper_SFH05_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for heating to increase PV consumption
                else:
                    lower["lower_SFH05_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to reduce electric consumption from the grid
        for i in range(n_SFH90_air):
            for t in range(samples):
                if P_ele_pv[t] > 0:
                    upper["upper_SFH90_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for heating to increase PV consumption
                else:
                    lower["lower_SFH90_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to reduce electric consumption from the grid
        for i in range(n_SFH60_air):
            for t in range(samples):
                if P_ele_pv[t] > 0:
                    upper["upper_SFH60_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for heating to increase PV consumption
                else:
                    lower["lower_SFH60_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to reduce electric consumption from the grid
        for i in range(n_SFH06_floor):
            for t in range(samples):
                if P_ele_pv[t] > 0:
                    upper["upper_SFH06_{0}_floor".format(i+1)][t] += float(th_tolerance_DR_sup_floor) # Increase in indoor air temperature set point for heating to increase PV consumption
                else:
                    lower["lower_SFH06_{0}_floor".format(i+1)][t] -= float(th_tolerance_DR_min_floor) # Reduction of indoor air temperature set point for cooling to reduce electric consumption from the grid
    
    """In the following, thermal loads provided as input to the state-space model formulation are calculated via linear programming. 
           Specifically, the thermal load of each building is minimized to maintain the indoor temperature set-point, 
           and the electricity demand of the cluster of buildings is constrained to remain below a certain 
           percentage reduction (Eele_max) for the duration of the peak-shaving event (duration_dr).
           Linear programming solves problems of the following form: 
               minimize c @ x
               such that:
               A_ub @ x <= b_ub
               A_eq @ x == b_eq
               lb <= x <= ub
               
               For more information see the documentation of the python library scipy.optimize.linprog"""

    """Number of decision variables.
            variable 1 = thermal load
            variable 2 = PV utilization factor"""
    Nv = 2 # number of decision variables
    
    """Definition of the number of constraints at the individual building and cluster level. """
    Nc_SFH06_air_SS  = (2)*n_SFH06_air          # 2 for indoor air temperature range 
    Nc_SFH06_air_PV  = (1)*n_SFH06_air          # 1 for PV consumption

    Nc_SFH05_air_SS  = (2)*n_SFH05_air          # 2 for indoor air temperature range 
    Nc_SFH05_air_PV  = (1)*n_SFH05_air          # 1 for PV consumption

    Nc_SFH90_air_SS  = (2)*n_SFH90_air          # 2 for indoor air temperature range 
    Nc_SFH90_air_PV  = (1)*n_SFH90_air          # 1 for PV consumption

    Nc_SFH60_air_SS  = (2)*n_SFH60_air          # 2 for indoor air temperature range 
    Nc_SFH60_air_PV  = (1)*n_SFH60_air          # 1 for PV consumption

    Nc_06_floor_SS   = (4)*n_SFH06_floor        # 2 for indoor air temperature range + 2 for temperature range of underfloor heating system
    Nc_06_floor_PV   = (1)*n_SFH06_floor        # 1 for PV consumption

    Nc_bui_SS     = Nc_SFH06_air_SS + Nc_SFH05_air_SS + Nc_SFH90_air_SS + Nc_SFH60_air_SS + Nc_06_floor_SS 
    Nc_bui_PV     = Nc_SFH06_air_PV + Nc_SFH05_air_PV + Nc_SFH90_air_PV + Nc_SFH60_air_PV + Nc_06_floor_PV 
    
    Nc_cluster = 2 # PV consumption, e.g., sum(Upv) <= 1, and max electricity consumption at cluster level
    Nc_tot     = int(Nc_bui_SS + Nc_bui_PV + Nc_cluster)  # Total number of constrains
    
    
    bounds = []                         # Sequence of pairs (min, max) for each element in x, defining the minimum and maximum values of the decision variable
    c = np.zeros([1,n_Bui*Nv*samples])  # Coefficients of the linear objective function to be minimized (electric consumption)
    
    """VARIABLE 1: thermal loads. """
    """Setting the bounds and coefficients of 2006-today archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH06_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH06_air[r]))              # bounds of the decision variable (thermal power, Qth)
            c[0,r+(i*samples)] = float(1/COP_SFH06_air[r]) # having x1 = Qth and c = 1/COP minimizes the power consumption
    
    """Setting the bounds and coefficients of 1991-2005 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH05_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH05_air[r]))                          # bounds of the decision variable (thermal power, Qth)
            c[0,r+(i+n_SFH06_air)*samples] = float(1/COP_SFH05_air[r]) # having x1 = Qth and c = 1/COP minimizes the power consumption
    
    """Setting the bounds and coefficients of 1976-1990 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH90_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH90_air[r]))                                      # bounds of the decision variable (thermal power, Qth)
            c[0,r+(i+n_SFH06_air+n_SFH05_air)*samples] = float(1/COP_SFH90_air[r]) # having x1 = Qth and c = 1/COP minimizes the power consumption
    
    """Setting the bounds and coefficients of 1946-1960 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH60_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH60_air[r]))                                                  # bounds of the decision variable (thermal power, Qth)
            c[0,r+(i+n_SFH06_air+n_SFH05_air+n_SFH90_air)*samples] = float(1/COP_SFH60_air[r]) # having x1 = Qth and c = 1/COP minimizes the power consumption
    
    """Setting the bounds and coefficients of 2006-today archetypes with underfloor heating system. """
    for i in range(0,n_SFH06_floor):
        for r in range(0,samples):
            bounds.append((0,Q_SFH06_floor[r]))                                                              # bounds of the decision variable (thermal power, Qth)
            c[0,r+(i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air)*samples] = float(1/COP_SFH06_floor[r]) # having x1 = Qth and c = 1/COP minimizes the power consumption
    
    """VARIABLE 2: PV utilization factor. """
    """Setting the bounds and coefficients of 2006-today archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH06_air):
        for r in range(0,samples):
            bounds.append((0,1))                                    # bounds of PV utilization factor
            c[0,r+((i+n_Bui)*samples)] = float((-1)*P_ele_pv[r]) # having x2 = Upv and c = (-1000)*P_ele_pv maximizes PV consumption
    
    """Setting the bounds and coefficients of 1991-2005 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH05_air):
        for r in range(0,samples):
            bounds.append((0,1))                                              # bounds of PV utilization factor
            c[0,r+(i+n_SFH06_air+n_Bui)*samples] = float((-1)*P_ele_pv[r]) # having x2 = Upv and c = (-1000)*P_ele_pv maximizes PV consumption
    
    """Setting the bounds and coefficients of 1976-1990 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH90_air):
        for r in range(0,samples):
            bounds.append((0,1))                                                          # bounds of PV utilization factor
            c[0,r+(i+n_SFH06_air+n_SFH05_air+n_Bui)*samples] = float((-1)*P_ele_pv[r]) # having x2 = Upv and c = (-1000)*P_ele_pv maximizes PV consumption
    
    """Setting the bounds and coefficients of 1946-1960 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH60_air):
        for r in range(0,samples):
            bounds.append((0,1))                                                                      # bounds of PV utilization factor
            c[0,r+(i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_Bui)*samples] = float((-1)*P_ele_pv[r]) # having x2 = Upv and c = (-1000)*P_ele_pv maximizes PV consumption
    
    """Setting the bounds and coefficients of 2006-today archetypes with underfloor heating system. """
    for i in range(0,n_SFH06_floor):
        for r in range(0,samples):
            bounds.append((0,1))                                                                                  # bounds of PV utilization factor
            c[0,r+(i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air+n_Bui)*samples] = float((-1)*P_ele_pv[r]) # having x2 = Upv and c = (-1000)*P_ele_pv maximizes PV consumption

    I         = np.identity(samples)
    
    """Setting matrix and vector size of inequality constraints. """
    Aub_opt = np.zeros([(Nc_tot)*samples,n_Bui*samples*Nv])
    Bub_opt = np.empty([(Nc_tot)*samples,1])

    """Compilation of the matrix and vector of inequality constraints. """
    for i in range (n_SFH06_air):
        for t in range (samples):
            """2006-today archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i)*2)*samples), ((i)*samples):((i + 1)*samples)] = np.array(Z["Z_SFH06_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i)*2)*samples), 0] = float(upper["upper_SFH06_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH06_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i)*2)+1)*samples), ((i)*samples):((i + 1)*samples)] = np.array(Z["Z_SFH06_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH06_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH06_{0}".format(i + 1)][t])
            
            """Photovoltaic consumption constraints. """
            Aub_opt[t+((i + Nc_bui_SS)*samples), ((i)*samples):((i + 1)*samples)] = I[t,:]*((-1)*(1/COP_SFH06_air[t]))
            Aub_opt[t+((i + Nc_bui_SS)*samples), ((i + n_Bui)*samples):((i + n_Bui + 1)*samples)] = I[t,:]*(P_ele_pv[t])
            Bub_opt[t+((i + Nc_bui_SS)*samples), 0] = 0
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot - 1)*samples), ((i)*samples):((i + 1)*samples)] = I[t,:]*(1/COP_SFH06_air[t])
    
    for i in range (n_SFH05_air):
        for t in range (samples):
            """1991-2005 archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air)*2)*samples), ((i + n_SFH06_air)*samples):((i + n_SFH06_air + 1)*samples)] = np.array(Z["Z_SFH05_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air)*2)*samples), 0] = float(upper["upper_SFH05_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH05_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air)*2)+1)*samples), ((i + n_SFH06_air)*samples):((i + n_SFH06_air + 1)*samples)] = np.array(Z["Z_SFH05_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH05_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH05_{0}".format(i + 1)][t])
            
            """Photovoltaic consumption constraints. """
            Aub_opt[t+((i + n_SFH06_air + Nc_bui_SS)*samples), ((i + n_SFH06_air)*samples):((i + n_SFH06_air + 1)*samples)] = I[t,:]*((-1)*(1/COP_SFH05_air[t]))
            Aub_opt[t+((i + n_SFH06_air + Nc_bui_SS)*samples), ((i + n_SFH06_air + n_Bui)*samples):((i + n_SFH06_air + n_Bui + 1)*samples)] = I[t,:]*(P_ele_pv[t])
            Bub_opt[t+((i + n_SFH06_air + Nc_bui_SS)*samples), 0] = 0
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot - 1)*samples), ((i + n_SFH06_air)*samples):((i + n_SFH06_air + 1)*samples)] = I[t,:]*(1/COP_SFH05_air[t])

    for i in range (n_SFH90_air):
        for t in range (samples):
            """1976-1990 archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air)*2)*samples), ((i + n_SFH06_air + n_SFH05_air)*samples):((i + n_SFH06_air + n_SFH05_air + 1)*samples)] = np.array(Z["Z_SFH90_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air)*2)*samples), 0] = float(upper["upper_SFH90_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH90_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air)*2)+1)*samples), ((i + n_SFH06_air + n_SFH05_air)*samples):((i + n_SFH06_air + n_SFH05_air + 1)*samples)] = np.array(Z["Z_SFH90_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH90_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH90_{0}".format(i + 1)][t])
            
            """Photovoltaic consumption constraints. """
            Aub_opt[t+((i + n_SFH06_air + n_SFH05_air + Nc_bui_SS)*samples), ((i + n_SFH06_air + n_SFH05_air)*samples):((i + n_SFH06_air+ n_SFH05_air + 1)*samples)] = I[t,:]*((-1)*(1/COP_SFH90_air[t]))
            Aub_opt[t+((i + n_SFH06_air + n_SFH05_air + Nc_bui_SS)*samples), ((i + n_SFH06_air + n_SFH05_air + n_Bui)*samples):((i + n_SFH06_air + n_SFH05_air + n_Bui + 1)*samples)] = I[t,:]*(P_ele_pv[t])
            Bub_opt[t+((i + n_SFH06_air + n_SFH05_air + Nc_bui_SS)*samples), 0] = 0
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot - 1)*samples), ((i + n_SFH06_air + n_SFH05_air)*samples):((i + n_SFH06_air + n_SFH05_air + 1)*samples)] = I[t,:]*(1/COP_SFH90_air[t])

    for i in range (n_SFH60_air):
        for t in range (samples):
            """1946-1960 archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = np.array(Z["Z_SFH60_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)*samples), 0] = float(upper["upper_SFH60_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH60_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)+1)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = np.array(Z["Z_SFH60_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH60_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH60_{0}".format(i + 1)][t])
            
            """Photovoltaic consumption constraints. """
            Aub_opt[t+((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + Nc_bui_SS)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air+ n_SFH05_air+ n_SFH90_air + 1)*samples)] = I[t,:]*((-1)*(1/COP_SFH60_air[t]))
            Aub_opt[t+((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + Nc_bui_SS)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_Bui)*samples):((i + n_SFH06_air+ n_SFH05_air+ n_SFH90_air + n_Bui + 1)*samples)] = I[t,:]*(P_ele_pv[t])
            Bub_opt[t+((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + Nc_bui_SS)*samples), 0] = 0
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot - 1)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = I[t,:]*(1/COP_SFH60_air[t])

    for i in range (n_SFH06_floor):
        for t in range (samples):
            """2006-today archetype with underfloor heating system indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air+ 1)*samples)] = np.array(Z["Z_SFH06_{0}_floor".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2)*samples), 0] = float(upper["upper_SFH06_{0}_floor".format(i + 1)][t]) - float(Mmax["Mmax_SFH06_{0}_floor".format(i + 1)][t][0])
    
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+1)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air+ 1)*samples)] = np.array(Z["Z_SFH06_{0}_floor".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH06_{0}_floor".format(i + 1)][t][0]) - float(lower["lower_SFH06_{0}_floor".format(i + 1)][t])
            
            """2006-today archetype with underfloor heating system floor temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2 + 2)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air+ 1)*samples)] = np.array(Z_floor["Z_SFH06_{0}_f".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2 + 2)*samples), 0] = float(29) - float(Mmax_floor["Mmax_SFH06_{0}_f".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+ 2 + 1)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = np.array(Z_floor["Z_SFH06_{0}_f".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+ 2 + 1)*samples), 0] = float(Mmax_floor["Mmax_SFH06_{0}_f".format(i+ 1)][t][0]) - float(15)

            """Photovoltaic consumption constraints. """
            Aub_opt[t+((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air + Nc_bui_SS)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air+ n_SFH05_air+ n_SFH90_air + n_SFH60_air+ 1)*samples)] = I[t,:]*((-1)*(1/COP_SFH06_floor[t]))
            Aub_opt[t+((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air + Nc_bui_SS)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air + n_Bui)*samples):((i + n_SFH06_air+ n_SFH05_air+ n_SFH90_air + n_SFH60_air+ n_Bui + 1)*samples)] = I[t,:]*(P_ele_pv[t])
            Bub_opt[t+((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air + Nc_bui_SS)*samples), 0] = 0
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot - 1)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air + 1)*samples)] = I[t,:]*(1/COP_SFH06_floor[t])
    
    """Cluster photovoltaic consumption constraint. """
    for i in range(n_Bui):
        for t in range(0, samples):
            Aub_opt[t+((Nc_tot - Nc_cluster)*samples), ((i + n_Bui)*samples):((i + n_Bui + 1)*samples)] = I[t,:]

    for t in range(0, samples):
        Bub_opt[t+((Nc_tot - Nc_cluster)*samples), 0] = 1 # sum of photovoltaic utilization factors Upv <= 1
    
    """Cluster max electricity constraint. """
    for t in range(0, samples):
        Bub_opt[t+((Nc_tot - 1)*samples), 0] = Eele_max[t]      # constrained electric consumption at cluster level

    """Optimization result (Demand-Response scenario). """
    res = linprog(c=c, A_ub=Aub_opt, b_ub=Bub_opt, bounds=bounds, method='highs-ipm', options={'maxiter': 1000000})
    print("Cluster level:")
    print(res.message + "(Demand-Response scenario)")
    if res.success == False:
        print("\n              ")
        print("\nWARNING")
        print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
        print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
        print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
        print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\nWhereas, the PV utilization factor (decision variable 2) is constrained according to available generation. ")
    
    """Calculation of thermal loads (decision variable x1), PV utilization factor (decision variable x2), 
            indoor air temperatures and building electricity consumption during the DR scenario. """

    Qres_dr   = {}  # Values of the decision variable x1 (thermal demand in Wh)
    Tair_dr   = {}  # Indoor air temeperature (°C) - Demand-Response scenario
    Q_dr   = {}     # Thermal demand (kWh) - Demand-Response scenario
    E_dr   = {}     # Electric demand (kWh) - Demand-Response scenario
    
    Upv   = {}      # Values of the decision variable x2 (photovoltaic utilization factor)
    Ele_pv   = {}   # Photovoltaic consumption (kWh) - Demand-Response scenario
    
    for i in range(n_SFH06_air):
        """Values of the decision variables x1 (i.e., thermal load). """
        Qres_dr["Qres_dr_SFH06_{0}".format(i+1)] = (res.x[(i*samples):((i + 1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH06_{0}".format(i+1)] = list(M["M_SFH06_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}".format(i+1)], Qres_dr["Qres_dr_SFH06_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH06_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH06_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption drawn from the grid (DR scenario). """
        E_dr["E_dr_SFH06_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH06_{0}".format(i+1)][r]/COP_SFH06_air[r] for r in range(0, samples)]
        
        """Values of the decision variables x2 (i.e., PV utilization factor). """
        Upv["Upv_SFH06_{0}".format(i+1)]         = (res.x[((i + n_Bui)*samples):((i + n_Bui + 1)*samples)].reshape(samples,1))
        
        """Calculating the PV consumption. """
        Ele_pv["Ele_pv_SFH06_{0}".format(i+1)]   = [Upv["Upv_SFH06_{0}".format(i+1)][r, 0]*P_ele_pv[r]/1000 for r in range(0, samples)]
    
    for i in range(n_SFH05_air):
        """Values of the decision variables x1 (i.e., thermal load). """
        Qres_dr["Qres_dr_SFH05_{0}".format(i+1)] = (res.x[((i+n_SFH06_air)*samples):((i+n_SFH06_air+1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH05_{0}".format(i+1)] = list(M["M_SFH05_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH05_{0}".format(i+1)], Qres_dr["Qres_dr_SFH05_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH05_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH05_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption drawn from the grid (DR scenario). """
        E_dr["E_dr_SFH05_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH05_{0}".format(i+1)][r]/COP_SFH05_air[r] for r in range(0, samples)]
        
        """Values of the decision variables x2 (i.e., PV utilization factor). """
        Upv["Upv_SFH05_{0}".format(i+1)] = (res.x[((i + n_Bui + n_SFH06_air)*samples):((i + n_Bui + n_SFH06_air + 1)*samples)].reshape(samples,1))
        
        """Calculating the PV consumption. """
        Ele_pv["Ele_pv_SFH05_{0}".format(i+1)]   = [Upv["Upv_SFH05_{0}".format(i+1)][r, 0]*P_ele_pv[r]/1000 for r in range(0, samples)]
        
    for i in range(n_SFH90_air):
        """Values of the decision variables x1 (i.e., thermal load). """
        Qres_dr["Qres_dr_SFH90_{0}".format(i+1)] = (res.x[((i+n_SFH06_air+n_SFH05_air)*samples):((i+n_SFH06_air+n_SFH05_air+1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH90_{0}".format(i+1)] = list(M["M_SFH90_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH90_{0}".format(i+1)], Qres_dr["Qres_dr_SFH90_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH90_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH90_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption drawn from the grid (DR scenario). """
        E_dr["E_dr_SFH90_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH90_{0}".format(i+1)][r]/COP_SFH90_air[r] for r in range(0, samples)]
        
        """Values of the decision variables x2 (i.e., PV utilization factor). """
        Upv["Upv_SFH90_{0}".format(i+1)] = (res.x[((i + n_Bui + n_SFH06_air + n_SFH05_air)*samples):((i + n_Bui + n_SFH06_air + n_SFH05_air + 1)*samples)].reshape(samples,1))
        
        """Calculating the PV consumption. """
        Ele_pv["Ele_pv_SFH90_{0}".format(i+1)]   = [Upv["Upv_SFH90_{0}".format(i+1)][r, 0]*P_ele_pv[r]/1000 for r in range(0, samples)]
    
    for i in range(n_SFH60_air):
        """Values of the decision variables x1 (i.e., thermal load). """
        Qres_dr["Qres_dr_SFH60_{0}".format(i+1)] = (res.x[((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH60_{0}".format(i+1)] = list(M["M_SFH60_{0}".format(i+1)][3, :].reshape(samples, 1) + np.dot(Z["Z_SFH60_{0}".format(i+1)], Qres_dr["Qres_dr_SFH60_{0}".format(i+1)]))
        Q_dr["Q_dr_SFH60_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH60_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption drawn from the grid (DR scenario). """
        E_dr["E_dr_SFH60_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH60_{0}".format(i+1)][r]/COP_SFH60_air[r] for r in range(0, samples)]
        
        """Values of the decision variables x2 (i.e., PV utilization factor). """
        Upv["Upv_SFH60_{0}".format(i+1)] = (res.x[((i + n_Bui + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_Bui + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)].reshape(samples,1))
        
        """Calculating the PV consumption. """
        Ele_pv["Ele_pv_SFH60_{0}".format(i+1)]   = [Upv["Upv_SFH60_{0}".format(i+1)][r, 0]*P_ele_pv[r]/1000 for r in range(0, samples)]

    for i in range(n_SFH06_floor):
        """Values of the decision variables x1 (i.e., thermal load). """
        Qres_dr["Qres_dr_SFH06_{0}_floor".format(i+1)] = (res.x[((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air + 1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH06_{0}_floor".format(i+1)] = list(M["M_SFH06_{0}_floor".format(i+1)][7, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}_floor".format(i+1)], Qres_dr["Qres_dr_SFH06_{0}_floor".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH06_{0}_floor".format(i+1)]       = [Qres_dr["Qres_dr_SFH06_{0}_floor".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption drawn from the grid (DR scenario). """
        E_dr["E_dr_SFH06_{0}_floor".format(i+1)]       = [Q_dr["Q_dr_SFH06_{0}_floor".format(i+1)][r]/COP_SFH06_floor[r] for r in range(0, samples)]
        
        """Values of the decision variables x2 (i.e., PV utilization factor). """
        Upv["Upv_SFH06_{0}_floor".format(i+1)] = (res.x[((i + n_Bui + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*samples):((i + n_Bui + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air + 1)*samples)].reshape(samples,1))
        
        """Calculating the PV consumption. """
        Ele_pv["Ele_pv_SFH06_{0}_floor".format(i+1)]   = [Upv["Upv_SFH06_{0}_floor".format(i+1)][r, 0]*P_ele_pv[r]/1000 for r in range(0, samples)]
        
    """Conversion of dictionary to data frame for calculation of aggregate electric power consumption. """
    E_dr = pd.DataFrame.from_dict(E_dr)
    
    """Calculating electric consumption at the aggregate level during DR scenario. """
    Etot_dr = E_dr.sum(axis=1)
    Edr_sum = float(Etot_dr.sum())
    Etot_dr = [Etot_dr[i] for i in range(0,samples)]  
    
    """Conversion of dictionary to data frame for calculation of aggregate PV consumption. """
    E_pv = pd.DataFrame.from_dict(Ele_pv)
    
    """Calculating PV consumption at the aggregate level. """
    Etot_pv = E_pv.sum(axis=1)
    Epv_dr_sum = float(Etot_pv.sum())
    Etot_pv = [Etot_pv[i] for i in range(0,samples)]  
    
    PVprod =[(P_ele_pv[r])/1000 for r in range(0,samples)]
    
    """Calculation of solar PV consumed during BL scenario useful for comparison. """
    PVdiff_bl   = []
    Etot_pv_BL = []
    for t in range(samples):
        PVdiff_bl.append(PVprod[t]-Etot_bl[t])
        if PVprod[t] > 0:
            if PVdiff_bl[t] > 0: # PV surplus
                Etot_pv_BL.append(Etot_bl[t])
            else:
                Etot_pv_BL.append(PVprod[t])
        else:
            Etot_pv_BL.append(0)
            
    PVprod = np.array(PVprod)
    PVprod_sum = float(PVprod.sum())
    
    Etot_pv_BL = np.array(Etot_pv_BL)
    Epv_bl_sum = float(Etot_pv_BL.sum())
    
    import matplotlib.pyplot as plt
    time = [i*timestep for i in range(0, samples)] # time array
    
    """Do plot. """
    if do_plot_E_dr_cluster:
        """If True, plot cluster electric and PV consumption trends. """
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Time (hours)', fontsize=14)
        ax1.set_ylabel('$\.P_\mathrm{cluster}$ (kW$_\mathrm{el}$)', fontsize=14)
        ax1.plot(time, Etot_bl, 'k', label="BL_scenario", linewidth=1.0)
        ax1.plot(time, Etot_dr, 'r', label="DR_scenario", linewidth=1.0)
        ax1.plot(time, PVprod, 'g--', label="PV_gen", linewidth=1.0)
        ax1.set_xlim(0, duration-1)
        ax1.grid(alpha=0.3)
        plt.legend(frameon=False, fontsize=12)
        fig.set_size_inches(w=10, h=4)
        fig.tight_layout()
        plt.savefig("Eele_cluster", dpi=200)
    
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Time (hours)', fontsize=14)
        ax1.set_ylabel('$\.P_\mathrm{PV cluster}$ (kW$_\mathrm{el}$)', fontsize=14)
        ax1.plot(time, Etot_pv_BL, 'k', label="PV_cons_BL", linewidth=1.0)
        ax1.plot(time, Etot_pv, 'r', label="PV_cons_DR", linewidth=1.0)
        ax1.plot(time, PVprod, 'g--', label="PV_gen", linewidth=1.0)
        ax1.set_xlim(0, duration-1)
        ax1.grid(alpha=0.3)
        plt.legend(frameon=False, fontsize=12)
        fig.set_size_inches(w=10, h=4)
        fig.tight_layout()
        plt.savefig("pv_cluster", dpi=200)
        
        plt.show()
        plt.close('all')
    
    if do_plot_E_dr_single:
        """If True, plot single buildings electric and PV consumption trends. """
        for i in range (n_SFH06_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            ax1.plot(time,PVprod,'g--', linewidth=1, label = "PV_gen")
            ax1.plot(time,E_dr["E_dr_SFH06_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH06_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH06_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH06_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
            
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH05_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            ax1.plot(time,PVprod,'g--', linewidth=1, label = "PV_gen")
            ax1.plot(time,E_dr["E_dr_SFH05_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH05_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH05_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH05_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH90_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            ax1.plot(time,PVprod,'g--', linewidth=1, label = "PV_gen")
            ax1.plot(time,E_dr["E_dr_SFH90_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH90_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH90_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH90_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH60_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            ax1.plot(time,PVprod,'g--', linewidth=1, label = "PV_gen")
            ax1.plot(time,E_dr["E_dr_SFH60_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH60_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH60_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH60_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
    
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH06_floor):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            ax1.plot(time,PVprod,'g--', linewidth=1, label = "PV_gen")
            ax1.plot(time,E_dr["E_dr_SFH06_{0}_floor".format(i+1)],'r', linewidth=1, label = "E_dr_SFH06_{0}_floor".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH06_{0}_floor".format(i+1)],'k', linewidth=1, label = "E_bl_SFH06_{0}_floor".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
            
        plt.show()
        plt.close('all')
        
    if do_plot_Tair_dr_single:
        """If True, plot single buildings indoor air temperature trends. """
        for i in range (n_SFH06_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            ax1.plot(time,Tair_dr["Tair_dr_SFH06_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH06_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH06_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH06_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
            
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH05_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            ax1.plot(time,Tair_dr["Tair_dr_SFH05_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH05_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH05_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH05_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH90_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            ax1.plot(time,Tair_dr["Tair_dr_SFH90_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH90_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH90_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH90_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH60_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            ax1.plot(time,Tair_dr["Tair_dr_SFH60_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH60_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH60_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH60_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
    
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH06_floor):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            ax1.plot(time,Tair_dr["Tair_dr_SFH06_{0}_floor".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH06_{0}_floor".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH06_{0}_floor".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH06_{0}_floor".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
    """Data export to excel. """
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writerEele = pd.ExcelWriter(name_file + '.xlsx', engine='xlsxwriter')

    # Write each dataframe to a different worksheet.
    E_bl.to_excel(writerEele, sheet_name='Eele_BL')
    E_dr.to_excel(writerEele, sheet_name='Eele_DR')
    E_pv.to_excel(writerEele, sheet_name='Epv_cons_DR')

    # Close the Pandas Excel writer and output the Excel file.
    writerEele.close()
    
    """Electricity consumption during BL and DR scenarios. """
    print("\nCluster electricity consumption of about " +
          str(round(Ebl_sum,2)) + " kWh during baseline (BL) scenario")
    print("\n------------------------------------------------------------------------")
    
    Diff   = ((Edr_sum - Ebl_sum)/Ebl_sum)*100
    
    if Diff >= 0:
        print("\nCluster electricity consumption of about " +
              str(round(Edr_sum,2)) + " kWh during DR scenario with centralized PV generation,"
              "\nwith a " + str(round(Diff,2)) + "% increase compared to BL scenario")
        print("\n------------------------------------------------------------------------")
    else:
        print("\nCluster electricity consumption of about " +
              str(round(Edr_sum,2)) + " kWh during DR scenario with centralized PV generation,"
              "\nwith a " + str(round(Diff*(-1),2)) + "% reduction compared to BL scenario")
        print("\n------------------------------------------------------------------------")
    
    """PV consumption during BL and DR scenarios. """
    print("\nCluster PV consumption of about " +
          str(round(Epv_bl_sum,2)) + " kWh during baseline (BL) scenario" + 
          "\n(consumption of " + str(round(Epv_bl_sum/PVprod_sum*100,2)) + "%)")
    print("\n------------------------------------------------------------------------")
   
    Diff_pv = ((Epv_dr_sum-Epv_bl_sum)/Epv_bl_sum)*100
    
    if Diff_pv >= 0:
        print("\nCluster PV consumption of about " +
              str(round(Epv_dr_sum,2)) + " kWh during DR scenario with centralized PV generation" + 
              "\n(consumption of " + str(round(Epv_dr_sum/PVprod_sum*100,2)) + "%)," 
              "\nwith a " + str(round(Diff_pv,2)) + "% increase compared to BL scenario")
        print("\n------------------------------------------------------------------------")
    else:
        print("\nCluster PV consumption of about " +
              str(round(Epv_dr_sum,2)) + " kWh during DR scenario with centralized PV generation" + 
              "\n(consumption of " + str(round(Epv_dr_sum/PVprod_sum*100,2)) + "%)," 
              "\nwith a " + str(round(Diff_pv,2)) + "% decrease compared to BL scenario")
        print("\n------------------------------------------------------------------------")
    
    
    return(Etot_dr)

def pv_distributed(Locality, day, month, duration, timestep, Thermal_load, 
                            n_SFH06_air, n_SFH05_air, n_SFH90_air, n_SFH60_air,n_SFH06_floor, 
                            HP_data_type, 
                            Supply_air, rated_cap_06_air, rated_cap_05_air, rated_cap_90_air, rated_cap_60_air, rated_COP_06_air, rated_COP_05_air, rated_COP_90_air, rated_COP_60_air,
                            Supply_floor, rated_cap_06_floor, rated_COP_06_floor,
                            Setpoint_type, setpoint, th_tolerance_BL_sup_air, th_tolerance_BL_min_air, th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor,
                            th_tolerance_DR_sup_air, th_tolerance_DR_min_air, th_tolerance_DR_sup_floor, th_tolerance_DR_min_floor,
                            pv_size, panel_dim, pv_rp_06_air, pv_rp_05_air, pv_rp_90_air, pv_rp_60_air, pv_rp_06_floor, f_limit,
                            do_plot_E_dr_single, do_plot_E_dr_cluster,do_plot_Tair_dr_single, name_file):
    """pv_distributed is a function based on linear programming that minimizes the consumption of building electricity drawn from the grid 
            during the application of a strategy of load shifting to distributed photovoltaic generation (i.e., each household is provided with its own PV generation). """
    
    # number of buildings:
    n_SFH06_air = int(n_SFH06_air)
    n_SFH05_air = int(n_SFH05_air)
    n_SFH90_air = int(n_SFH90_air)
    n_SFH60_air = int(n_SFH60_air)
    
    if Thermal_load == "heating":
        n_SFH06_floor = int(n_SFH06_floor) # Radiant floor emission system for the SFH06 archetype (age class 2006-today) only during heating
    else: 
        n_SFH06_floor = 0
    
    n_Bui       = n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air + n_SFH06_floor # Total number of buildings composing the cluster
    
    samples = int(duration/timestep) # Calculation of total steps

    T_est = MeteoFile(Locality, day, month, duration, timestep)[0] # Outdoor air temperature (°C)
    
    """Call of the performance characteristics of heat pumps based on the definition of heat loads 
           and attribute definition for calling temperature setpoint functions"""
    if Thermal_load == 'cooling':
        """Call performance characteristics and capacity of heat pumps. """
        if n_SFH06_air > 0:
            COP_SFH06_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH06_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH05_air > 0:
            COP_SFH05_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH05_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH90_air > 0:
            COP_SFH90_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH90_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH60_air > 0:
            COP_SFH60_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH60_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[1] for i in range(0,samples)]
            
        """Thermostat mode setting in space cooling"""
        thermostat = 'th_cooling'
    elif Thermal_load == 'heating':
        if n_SFH06_air > 0:
            COP_SFH06_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH06_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH05_air > 0:
            COP_SFH05_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH05_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH90_air > 0:
            COP_SFH90_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH90_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH60_air > 0:
            COP_SFH60_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH60_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH06_floor > 0:
            COP_SFH06_floor = [HPperformance.HP_heating(T_est[i],Supply_floor,rated_cap_06_air,rated_COP_06_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH06_floor = [HPperformance.HP_heating(T_est[i],Supply_floor,rated_cap_06_air,rated_COP_06_air,HP_data_type)[1] for i in range(0,samples)]
        
        """Thermostat mode setting in space heating"""
        thermostat = 'th_heating'
    else:
        print('Thermal load not available')
    
    """Definition of occupancy profiles, relative internal gains (people, equipment and lighting) and setpoint temperature for each building 
           (useful for diversification of consumption loads within the cluster)"""
    occ   = {} # occupancy
    gains = {} # internal gains (W)
    tsp   = {} # thermostat set-points (°C)
    
    if Setpoint_type == "pre_defined":
        """Pre-defined temperature set-point definition
               For adding profiles or modifying existing ones, go to the User_pattern.py file and see the functions 
               th_heating (set-points for space heating) or th_cooling (set-points for space cooling)"""
        for i in range(n_SFH06_air):
            tsp["tspSFH06_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i]
        for i in range(n_SFH05_air):
            tsp["tspSFH05_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air]
        for i in range(n_SFH90_air):
            tsp["tspSFH90_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air]
        for i in range(n_SFH60_air):
            tsp["tspSFH60_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air]
        for i in range(n_SFH06_floor):
            tsp["tspSFH06_{0}_floor".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air]
    
    elif Setpoint_type == "user_defined":
        """User-defined temperature set-point definition"""
        t_setpoint = []
        for i in range (samples):
            t_setpoint.append(setpoint)
        for i in range(n_SFH06_air):
            tsp["tspSFH06_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH05_air):
            tsp["tspSFH05_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH90_air):
            tsp["tspSFH90_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH60_air):
            tsp["tspSFH60_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH06_floor):
            tsp["tspSFH06_{0}_floor".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
    
    """Occupancy patterns and related internal gain profiles definition"""
    for i in range(n_SFH06_air):
        occ["occupancySFH06_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i]
        gains["gainsSFH06_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i]
    for i in range(n_SFH05_air):
        occ["occupancySFH05_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air]
        gains["gainsSFH05_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air]
    for i in range(n_SFH90_air):
        occ["occupancySFH90_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air]
        gains["gainsSFH90_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air]
    for i in range(n_SFH60_air):
        occ["occupancySFH60_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air]
        gains["gainsSFH60_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air]
    for i in range(n_SFH06_floor):
        occ["occupancySFH06_{0}_floor".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air]
        gains["gainsSFH06_{0}_floor".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air]

    """In the following, thermal loads provided as input to the state-space model formulation are calculated via linear programming . 
           Specifically, the thermal load of each building is minimized to maintain the indoor temperature set-point.
           Linear programming solves problems of the following form: 
               minimize c @ x
               such that:
               A_ub @ x <= b_ub
               A_eq @ x == b_eq
               lb <= x <= ub
               
               For more information see the documentation of the python library scipy.optimize.linprog"""
           
    c      = {}     # Coefficients of the linear objective function to be minimized
    Aub    = {}     # Matrix of inequality constraints
    Bub    = {}     # Vector of inequality constraints
    bounds = {}     # Sequence of pairs (min, max) for each element in x, defining the minimum and maximum values of the decision variable
    
    Z      = {}     # Variable describing the state of the system containing the decision variable (e.g., thermal demand)
    M      = {}     # Variable describing the state of the indoor air temperture
    Mmax      = {}  # Variable describing the state of the indoor air temperture
    res_opt   = {}  # Optimization result - baseline scenario
    Qres_bl   = {}  # Values of the decision variable (thermal demand in Wh)
    Tair_bl   = {}  # Indoor air temeperature (°C) - baseline scenario
    Q_bl   = {}     # Thermal demand (kWh) - baseline scenario
    E_bl   = {}     # Electric demand (kWh) - baseline scenario
    upper   = {}    # Upper range thermostat set-point profile (°C) - baseline scenario
    lower   = {}    # Lower range thermostat set-point profile (°C) - baseline scenario
    
    Mmax_floor = {} # Variable describing the state of the floor temperture
    Z_floor    = {} # Variable describing the state of the floor system containing the decision variable (e.g., thermal demand)

    """Calling archetype (SFHs) functions from the "Archetypes.py" module. """
    for i in range(n_SFH06_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH06 buildings (age class 2006-today) with air emission system."""
        c["c_SFH06_{0}".format(i+1)]             = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[0]
        Aub["Aub_SFH06_{0}".format(i+1)]         = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[1]
        Bub["Bub_SFH06_{0}".format(i+1)]         = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[2]
        bounds["bounds_SFH06_{0}".format(i+1)]   = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH06_{0}".format(i+1)]             = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[4]
        Z["Z_SFH06_{0}".format(i+1)]             = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[5]
        Mmax["Mmax_SFH06_{0}".format(i+1)]       = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH06_{0}".format(i+1)]     = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[7]
        lower["lower_SFH06_{0}".format(i+1)]     = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH06_{0}".format(i+1)] = linprog(c=c["c_SFH06_{0}".format(i+1)], A_ub=Aub["Aub_SFH06_{0}".format(i+1)], b_ub=Bub["Bub_SFH06_{0}".format(i+1)],bounds=bounds["bounds_SFH06_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 2006-today building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH06_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH06_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH06_{0}".format(i+1)] = (res_opt["res_opt_SFH06_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH06_{0}".format(i+1)] = list(M["M_SFH06_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}".format(i+1)], Qres_bl["Qres_bl_SFH06_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH06_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH06_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH06_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH06_{0}".format(i+1)][r]/COP_SFH06_air[r] for r in range(0, samples)]
    for i in range(n_SFH05_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH05 buildings (age class 1991-2005) with air emission system."""
        c["c_SFH05_{0}".format(i+1)]             = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[0]
        Aub["Aub_SFH05_{0}".format(i+1)]         = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[1]
        Bub["Bub_SFH05_{0}".format(i+1)]         = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[2]
        bounds["bounds_SFH05_{0}".format(i+1)]   = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH05_{0}".format(i+1)]             = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[4]
        Z["Z_SFH05_{0}".format(i+1)]             = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[5]
        Mmax["Mmax_SFH05_{0}".format(i+1)]       = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH05_{0}".format(i+1)] = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[7]
        lower["lower_SFH05_{0}".format(i+1)] = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH05_{0}".format(i+1)] = linprog(c=c["c_SFH05_{0}".format(i+1)], A_ub=Aub["Aub_SFH05_{0}".format(i+1)], b_ub=Bub["Bub_SFH05_{0}".format(i+1)],bounds=bounds["bounds_SFH05_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 1990-2005 building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH05_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH05_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH05_{0}".format(i+1)] = (res_opt["res_opt_SFH05_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH05_{0}".format(i+1)] = list(M["M_SFH05_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH05_{0}".format(i+1)], Qres_bl["Qres_bl_SFH05_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH05_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH05_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH05_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH05_{0}".format(i+1)][r]/COP_SFH05_air[r] for r in range(0, samples)]
    for i in range(n_SFH90_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH90 buildings (age class 1976-1990) with air emission system."""
        c["c_SFH90_{0}".format(i+1)]             = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[0]
        Aub["Aub_SFH90_{0}".format(i+1)]         = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[1]
        Bub["Bub_SFH90_{0}".format(i+1)]         = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[2]
        bounds["bounds_SFH90_{0}".format(i+1)]   = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH90_{0}".format(i+1)]             = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[4]
        Z["Z_SFH90_{0}".format(i+1)]             = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[5]
        Mmax["Mmax_SFH90_{0}".format(i+1)]       = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH90_{0}".format(i+1)]     = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[7]
        lower["lower_SFH90_{0}".format(i+1)]     = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[8]
        
        """Optimization result (baseline scenario)"""
        res_opt["res_opt_SFH90_{0}".format(i+1)] = linprog(c=c["c_SFH90_{0}".format(i+1)], A_ub=Aub["Aub_SFH90_{0}".format(i+1)], b_ub=Bub["Bub_SFH90_{0}".format(i+1)],bounds=bounds["bounds_SFH90_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 1976-1990 building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH90_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH90_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH90_{0}".format(i+1)] = (res_opt["res_opt_SFH90_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH90_{0}".format(i+1)] = list(M["M_SFH90_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH90_{0}".format(i+1)], Qres_bl["Qres_bl_SFH90_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH90_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH90_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH90_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH90_{0}".format(i+1)][r]/COP_SFH90_air[r] for r in range(0, samples)]
    for i in range(n_SFH60_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH60 buildings (age class 1946-1960) with air emission system."""
        c["c_SFH60_{0}".format(i+1)]             = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[0]
        Aub["Aub_SFH60_{0}".format(i+1)]         = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[1]
        Bub["Bub_SFH60_{0}".format(i+1)]         = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[2]
        bounds["bounds_SFH60_{0}".format(i+1)]   = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH60_{0}".format(i+1)]             = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[4]
        Z["Z_SFH60_{0}".format(i+1)]             = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[5]
        Mmax["Mmax_SFH60_{0}".format(i+1)]       = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH60_{0}".format(i+1)] = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[7]
        lower["lower_SFH60_{0}".format(i+1)] = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH60_{0}".format(i+1)] = linprog(c=c["c_SFH60_{0}".format(i+1)], A_ub=Aub["Aub_SFH60_{0}".format(i+1)], b_ub=Bub["Bub_SFH60_{0}".format(i+1)],bounds=bounds["bounds_SFH60_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 1946-1960 building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH60_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH60_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")

        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH60_{0}".format(i+1)] = (res_opt["res_opt_SFH60_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH60_{0}".format(i+1)] = list(M["M_SFH60_{0}".format(i+1)][3, :].reshape(samples, 1) + np.dot(Z["Z_SFH60_{0}".format(i+1)], Qres_bl["Qres_bl_SFH60_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH60_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH60_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH60_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH60_{0}".format(i+1)][r]/COP_SFH60_air[r] for r in range(0, samples)]
    for i in range(n_SFH06_floor):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH06 buildings (age class 2006-today) with underfloor heating system."""
        c["c_SFH06_{0}_floor".format(i+1)]             = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[0]
        Aub["Aub_SFH06_{0}_floor".format(i+1)]         = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[1]
        Bub["Bub_SFH06_{0}_floor".format(i+1)]         = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[2]
        bounds["bounds_SFH06_{0}_floor".format(i+1)]   = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[3]
        
        """Calling of variables referring to the state of indoor air temperature node. """
        M["M_SFH06_{0}_floor".format(i+1)]             = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[4]
        Z["Z_SFH06_{0}_floor".format(i+1)]             = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[5]
        Mmax["Mmax_SFH06_{0}_floor".format(i+1)]       = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[6]
        
        """Calling of variables referring to the state of floor node. """
        Z_floor["Z_SFH06_{0}_f".format(i+1)]           = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[9]
        Mmax_floor["Mmax_SFH06_{0}_f".format(i+1)]     = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[10]
        
        """Calling of variables referring to the state of indoor air temperature. """
        upper["upper_SFH06_{0}_floor".format(i+1)]     = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[7]
        lower["lower_SFH06_{0}_floor".format(i+1)]     = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH06_{0}_floor".format(i+1)] = linprog(c=c["c_SFH06_{0}_floor".format(i+1)], A_ub=Aub["Aub_SFH06_{0}_floor".format(i+1)], b_ub=Bub["Bub_SFH06_{0}_floor".format(i+1)],bounds=bounds["bounds_SFH06_{0}_floor".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype building {0}".format(i+1) + " with underfloor heating system:")
        print(res_opt["res_opt_SFH06_{0}_floor".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH06_{0}_floor".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH06_{0}_floor".format(i+1)] = (res_opt["res_opt_SFH06_{0}_floor".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH06_{0}_floor".format(i+1)] = list(M["M_SFH06_{0}_floor".format(i+1)][7, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}_floor".format(i+1)], Qres_bl["Qres_bl_SFH06_{0}_floor".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH06_{0}_floor".format(i+1)]       = [Qres_bl["Qres_bl_SFH06_{0}_floor".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH06_{0}_floor".format(i+1)]       = [Q_bl["Q_bl_SFH06_{0}_floor".format(i+1)][r]/COP_SFH06_floor[r] for r in range(0, samples)]
    
    """Conversion of dictionary to data frame for calculation of aggregate electric power consumption. """
    E_bl = pd.DataFrame.from_dict(E_bl)
    
    """Calculating base electrical loads at the aggregate level. """
    Etot_bl = E_bl.sum(axis=1)
    Ebl_sum = float(Etot_bl.sum())
    Etot_bl = [Etot_bl[i] for i in range(0,samples)]  
    
    """Definition of distributed photovoltaic generation profiles by applying the python pvlib library. 
            Specifically, this distributed generation scenario represents the scenario in which each building is provided with its own photovoltaic system. 
            As a result, solar PV that is not self-consumed is released to the power grid. """
    
    start      = MeteoFile(Locality,int(day),int(month),duration,timestep)[10] # start time to give as input to the PV_load function
    
    if pv_size == "calculated":
        """If "calculated", the rated power of the PV systems available to archetypes are calculated according to the Italian regulation Dlgs. 28/2011. """
        panel = float(panel_dim)                                  # W
     
        PV_06_air    = (Archetypes.floor_area_SFH2006/50)*1000    # W
        PV_05_air    = (Archetypes.floor_area_SFH2005/50)*1000    # W
        PV_90_air    = (Archetypes.floor_area_SFH1990/50)*1000    # W
        PV_60_air    = (Archetypes.floor_area_SFH1960/50)*1000    # W
        PV_06_floor  = (Archetypes.floor_area_SFH2006/50)*1000    # W
        
        n_panels_06_air   = math.ceil(PV_06_air/panel)            # Number of panels defining the photovoltaic system of the 2006-today archetype with air system
        n_panels_05_air   = math.ceil(PV_05_air/panel)            # Number of panels defining the photovoltaic system of the 1991-2005 archetype with air system
        n_panels_90_air   = math.ceil(PV_90_air/panel)            # Number of panels defining the photovoltaic system of the 1976-1990 archetype with air system
        n_panels_60_air   = math.ceil(PV_60_air/panel)            # Number of panels defining the photovoltaic system of the 1946-1960 archetype with air system
        n_panels_06_floor = math.ceil(PV_06_floor/panel)          # Number of panels defining the photovoltaic system of the 2006-today archetype with underfloor heating system
        
        PV_rp_06_air   = n_panels_06_air*panel/1000               # PV array's rated power (kW)
        PV_rp_05_air   = n_panels_05_air*panel/1000               # PV array's rated power (kW)
        PV_rp_90_air   = n_panels_90_air*panel/1000               # PV array's rated power (kW)
        PV_rp_60_air   = n_panels_60_air*panel/1000               # PV array's rated power (kW)
        PV_rp_06_floor = n_panels_06_floor*panel/1000             # PV array's rated power (kW)
    elif pv_size == "user_defined":
        """If "user-defined", the rated power of the PV systems are set by the user. """
        PV_rp_06_air   = pv_rp_06_air   # PV array's rated power (kW)
        PV_rp_05_air   = pv_rp_05_air   # PV array's rated power (kW)
        PV_rp_90_air   = pv_rp_90_air   # PV array's rated power (kW)
        PV_rp_60_air   = pv_rp_60_air   # PV array's rated power (kW)
        PV_rp_06_floor = pv_rp_06_floor # PV array's rated power (kW)
    
    """Calculation of distributed photovoltaic generation profiles. """
    P_ele_06_air   = PV_load(PV_rp_06_air, start, duration, timestep, Locality)   # (W)
    P_ele_05_air   = PV_load(PV_rp_05_air, start, duration, timestep, Locality)   # (W)
    P_ele_90_air   = PV_load(PV_rp_90_air, start, duration, timestep, Locality)   # (W)
    P_ele_60_air   = PV_load(PV_rp_60_air, start, duration, timestep, Locality)   # (W)
    P_ele_06_floor = PV_load(PV_rp_06_floor, start, duration, timestep, Locality) # (W)
    
    """------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    The Demand-Response (DR) event is modeled below. 
    From the definition of PV generation profiles and electricity demand during the baseline scenario, electricity demand 
    during the DR scenario is managed by shifting building loads to their own PV generation hours.
    Specifically, electricity for each building drawn from the grid is minimized:
        min (Q_dr/COP - Upv * P_ele_pv), 
        with
        Q_dr = thermal load (variable 1),
        Upv  = PV utilization factor (variable 2).
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
    """For each simulation day, calculation of the aggregate peak load. """
    n_days        = int(duration/24) # number of simulation days
    
    Max         = []
    Max_dr      = []

    for n in range (n_days):
        Max.append(np.amax(np.array(Etot_bl[int(samples/n_days)*n:int(samples/n_days)*(n + 1)])))
        for i in range (int(samples/n_days)*n,int(samples/n_days)*(n + 1)):
            Max_dr.append(Max[n]) # maximum electricity consumption during the baseline scenario

    """Definition of the electricity demand of the cluster of buildings constrained to remain below a certain limit to reduce rebound effects.. """
    Eele_max = []

    for t in range(samples):
            Eele_max.append(float(Max_dr[t]*1000)*(float(f_limit)))
    
    """Definition of temperature profiles during the DR event to unlock building energy flexibility. """
    if Thermal_load == 'cooling':
        for i in range(n_SFH06_air):
            for t in range(samples):
                if P_ele_06_air[t] > 0:
                    lower["lower_SFH06_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to increase PV consumption
                else:
                    upper["upper_SFH06_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for cooling to reduce electric consumption from the grid
        for i in range(n_SFH05_air):
            for t in range(samples):
                if P_ele_05_air[t] > 0:
                    lower["lower_SFH05_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to increase PV consumption
                else:
                    upper["upper_SFH05_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for cooling to reduce electric consumption from the grid
        for i in range(n_SFH90_air):
            for t in range(samples):
                if P_ele_90_air[t] > 0:
                    lower["lower_SFH90_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to increase PV consumption
                else:
                    upper["upper_SFH90_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for cooling to reduce electric consumption from the grid
        for i in range(n_SFH60_air):
            for t in range(samples):
                if P_ele_60_air[t] > 0:
                    lower["lower_SFH60_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to increase PV consumption
                else:
                    upper["upper_SFH60_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for cooling to reduce electric consumption from the grid

    elif Thermal_load == 'heating':
        for i in range(n_SFH06_air):
            for t in range(samples):
                if P_ele_06_air[t] > 0:
                    upper["upper_SFH06_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for heating to increase PV consumption
                else:
                    lower["lower_SFH06_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to reduce electric consumption from the grid
        for i in range(n_SFH05_air):
            for t in range(samples):
                if P_ele_05_air[t] > 0:
                    upper["upper_SFH05_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for heating to increase PV consumption
                else:
                    lower["lower_SFH05_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to reduce electric consumption from the grid
        for i in range(n_SFH90_air):
            for t in range(samples):
                if P_ele_90_air[t] > 0:
                    upper["upper_SFH90_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for heating to increase PV consumption
                else:
                    lower["lower_SFH90_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to reduce electric consumption from the grid
        for i in range(n_SFH60_air):
            for t in range(samples):
                if P_ele_60_air[t] > 0:
                    upper["upper_SFH60_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for heating to increase PV consumption
                else:
                    lower["lower_SFH60_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to reduce electric consumption from the grid
        for i in range(n_SFH06_floor):
            for t in range(samples):
                if P_ele_06_floor[t] > 0:
                    upper["upper_SFH06_{0}_floor".format(i+1)][t] += float(th_tolerance_DR_sup_floor) # Increase in indoor air temperature set point for heating to increase PV consumption
                else:
                    lower["lower_SFH06_{0}_floor".format(i+1)][t] -= float(th_tolerance_DR_min_floor) # Reduction of indoor air temperature set point for cooling to reduce electric consumption from the grid
    
    """In the following, thermal loads provided as input to the state-space model formulation are calculated via linear programming. 
           Specifically, the thermal load of each building is minimized to maintain the indoor temperature set-point, 
           and the electricity demand of the cluster of buildings is constrained to remain below a certain 
           percentage reduction (Eele_max) for the duration of the peak-shaving event (duration_dr).
           Linear programming solves problems of the following form: 
               minimize c @ x
               such that:
               A_ub @ x <= b_ub
               A_eq @ x == b_eq
               lb <= x <= ub
               
               For more information see the documentation of the python library scipy.optimize.linprog"""

    """Number of decision varibales.
            variable 1 = thermal load
            variable 2 = PV utilization factor"""
    Nv = 2 # number of decision variables
    
    """Definition of the number of constraints at the individual building and cluster level. """
    Nc_SFH06_air_SS  = (2)*n_SFH06_air          # 2 for indoor air temperature range 
    Nc_SFH06_air_PV  = (1)*n_SFH06_air          # 1 for PV consumption

    Nc_SFH05_air_SS  = (2)*n_SFH05_air          # 2 for indoor air temperature range 
    Nc_SFH05_air_PV  = (1)*n_SFH05_air          # 1 for PV consumption

    Nc_SFH90_air_SS  = (2)*n_SFH90_air          # 2 for indoor air temperature range 
    Nc_SFH90_air_PV  = (1)*n_SFH90_air          # 1 for PV consumption

    Nc_SFH60_air_SS  = (2)*n_SFH60_air          # 2 for indoor air temperature range 
    Nc_SFH60_air_PV  = (1)*n_SFH60_air          # 1 for PV consumption

    Nc_06_floor_SS   = (4)*n_SFH06_floor        # 2 for indoor air temperature range + 2 for temperature range of underfloor heating system
    Nc_06_floor_PV   = (1)*n_SFH06_floor        # 1 for PV consumption

    Nc_bui_SS     = Nc_SFH06_air_SS + Nc_SFH05_air_SS + Nc_SFH90_air_SS + Nc_SFH60_air_SS + Nc_06_floor_SS 
    Nc_bui_PV     = Nc_SFH06_air_PV + Nc_SFH05_air_PV + Nc_SFH90_air_PV + Nc_SFH60_air_PV + Nc_06_floor_PV 
    
    Nc_cluster = 1 # max electricity consumption at cluster level
    Nc_tot     = int(Nc_bui_SS + Nc_bui_PV + Nc_cluster)  # Total number of constrains
    
    bounds = []                         # Sequence of pairs (min, max) for each element in x, defining the minimum and maximum values of the decision variable
    c = np.zeros([1,n_Bui*Nv*samples])  # Coefficients of the linear objective function to be minimized (electric consumption)
    
    """VARIABLE 1: thermal loads. """
    """Setting the bounds and coefficients of 2006-today archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH06_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH06_air[r]))              # bounds of the decision variable (thermal power, Qth)
            c[0,r+(i*samples)] = float(1/COP_SFH06_air[r]) # having x1 = Qth and c = 1/COP minimizes the power consumption
    
    """Setting the bounds and coefficients of 1991-2005 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH05_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH05_air[r]))                          # bounds of the decision variable (thermal power, Qth)
            c[0,r+(i+n_SFH06_air)*samples] = float(1/COP_SFH05_air[r]) # having x1 = Qth and c = 1/COP minimizes the power consumption
    
    """Setting the bounds and coefficients of 1976-1990 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH90_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH90_air[r]))                                      # bounds of the decision variable (thermal power, Qth)
            c[0,r+(i+n_SFH06_air+n_SFH05_air)*samples] = float(1/COP_SFH90_air[r]) # having x1 = Qth and c = 1/COP minimizes the power consumption
    
    """Setting the bounds and coefficients of 1946-1960 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH60_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH60_air[r]))                                                  # bounds of the decision variable (thermal power, Qth)
            c[0,r+(i+n_SFH06_air+n_SFH05_air+n_SFH90_air)*samples] = float(1/COP_SFH60_air[r]) # having x1 = Qth and c = 1/COP minimizes the power consumption
    
    """Setting the bounds and coefficients of 2006-today archetypes with underfloor heating system. """
    for i in range(0,n_SFH06_floor):
        for r in range(0,samples):
            bounds.append((0,Q_SFH06_floor[r]))                                                              # bounds of the decision variable (thermal power, Qth)
            c[0,r+(i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air)*samples] = float(1/COP_SFH06_floor[r]) # having x1 = Qth and c = 1/COP minimizes the power consumption
    
    """VARIABLE 2: PV utilization factor. """
    """Setting the bounds and coefficients of 2006-today archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH06_air):
        for r in range(0,samples):
            bounds.append((0,1)) # bounds of PV 
            c[0,r+((i+n_Bui)*samples)] = float((-1)*P_ele_06_air[r])
    
    """Setting the bounds and coefficients of 1991-2005 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH05_air):
        for r in range(0,samples):
            bounds.append((0,1)) # bounds of PV 
            c[0,r+(i+n_SFH06_air+n_Bui)*samples] = float((-1)*P_ele_05_air[r])
    
    """Setting the bounds and coefficients of 1976-1990 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH90_air):
        for r in range(0,samples):
            bounds.append((0,1)) # bounds of PV 
            c[0,r+(i+n_SFH06_air+n_SFH05_air+n_Bui)*samples] = float((-1)*P_ele_90_air[r])
    
    """Setting the bounds and coefficients of 1946-1960 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH60_air):
        for r in range(0,samples):
            bounds.append((0,1)) # bounds of PV 
            c[0,r+(i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_Bui)*samples] = float((-1)*P_ele_60_air[r])
    
    """Setting the bounds and coefficients of 2006-today archetypes with underfloor heating system. """
    for i in range(0,n_SFH06_floor):
        for r in range(0,samples):
            bounds.append((0,1))
            c[0,r+(i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air+n_Bui)*samples] = float((-1)*P_ele_06_floor[r])

    I         = np.identity(samples)
    
    """Setting matrix and vector size of inequality constraints. """
    Aub_opt = np.zeros([(Nc_tot)*samples,n_Bui*samples*Nv])
    Bub_opt = np.empty([(Nc_tot)*samples,1])

    """Compilation of the matrix and vector of inequality constraints. """
    for i in range (n_SFH06_air):
        for t in range (samples):
            """2006-today archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i)*2)*samples), ((i)*samples):((i + 1)*samples)] = np.array(Z["Z_SFH06_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i)*2)*samples), 0] = float(upper["upper_SFH06_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH06_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i)*2)+1)*samples), ((i)*samples):((i + 1)*samples)] = np.array(Z["Z_SFH06_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH06_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH06_{0}".format(i + 1)][t])
            
            """Photovoltaic consumption constraints. """
            Aub_opt[t+((i + Nc_bui_SS)*samples), ((i)*samples):((i + 1)*samples)] = I[t,:]*((-1)*(1/COP_SFH06_air[t]))
            Aub_opt[t+((i + Nc_bui_SS)*samples), ((i + n_Bui)*samples):((i + n_Bui + 1)*samples)] = I[t,:]*(P_ele_06_air[t])
            Bub_opt[t+((i + Nc_bui_SS)*samples), 0] = 0
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot-Nc_cluster)*samples), ((i)*samples):((i + 1)*samples)] = I[t,:]*(1/COP_SFH06_air[t])
    
    for i in range (n_SFH05_air):
        for t in range (samples):
            """1991-2005 archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air)*2)*samples), ((i + n_SFH06_air)*samples):((i + n_SFH06_air + 1)*samples)] = np.array(Z["Z_SFH05_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air)*2)*samples), 0] = float(upper["upper_SFH05_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH05_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air)*2)+1)*samples), ((i + n_SFH06_air)*samples):((i + n_SFH06_air + 1)*samples)] = np.array(Z["Z_SFH05_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH05_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH05_{0}".format(i + 1)][t])
            
            """Photovoltaic consumption constraints. """
            Aub_opt[t+((i + n_SFH06_air + Nc_bui_SS)*samples), ((i + n_SFH06_air)*samples):((i + n_SFH06_air + 1)*samples)] = I[t,:]*((-1)*(1/COP_SFH05_air[t]))
            Aub_opt[t+((i + n_SFH06_air + Nc_bui_SS)*samples), ((i + n_SFH06_air + n_Bui)*samples):((i + n_SFH06_air + n_Bui + 1)*samples)] = I[t,:]*(P_ele_05_air[t])
            Bub_opt[t+((i + n_SFH06_air + Nc_bui_SS)*samples), 0] = 0
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot-Nc_cluster)*samples), ((i + n_SFH06_air)*samples):((i + n_SFH06_air + 1)*samples)] = I[t,:]*(1/COP_SFH05_air[t])

    for i in range (n_SFH90_air):
        for t in range (samples):
            """1976-1990 archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air)*2)*samples), ((i + n_SFH06_air + n_SFH05_air)*samples):((i + n_SFH06_air + n_SFH05_air + 1)*samples)] = np.array(Z["Z_SFH90_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air)*2)*samples), 0] = float(upper["upper_SFH90_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH90_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air)*2)+1)*samples), ((i + n_SFH06_air + n_SFH05_air)*samples):((i + n_SFH06_air + n_SFH05_air + 1)*samples)] = np.array(Z["Z_SFH90_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH90_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH90_{0}".format(i + 1)][t])
            
            """Photovoltaic consumption constraints. """
            Aub_opt[t+((i + n_SFH06_air + n_SFH05_air + Nc_bui_SS)*samples), ((i + n_SFH06_air + n_SFH05_air)*samples):((i + n_SFH06_air+ n_SFH05_air + 1)*samples)] = I[t,:]*((-1)*(1/COP_SFH90_air[t]))
            Aub_opt[t+((i + n_SFH06_air + n_SFH05_air + Nc_bui_SS)*samples), ((i + n_SFH06_air + n_SFH05_air + n_Bui)*samples):((i + n_SFH06_air + n_SFH05_air + n_Bui + 1)*samples)] = I[t,:]*(P_ele_90_air[t])
            Bub_opt[t+((i + n_SFH06_air + n_SFH05_air + Nc_bui_SS)*samples), 0] = 0
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot-Nc_cluster)*samples), ((i + n_SFH06_air + n_SFH05_air)*samples):((i + n_SFH06_air + n_SFH05_air + 1)*samples)] = I[t,:]*(1/COP_SFH90_air[t])

    for i in range (n_SFH60_air):
        for t in range (samples):
            """1946-1960 archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = np.array(Z["Z_SFH60_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)*samples), 0] = float(upper["upper_SFH60_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH60_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)+1)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = np.array(Z["Z_SFH60_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH60_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH60_{0}".format(i + 1)][t])
            
            """Photovoltaic consumption constraints. """
            Aub_opt[t+((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + Nc_bui_SS)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air+ n_SFH05_air+ n_SFH90_air + 1)*samples)] = I[t,:]*((-1)*(1/COP_SFH60_air[t]))
            Aub_opt[t+((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + Nc_bui_SS)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_Bui)*samples):((i + n_SFH06_air+ n_SFH05_air+ n_SFH90_air + n_Bui + 1)*samples)] = I[t,:]*(P_ele_60_air[t])
            Bub_opt[t+((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + Nc_bui_SS)*samples), 0] = 0
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot-Nc_cluster)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = I[t,:]*(1/COP_SFH60_air[t])

    for i in range (n_SFH06_floor):
        for t in range (samples):
            """2006-today archetype with underfloor heating system indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air+ 1)*samples)] = np.array(Z["Z_SFH06_{0}_floor".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2)*samples), 0] = float(upper["upper_SFH06_{0}_floor".format(i + 1)][t]) - float(Mmax["Mmax_SFH06_{0}_floor".format(i + 1)][t][0])
    
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+1)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air+ 1)*samples)] = np.array(Z["Z_SFH06_{0}_floor".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH06_{0}_floor".format(i + 1)][t][0]) - float(lower["lower_SFH06_{0}_floor".format(i + 1)][t])
            
            """2006-today archetype with underfloor heating system floor temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2 + 2)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air+ 1)*samples)] = np.array(Z_floor["Z_SFH06_{0}_f".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2 + 2)*samples), 0] = float(29) - float(Mmax_floor["Mmax_SFH06_{0}_f".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+ 2 + 1)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = np.array(Z_floor["Z_SFH06_{0}_f".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+ 2 + 1)*samples), 0] = float(Mmax_floor["Mmax_SFH06_{0}_f".format(i+ 1)][t][0]) - float(15)

            """Photovoltaic consumption constraints. """
            Aub_opt[t+((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air + Nc_bui_SS)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air+ n_SFH05_air+ n_SFH90_air + n_SFH60_air+ 1)*samples)] = I[t,:]*((-1)*(1/COP_SFH06_floor[t]))
            Aub_opt[t+((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air + Nc_bui_SS)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air + n_Bui)*samples):((i + n_SFH06_air+ n_SFH05_air+ n_SFH90_air + n_SFH60_air+ n_Bui + 1)*samples)] = I[t,:]*(P_ele_06_floor[t])
            Bub_opt[t+((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air + Nc_bui_SS)*samples), 0] = 0

            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot-Nc_cluster)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air + 1)*samples)] = I[t,:]*(1/COP_SFH06_floor[t])
    
    """Cluster max electricity constraint. """
    for t in range(0, samples):
        Bub_opt[t+((Nc_tot - Nc_cluster)*samples), 0] = Eele_max[t]      # constrained electric consumption at cluster level
    
    """Optimization result (Demand-Response scenario). """
    res = linprog(c=c, A_ub=Aub_opt, b_ub=Bub_opt, bounds=bounds, method='highs-ipm', options={'maxiter': 1000000})
    print("Cluster level:")
    print(res.message + "(Demand-Response scenario)")
    if res.success == False:
        print("\n              ")
        print("\nWARNING")
        print("\nThe thermal demand of the building (decision variable 1) depends on the state of the indoor air temperature Tair.")
        print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
        print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
        print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\nWhereas, the PV utilization factor (decision variable 2) is constrained according to available generation. ")
    
    """Calculation of thermal loads (decision variable x1), PV utilization factos (decision variable x2), 
            indoor air temperatures and building electricity consumption during the DR scenario. """
            
    Qres_dr   = {}  # Values of the decision variable x1 (thermal demand in Wh)
    Tair_dr   = {}  # Indoor air temeperature (°C) - Demand-Response scenario
    Q_dr   = {}     # Thermal demand (kWh) - Demand-Response scenario
    E_dr   = {}     # Electric demand (kWh) - Demand-Response scenario
    
    Upv   = {}      # Values of the decision variable x2 (photovoltaic utilization factor)
    Ele_pv   = {}   # Photovoltaic consumption (kWh) - Demand-Response scenario
    
    for i in range(n_SFH06_air):
        """Values of the decision variables x1 (i.e., thermal load). """
        Qres_dr["Qres_dr_SFH06_{0}".format(i+1)] = (res.x[(i*samples):((i + 1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH06_{0}".format(i+1)] = list(M["M_SFH06_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}".format(i+1)], Qres_dr["Qres_dr_SFH06_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH06_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH06_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption drawn from the grid (DR scenario). """
        E_dr["E_dr_SFH06_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH06_{0}".format(i+1)][r]/COP_SFH06_air[r] for r in range(0, samples)]
        
        """Values of the decision variables x2 (i.e., PV utilization factor). """
        Upv["Upv_SFH06_{0}".format(i+1)]         = (res.x[((i + n_Bui)*samples):((i + n_Bui + 1)*samples)].reshape(samples,1))
        
        """Calculating the PV consumption. """
        Ele_pv["Ele_pv_SFH06_{0}".format(i+1)]   = [Upv["Upv_SFH06_{0}".format(i+1)][r, 0]*P_ele_06_air[r]/1000 for r in range(0, samples)]
    
    for i in range(n_SFH05_air):
        """Values of the decision variables x1 (i.e., thermal load). """
        Qres_dr["Qres_dr_SFH05_{0}".format(i+1)] = (res.x[((i+n_SFH06_air)*samples):((i+n_SFH06_air+1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH05_{0}".format(i+1)] = list(M["M_SFH05_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH05_{0}".format(i+1)], Qres_dr["Qres_dr_SFH05_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH05_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH05_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption drawn from the grid (DR scenario). """
        E_dr["E_dr_SFH05_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH05_{0}".format(i+1)][r]/COP_SFH05_air[r] for r in range(0, samples)]
        
        """Values of the decision variables x2 (i.e., PV utilization factor). """
        Upv["Upv_SFH05_{0}".format(i+1)] = (res.x[((i + n_Bui + n_SFH06_air)*samples):((i + n_Bui + n_SFH06_air + 1)*samples)].reshape(samples,1))
        
        """Calculating the PV consumption. """
        Ele_pv["Ele_pv_SFH05_{0}".format(i+1)]   = [Upv["Upv_SFH05_{0}".format(i+1)][r, 0]*P_ele_05_air[r]/1000 for r in range(0, samples)]
        
    for i in range(n_SFH90_air):
        """Values of the decision variables x1 (i.e., thermal load). """
        Qres_dr["Qres_dr_SFH90_{0}".format(i+1)] = (res.x[((i+n_SFH06_air+n_SFH05_air)*samples):((i+n_SFH06_air+n_SFH05_air+1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH90_{0}".format(i+1)] = list(M["M_SFH90_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH90_{0}".format(i+1)], Qres_dr["Qres_dr_SFH90_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH90_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH90_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption drawn from the grid (DR scenario). """
        E_dr["E_dr_SFH90_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH90_{0}".format(i+1)][r]/COP_SFH90_air[r] for r in range(0, samples)]
        
        """Values of the decision variables x2 (i.e., PV utilization factor). """
        Upv["Upv_SFH90_{0}".format(i+1)] = (res.x[((i + n_Bui + n_SFH06_air + n_SFH05_air)*samples):((i + n_Bui + n_SFH06_air + n_SFH05_air + 1)*samples)].reshape(samples,1))
        
        """Calculating the PV consumption. """
        Ele_pv["Ele_pv_SFH90_{0}".format(i+1)]   = [Upv["Upv_SFH90_{0}".format(i+1)][r, 0]*P_ele_90_air[r]/1000 for r in range(0, samples)]
    
    for i in range(n_SFH60_air):
        """Values of the decision variables x1 (i.e., thermal load). """
        Qres_dr["Qres_dr_SFH60_{0}".format(i+1)] = (res.x[((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH60_{0}".format(i+1)] = list(M["M_SFH60_{0}".format(i+1)][3, :].reshape(samples, 1) + np.dot(Z["Z_SFH60_{0}".format(i+1)], Qres_dr["Qres_dr_SFH60_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH60_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH60_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption drawn from the grid (DR scenario). """
        E_dr["E_dr_SFH60_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH60_{0}".format(i+1)][r]/COP_SFH60_air[r] for r in range(0, samples)]
        
        """Values of the decision variables x2 (i.e., PV utilization factor). """
        Upv["Upv_SFH60_{0}".format(i+1)] = (res.x[((i + n_Bui + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_Bui + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)].reshape(samples,1))
        
        """Calculating the PV consumption. """
        Ele_pv["Ele_pv_SFH60_{0}".format(i+1)]   = [Upv["Upv_SFH60_{0}".format(i+1)][r, 0]*P_ele_60_air[r]/1000 for r in range(0, samples)]

    for i in range(n_SFH06_floor):
        """Values of the decision variables x1 (i.e., thermal load). """
        Qres_dr["Qres_dr_SFH06_{0}_floor".format(i+1)] = (res.x[((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air + 1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH06_{0}_floor".format(i+1)] = list(M["M_SFH06_{0}_floor".format(i+1)][7, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}_floor".format(i+1)], Qres_dr["Qres_dr_SFH06_{0}_floor".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH06_{0}_floor".format(i+1)]       = [Qres_dr["Qres_dr_SFH06_{0}_floor".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption drawn from the grid (DR scenario). """
        E_dr["E_dr_SFH06_{0}_floor".format(i+1)]       = [Q_dr["Q_dr_SFH06_{0}_floor".format(i+1)][r]/COP_SFH06_floor[r] for r in range(0, samples)]

        Upv["Upv_SFH06_{0}_floor".format(i+1)] = (res.x[((i + n_Bui + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*samples):((i + n_Bui + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air + 1)*samples)].reshape(samples,1))
        
        """Calculating the PV consumption. """
        Ele_pv["Ele_pv_SFH06_{0}_floor".format(i+1)]   = [Upv["Upv_SFH06_{0}_floor".format(i+1)][r, 0]*P_ele_06_floor[r]/1000 for r in range(0, samples)]

    """Conversion of dictionary to data frame for calculation of aggregate electric power consumption. """
    E_dr = pd.DataFrame.from_dict(E_dr)
    
    """Calculating electric consumption at the aggregate level during DR scenario. """
    Etot_dr = E_dr.sum(axis=1)
    Edr_sum = float(Etot_dr.sum())
    Etot_dr = [Etot_dr[i] for i in range(0,samples)]  
    
    """Conversion of dictionary to data frame for calculation of aggregate PV consumption. """
    E_pv = pd.DataFrame.from_dict(Ele_pv)
    
    """Calculating PV consumption at the aggregate level. """
    Etot_pv = E_pv.sum(axis=1)
    Epv_dr_sum = float(Etot_pv.sum())
    Etot_pv = [Etot_pv[i] for i in range(0,samples)]  
    
    PVprod =[(P_ele_06_air[r]*n_SFH06_air + P_ele_05_air[r]*n_SFH05_air + P_ele_90_air[r]*n_SFH90_air + P_ele_60_air[r]*n_SFH60_air + P_ele_06_floor[r]*n_SFH06_floor)/1000 for r in range(0,samples)]
    
    """Calculation of solar PV consumed during BL scenario useful for comparison. """
    PVdiff_bl   = []
    Etot_pv_BL = []
    for t in range(samples):
        PVdiff_bl.append(PVprod[t]-Etot_bl[t])
        if PVprod[t] > 0:
            if PVdiff_bl[t] > 0: # PV surplus
                Etot_pv_BL.append(Etot_bl[t])
            else:
                Etot_pv_BL.append(PVprod[t])
        else:
            Etot_pv_BL.append(0)
            
    PVprod = np.array(PVprod)
    PVprod_sum = float(PVprod.sum())
    
    Etot_pv_BL = np.array(Etot_pv_BL)
    Epv_bl_sum = float(Etot_pv_BL.sum())
    
    import matplotlib.pyplot as plt
    time = [i*timestep for i in range(0, samples)] # time array
    
    """Do plot. """
    if do_plot_E_dr_cluster:
        """If True, plot cluster electric and PV consumption trends. """
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Time (hours)', fontsize=14)
        ax1.set_ylabel('$\.P_\mathrm{cluster}$ (kW$_\mathrm{el}$)', fontsize=14)
        ax1.plot(time, Etot_bl, 'k', label="BL_scenario", linewidth=1.0)
        ax1.plot(time, Etot_dr, 'r', label="DR_scenario", linewidth=1.0)
        ax1.plot(time, PVprod, 'g--', label="PV_gen", linewidth=1.0)
        ax1.set_xlim(0, duration-1)
        ax1.grid(alpha=0.3)
        plt.legend(frameon=False, fontsize=12)
        fig.set_size_inches(w=10, h=4)
        fig.tight_layout()
        plt.savefig("Eele_cluster", dpi=200)
    
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Time (hours)', fontsize=14)
        ax1.set_ylabel('$\.P_\mathrm{PV cluster}$ (kW$_\mathrm{el}$)', fontsize=14)
        ax1.plot(time, Etot_pv_BL, 'k', label="PV_cons_BL", linewidth=1.0)
        ax1.plot(time, Etot_pv, 'r', label="PV_cons_DR", linewidth=1.0)
        ax1.plot(time, PVprod, 'g--', label="PV_gen", linewidth=1.0)
        ax1.set_xlim(0, duration-1)
        ax1.grid(alpha=0.3)
        plt.legend(frameon=False, fontsize=12)
        fig.set_size_inches(w=10, h=4)
        fig.tight_layout()
        plt.savefig("pv_cluster", dpi=200)
        
        plt.show()
        plt.close('all')

    if do_plot_E_dr_single:
        """If True, plot single buildings electric and PV consumption trends. """
        for i in range (n_SFH06_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            ax1.plot(time,PVprod,'g--', linewidth=1, label = "PV_gen")
            ax1.plot(time,E_dr["E_dr_SFH06_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH06_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH06_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH06_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
            
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH05_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            ax1.plot(time,PVprod,'g--', linewidth=1, label = "PV_gen")
            ax1.plot(time,E_dr["E_dr_SFH05_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH05_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH05_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH05_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH90_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            ax1.plot(time,PVprod,'g--', linewidth=1, label = "PV_gen")
            ax1.plot(time,E_dr["E_dr_SFH90_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH90_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH90_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH90_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH60_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            ax1.plot(time,PVprod,'g--', linewidth=1, label = "PV_gen")
            ax1.plot(time,E_dr["E_dr_SFH60_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH60_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH60_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH60_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
    
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH06_floor):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            ax1.plot(time,PVprod,'g--', linewidth=1, label = "PV_gen")
            ax1.plot(time,E_dr["E_dr_SFH06_{0}_floor".format(i+1)],'r', linewidth=1, label = "E_dr_SFH06_{0}_floor".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH06_{0}_floor".format(i+1)],'k', linewidth=1, label = "E_bl_SFH06_{0}_floor".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
            
        plt.show()
        plt.close('all')

    if do_plot_Tair_dr_single:
        """If True, plot single buildings indoor air temperature trends. """
        for i in range (n_SFH06_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            ax1.plot(time,Tair_dr["Tair_dr_SFH06_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH06_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH06_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH06_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
            
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH05_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            ax1.plot(time,Tair_dr["Tair_dr_SFH05_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH05_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH05_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH05_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH90_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            ax1.plot(time,Tair_dr["Tair_dr_SFH90_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH90_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH90_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH90_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH60_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            ax1.plot(time,Tair_dr["Tair_dr_SFH60_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH60_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH60_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH60_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
    
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH06_floor):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            ax1.plot(time,Tair_dr["Tair_dr_SFH06_{0}_floor".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH06_{0}_floor".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH06_{0}_floor".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH06_{0}_floor".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
    """Data export to excel. """
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writerEele = pd.ExcelWriter(name_file + '.xlsx', engine='xlsxwriter')

    # Write each dataframe to a different worksheet.
    E_bl.to_excel(writerEele, sheet_name='Eele_BL')
    E_dr.to_excel(writerEele, sheet_name='Eele_DR')
    E_pv.to_excel(writerEele, sheet_name='Epv_cons_DR')

    # Close the Pandas Excel writer and output the Excel file.
    writerEele.close()
    
    """Electricity consumption during BL and DR scenarios. """
    print("\nCluster electricity consumption of about " +
          str(round(Ebl_sum,2)) + " kWh during baseline (BL) scenario")
    print("\n------------------------------------------------------------------------")
    
    Diff   = ((Edr_sum - Ebl_sum)/Ebl_sum)*100
    
    if Diff >= 0:
        print("\nCluster electricity consumption of about " +
              str(round(Edr_sum,2)) + " kWh during DR scenario with centralized PV generation,"
              "\nwith a " + str(round(Diff,2)) + "% increase compared to BL scenario")
        print("\n------------------------------------------------------------------------")
    else:
        print("\nCluster electricity consumption of about " +
              str(round(Edr_sum,2)) + " kWh during DR scenario with centralized PV generation,"
              "\nwith a " + str(round(Diff*(-1),2)) + "% reduction compared to BL scenario")
        print("\n------------------------------------------------------------------------")
    
    """PV consumption during BL and DR scenarios. """
    print("\nCluster PV consumption of about " +
          str(round(Epv_bl_sum,2)) + " kWh during baseline (BL) scenario" + 
          "\n(consumption of " + str(round(Epv_bl_sum/PVprod_sum*100,2)) + "%)")
    print("\n------------------------------------------------------------------------")
   
    Diff_pv = ((Epv_dr_sum-Epv_bl_sum)/Epv_bl_sum)*100
    
    if Diff_pv >= 0:
        print("\nCluster PV consumption of about " +
              str(round(Epv_dr_sum,2)) + " kWh during DR scenario with centralized PV generation" + 
              "\n(consumption of " + str(round(Epv_dr_sum/PVprod_sum*100,2)) + "%)," 
              "\nwith a " + str(round(Diff_pv,2)) + "% increase compared to BL scenario")
        print("\n------------------------------------------------------------------------")
    else:
        print("\nCluster PV consumption of about " +
              str(round(Epv_dr_sum,2)) + " kWh during DR scenario with centralized PV generation" + 
              "\n(consumption of " + str(round(Epv_dr_sum/PVprod_sum*100,2)) + "%)," 
              "\nwith a " + str(round(Diff_pv,2)) + "% decrease compared to BL scenario")
        print("\n------------------------------------------------------------------------")
    
    return(Etot_dr)

def load_shifting(Locality, day, month, duration, timestep, Thermal_load, 
                            n_SFH06_air, n_SFH05_air, n_SFH90_air, n_SFH60_air,n_SFH06_floor, 
                            HP_data_type, 
                            Supply_air, rated_cap_06_air, rated_cap_05_air, rated_cap_90_air, rated_cap_60_air, rated_COP_06_air, rated_COP_05_air, rated_COP_90_air, rated_COP_60_air,
                            Supply_floor, rated_cap_06_floor, rated_COP_06_floor,
                            Setpoint_type, setpoint, th_tolerance_BL_sup_air, th_tolerance_BL_min_air, th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor,
                            th_tolerance_DR_sup_air, th_tolerance_DR_min_air, th_tolerance_DR_sup_floor, th_tolerance_DR_min_floor,
                            price_ranges, f_limit,
                            do_plot_E_dr_single, do_plot_E_dr_cluster,do_plot_Tair_dr_single, name_file):
    """load_shifting is a function based on linear programming that minimizes aggregate electricity consumption 
            during the application of a strategy of load shifting under a bi-hourly price signal. """
    
    # number of buildings:
    n_SFH06_air = int(n_SFH06_air)
    n_SFH05_air = int(n_SFH05_air)
    n_SFH90_air = int(n_SFH90_air)
    n_SFH60_air = int(n_SFH60_air)
    
    if Thermal_load == "heating":
        n_SFH06_floor = int(n_SFH06_floor) # Radiant floor emission system for the SFH06 archetype (age class 2006-today) only during heating
    else: 
        n_SFH06_floor = 0
    
    n_Bui       = n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air + n_SFH06_floor # Total number of buildings composing the cluster
    
    samples = int(duration/timestep) # Calculation of total steps

    T_est = MeteoFile(Locality, day, month, duration, timestep)[0] # Outdoor air temperature (°C)
    
    """Call of the performance characteristics of heat pumps based on the definition of heat loads 
           and attribute definition for calling temperature setpoint functions"""
    if Thermal_load == 'cooling':
        """Call performance characteristics and capacity of heat pumps. """
        if n_SFH06_air > 0:
            COP_SFH06_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH06_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH05_air > 0:
            COP_SFH05_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH05_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH90_air > 0:
            COP_SFH90_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH90_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH60_air > 0:
            COP_SFH60_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH60_air = [HPperformance.HP_cooling(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[1] for i in range(0,samples)]
            
        """Thermostat mode setting in space cooling"""
        thermostat = 'th_cooling'
    elif Thermal_load == 'heating':
        if n_SFH06_air > 0:
            COP_SFH06_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH06_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_06_air,rated_COP_06_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH05_air > 0:
            COP_SFH05_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH05_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_05_air,rated_COP_05_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH90_air > 0:
            COP_SFH90_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH90_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_90_air,rated_COP_90_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH60_air > 0:
            COP_SFH60_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH60_air = [HPperformance.HP_heating(T_est[i],Supply_air,rated_cap_60_air,rated_COP_60_air,HP_data_type)[1] for i in range(0,samples)]
        if n_SFH06_floor > 0:
            COP_SFH06_floor = [HPperformance.HP_heating(T_est[i],Supply_floor,rated_cap_06_air,rated_COP_06_air,HP_data_type)[0] for i in range(0,samples)]
            Q_SFH06_floor = [HPperformance.HP_heating(T_est[i],Supply_floor,rated_cap_06_air,rated_COP_06_air,HP_data_type)[1] for i in range(0,samples)]
        
        """Thermostat mode setting in space heating"""
        thermostat = 'th_heating'
    else:
        print('Thermal load not available')
    
    """Definition of occupancy profiles, relative internal gains (people, equipment and lighting) and setpoint temperature for each building 
           (useful for diversification of consumption loads within the cluster)"""
    occ   = {} # occupancy
    gains = {} # internal gains (W)
    tsp   = {} # thermostat set-points (°C)
    
    if Setpoint_type == "pre_defined":
        """Pre-defined temperature set-point definition
               For adding profiles or modifying existing ones, go to the User_pattern.py file and see the functions 
               th_heating (set-points for space heating) or th_cooling (set-points for space cooling)"""
        for i in range(n_SFH06_air):
            tsp["tspSFH06_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i]
        for i in range(n_SFH05_air):
            tsp["tspSFH05_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air]
        for i in range(n_SFH90_air):
            tsp["tspSFH90_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air]
        for i in range(n_SFH60_air):
            tsp["tspSFH60_{0}".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air]
        for i in range(n_SFH06_floor):
            tsp["tspSFH06_{0}_floor".format(i+1)] = getattr(User_pattern,thermostat)(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air]
    
    elif Setpoint_type == "user_defined":
        """User-defined temperature set-point definition"""
        t_setpoint = []
        for i in range (samples):
            t_setpoint.append(setpoint)
        for i in range(n_SFH06_air):
            tsp["tspSFH06_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH05_air):
            tsp["tspSFH05_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH90_air):
            tsp["tspSFH90_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH60_air):
            tsp["tspSFH60_{0}".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
        for i in range(n_SFH06_floor):
            tsp["tspSFH06_{0}_floor".format(i+1)] = [t_setpoint[r] for r in range(0, samples)]
    
    """Occupancy patterns and related internal gain profiles definition"""
    for i in range(n_SFH06_air):
        occ["occupancySFH06_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i]
        gains["gainsSFH06_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i]
    for i in range(n_SFH05_air):
        occ["occupancySFH05_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air]
        gains["gainsSFH05_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air]
    for i in range(n_SFH90_air):
        occ["occupancySFH90_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air]
        gains["gainsSFH90_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air]
    for i in range(n_SFH60_air):
        occ["occupancySFH60_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air]
        gains["gainsSFH60_{0}".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air]
    for i in range(n_SFH06_floor):
        occ["occupancySFH06_{0}_floor".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air]
        gains["gainsSFH06_{0}_floor".format(i+1)] = intGains(n_Bui,day,month,duration,timestep)[i+n_SFH06_air+n_SFH05_air+n_SFH90_air+n_SFH60_air]
        
    """Defining the price signal. """
    time = np.arange(0,24,timestep)  # hour of the day
    n_days        = int(duration/24) # number of simulation days
    hour = []
    for n in range(n_days):
        for t in range(len(time)):
            hour.append(time[t])

    ele_grid_price = []   # price signal array

    for t in range(0,samples):
        if hour[t] >= price_ranges["hour_day"] and hour[t] < price_ranges["hour_night"]:
            ele_grid_price.append(price_ranges["price_day"]) # Eur/kWh
        else:
            ele_grid_price.append(price_ranges["price_night"]) # Eur/kWh
            
    """In the following, thermal loads provided as input to the state-space model formulation are calculated via linear programming . 
           Specifically, the thermal load of each building is minimized to maintain the indoor temperature set-point.
           Linear programming solves problems of the following form: 
               minimize c @ x
               such that:
               A_ub @ x <= b_ub
               A_eq @ x == b_eq
               lb <= x <= ub
               
               For more information see the documentation of the python library scipy.optimize.linprog"""
           
    c      = {}     # Coefficients of the linear objective function to be minimized
    Aub    = {}     # Matrix of inequality constraints
    Bub    = {}     # Vector of inequality constraints
    bounds = {}     # Sequence of pairs (min, max) for each element in x, defining the minimum and maximum values of the decision variable
    
    Z      = {}     # Variable describing the state of the system containing the decision variable (e.g., thermal demand)
    M      = {}     # Variable describing the state of the indoor air temperture
    Mmax      = {}  # Variable describing the state of the indoor air temperture
    res_opt   = {}  # Optimization result - baseline scenario
    Qres_bl   = {}  # Values of the decision variable (thermal demand in Wh)
    Tair_bl   = {}  # Indoor air temeperature (°C) - baseline scenario
    Q_bl   = {}     # Thermal demand (kWh) - baseline scenario
    E_bl   = {}     # Electric demand (kWh) - baseline scenario
    upper   = {}    # Upper range thermostat set-point profile (°C) - baseline scenario
    lower   = {}    # Lower range thermostat set-point profile (°C) - baseline scenario
    cost_bl = {}    # Electricity bill (Eur) - baseline scenario
    
    Mmax_floor = {} # Variable describing the state of the floor temperture
    Z_floor    = {} # Variable describing the state of the floor system containing the decision variable (e.g., thermal demand)

    """Calling archetype (SFHs) functions from the "Archetypes.py" module. """
    for i in range(n_SFH06_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH06 buildings (age class 2006-today) with air emission system."""
        c["c_SFH06_{0}".format(i+1)]             = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[0]
        Aub["Aub_SFH06_{0}".format(i+1)]         = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[1]
        Bub["Bub_SFH06_{0}".format(i+1)]         = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[2]
        bounds["bounds_SFH06_{0}".format(i+1)]   = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH06_{0}".format(i+1)]             = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[4]
        Z["Z_SFH06_{0}".format(i+1)]             = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[5]
        Mmax["Mmax_SFH06_{0}".format(i+1)]       = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH06_{0}".format(i+1)]     = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[7]
        lower["lower_SFH06_{0}".format(i+1)]     = SFH2006_air(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH06_{0}".format(i+1)], gains["gainsSFH06_{0}".format(i+1)],Q_SFH06_air)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH06_{0}".format(i+1)] = linprog(c=c["c_SFH06_{0}".format(i+1)], A_ub=Aub["Aub_SFH06_{0}".format(i+1)], b_ub=Bub["Bub_SFH06_{0}".format(i+1)],bounds=bounds["bounds_SFH06_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 2006-today building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH06_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH06_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH06_{0}".format(i+1)] = (res_opt["res_opt_SFH06_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH06_{0}".format(i+1)] = list(M["M_SFH06_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}".format(i+1)], Qres_bl["Qres_bl_SFH06_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH06_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH06_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH06_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH06_{0}".format(i+1)][r]/COP_SFH06_air[r] for r in range(0, samples)]
        
        """Calculating the cost of electiricty drawn from the grid (baseline scenario). """
        cost_bl["cost_bl_SFH06_{0}".format(i+1)] = [E_bl["E_bl_SFH06_{0}".format(i+1)][r]*ele_grid_price[r] for r in range(0, samples)]
        
    for i in range(n_SFH05_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH05 buildings (age class 1991-2005) with air emission system."""
        c["c_SFH05_{0}".format(i+1)]             = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[0]
        Aub["Aub_SFH05_{0}".format(i+1)]         = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[1]
        Bub["Bub_SFH05_{0}".format(i+1)]         = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[2]
        bounds["bounds_SFH05_{0}".format(i+1)]   = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH05_{0}".format(i+1)]             = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[4]
        Z["Z_SFH05_{0}".format(i+1)]             = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[5]
        Mmax["Mmax_SFH05_{0}".format(i+1)]       = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH05_{0}".format(i+1)] = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[7]
        lower["lower_SFH05_{0}".format(i+1)] = SFH2005_air(Locality, day, month, duration, timestep, tsp["tspSFH05_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH05_{0}".format(i+1)], gains["gainsSFH05_{0}".format(i+1)],Q_SFH05_air)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH05_{0}".format(i+1)] = linprog(c=c["c_SFH05_{0}".format(i+1)], A_ub=Aub["Aub_SFH05_{0}".format(i+1)], b_ub=Bub["Bub_SFH05_{0}".format(i+1)],bounds=bounds["bounds_SFH05_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 1990-2005 building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH05_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH05_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH05_{0}".format(i+1)] = (res_opt["res_opt_SFH05_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH05_{0}".format(i+1)] = list(M["M_SFH05_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH05_{0}".format(i+1)], Qres_bl["Qres_bl_SFH05_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH05_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH05_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH05_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH05_{0}".format(i+1)][r]/COP_SFH05_air[r] for r in range(0, samples)]
        
        """Calculating the cost of electiricty drawn from the grid (baseline scenario). """
        cost_bl["cost_bl_SFH05_{0}".format(i+1)] = [E_bl["E_bl_SFH05_{0}".format(i+1)][r]*ele_grid_price[r] for r in range(0, samples)]
        
    for i in range(n_SFH90_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH90 buildings (age class 1976-1990) with air emission system."""
        c["c_SFH90_{0}".format(i+1)]             = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[0]
        Aub["Aub_SFH90_{0}".format(i+1)]         = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[1]
        Bub["Bub_SFH90_{0}".format(i+1)]         = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[2]
        bounds["bounds_SFH90_{0}".format(i+1)]   = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH90_{0}".format(i+1)]             = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[4]
        Z["Z_SFH90_{0}".format(i+1)]             = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[5]
        Mmax["Mmax_SFH90_{0}".format(i+1)]       = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH90_{0}".format(i+1)]     = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[7]
        lower["lower_SFH90_{0}".format(i+1)]     = SFH1990_air(Locality, day, month, duration, timestep, tsp["tspSFH90_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH90_{0}".format(i+1)], gains["gainsSFH90_{0}".format(i+1)],Q_SFH90_air)[8]
        
        """Optimization result (baseline scenario)"""
        res_opt["res_opt_SFH90_{0}".format(i+1)] = linprog(c=c["c_SFH90_{0}".format(i+1)], A_ub=Aub["Aub_SFH90_{0}".format(i+1)], b_ub=Bub["Bub_SFH90_{0}".format(i+1)],bounds=bounds["bounds_SFH90_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 1976-1990 building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH90_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH90_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH90_{0}".format(i+1)] = (res_opt["res_opt_SFH90_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH90_{0}".format(i+1)] = list(M["M_SFH90_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH90_{0}".format(i+1)], Qres_bl["Qres_bl_SFH90_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH90_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH90_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH90_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH90_{0}".format(i+1)][r]/COP_SFH90_air[r] for r in range(0, samples)]
        
        """Calculating the cost of electiricty drawn from the grid (baseline scenario). """
        cost_bl["cost_bl_SFH90_{0}".format(i+1)] = [E_bl["E_bl_SFH90_{0}".format(i+1)][r]*ele_grid_price[r] for r in range(0, samples)]
        
    for i in range(n_SFH60_air):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH60 buildings (age class 1946-1960) with air emission system."""
        c["c_SFH60_{0}".format(i+1)]             = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[0]
        Aub["Aub_SFH60_{0}".format(i+1)]         = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[1]
        Bub["Bub_SFH60_{0}".format(i+1)]         = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[2]
        bounds["bounds_SFH60_{0}".format(i+1)]   = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[3]
        
        """Calling of variables referring to the state of indoor air temperature. """
        M["M_SFH60_{0}".format(i+1)]             = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[4]
        Z["Z_SFH60_{0}".format(i+1)]             = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[5]
        Mmax["Mmax_SFH60_{0}".format(i+1)]       = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[6]
        
        """Calling upper and lower temperature profiles that define the temperature range set in the tersmostat during the BL scenario. """
        upper["upper_SFH60_{0}".format(i+1)] = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[7]
        lower["lower_SFH60_{0}".format(i+1)] = SFH1960_air(Locality, day, month, duration, timestep, tsp["tspSFH60_{0}".format(i+1)], th_tolerance_BL_sup_air, th_tolerance_BL_min_air, Thermal_load, occ["occupancySFH60_{0}".format(i+1)], gains["gainsSFH60_{0}".format(i+1)],Q_SFH60_air)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH60_{0}".format(i+1)] = linprog(c=c["c_SFH60_{0}".format(i+1)], A_ub=Aub["Aub_SFH60_{0}".format(i+1)], b_ub=Bub["Bub_SFH60_{0}".format(i+1)],bounds=bounds["bounds_SFH60_{0}".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype 1946-1960 building {0}".format(i+1) + " with air " + Thermal_load + " system:")
        print(res_opt["res_opt_SFH60_{0}".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH60_{0}".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")

        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH60_{0}".format(i+1)] = (res_opt["res_opt_SFH60_{0}".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH60_{0}".format(i+1)] = list(M["M_SFH60_{0}".format(i+1)][3, :].reshape(samples, 1) + np.dot(Z["Z_SFH60_{0}".format(i+1)], Qres_bl["Qres_bl_SFH60_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH60_{0}".format(i+1)]       = [Qres_bl["Qres_bl_SFH60_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH60_{0}".format(i+1)]       = [Q_bl["Q_bl_SFH60_{0}".format(i+1)][r]/COP_SFH60_air[r] for r in range(0, samples)]
        
        """Calculating the cost of electiricty drawn from the grid (baseline scenario). """
        cost_bl["cost_bl_SFH60_{0}".format(i+1)] = [E_bl["E_bl_SFH60_{0}".format(i+1)][r]*ele_grid_price[r] for r in range(0, samples)]
        
    for i in range(n_SFH06_floor):
        """Calling the coefficients of the linear objective function (c), matrix of inequality constraints (Aub), 
        vector of inequality constraints (Bub), sequence of pairs defining the minimum and maximum values of the decision variable 
        (bounds) for calculating the thermal demand of SFH06 buildings (age class 2006-today) with underfloor heating system."""
        c["c_SFH06_{0}_floor".format(i+1)]             = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[0]
        Aub["Aub_SFH06_{0}_floor".format(i+1)]         = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[1]
        Bub["Bub_SFH06_{0}_floor".format(i+1)]         = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[2]
        bounds["bounds_SFH06_{0}_floor".format(i+1)]   = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[3]
        
        """Calling of variables referring to the state of indoor air temperature node. """
        M["M_SFH06_{0}_floor".format(i+1)]             = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[4]
        Z["Z_SFH06_{0}_floor".format(i+1)]             = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[5]
        Mmax["Mmax_SFH06_{0}_floor".format(i+1)]       = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[6]
        
        """Calling of variables referring to the state of floor node. """
        Z_floor["Z_SFH06_{0}_f".format(i+1)]           = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[9]
        Mmax_floor["Mmax_SFH06_{0}_f".format(i+1)]     = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[10]
        
        """Calling of variables referring to the state of indoor air temperature. """
        upper["upper_SFH06_{0}_floor".format(i+1)]     = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[7]
        lower["lower_SFH06_{0}_floor".format(i+1)]     = SFH2006_floor(Locality, day, month, duration, timestep, tsp["tspSFH06_{0}_floor".format(i+1)], th_tolerance_BL_sup_floor, th_tolerance_BL_min_floor, occ["occupancySFH06_{0}_floor".format(i+1)], gains["gainsSFH06_{0}_floor".format(i+1)],Q_SFH06_floor)[8]
        
        """Optimization result (baseline scenario). """
        res_opt["res_opt_SFH06_{0}_floor".format(i+1)] = linprog(c=c["c_SFH06_{0}_floor".format(i+1)], A_ub=Aub["Aub_SFH06_{0}_floor".format(i+1)], b_ub=Bub["Bub_SFH06_{0}_floor".format(i+1)],bounds=bounds["bounds_SFH06_{0}_floor".format(i+1)], method='highs-ipm', options={'maxiter': 1000000})
        print("Archetype building {0}".format(i+1) + " with underfloor heating system:")
        print(res_opt["res_opt_SFH06_{0}_floor".format(i+1)].message + "(baseline scenario)")
        if res_opt["res_opt_SFH06_{0}_floor".format(i+1)].success == False:
            print("\n              ")
            print("\nWARNING")
            print("\nThe thermal demand of the building (decision variable) depends on the state of the indoor air temperature Tair.")
            print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
            print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
            print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
        print("\n------------------------------------------------------------------------")
        
        """Base load - Values of the decision variables that minimize the objective function and satisfy the constraints (i.e., thermal power). """
        Qres_bl["Qres_bl_SFH06_{0}_floor".format(i+1)] = (res_opt["res_opt_SFH06_{0}_floor".format(i+1)].x[0:(samples)].reshape(samples, 1))
        
        """Indoor air temperature calculation (baseline scenario). """
        Tair_bl["Tair_bl_SFH06_{0}_floor".format(i+1)] = list(M["M_SFH06_{0}_floor".format(i+1)][7, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}_floor".format(i+1)], Qres_bl["Qres_bl_SFH06_{0}_floor".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_bl["Q_bl_SFH06_{0}_floor".format(i+1)]       = [Qres_bl["Qres_bl_SFH06_{0}_floor".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption (baseline scenario) to be managed during Demand-Response scenario. """
        E_bl["E_bl_SFH06_{0}_floor".format(i+1)]       = [Q_bl["Q_bl_SFH06_{0}_floor".format(i+1)][r]/COP_SFH06_floor[r] for r in range(0, samples)]
        
        """Calculating the cost of electiricty drawn from the grid (baseline scenario). """
        cost_bl["cost_bl_SFH06_{0}_floor".format(i+1)] = [E_bl["E_bl_SFH06_{0}_floor".format(i+1)][r]*ele_grid_price[r] for r in range(0, samples)]
        
    """Conversion of dictionary to data frame for calculation of aggregate electric power consumption. """
    E_bl = pd.DataFrame.from_dict(E_bl)
    
    """Calculating base electrical loads at the aggregate level. """
    Etot_bl = E_bl.sum(axis=1)
    Ebl_sum = float(Etot_bl.sum())
    Etot_bl = [Etot_bl[i] for i in range(0,samples)]  
    
    """Conversion of dictionary to data frame for calculation of aggregate electricity bill. """
    cost_bl = pd.DataFrame.from_dict(cost_bl)
    
    """Calculating base electricity bill at the aggregate level. """
    cost_tot_bl = cost_bl.sum(axis=1)
    cost_bl_sum = float(cost_tot_bl.sum())
    cost_tot_bl = [cost_tot_bl[i] for i in range(0,samples)]  
    
    """------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    The Demand-Response (DR) event is modeled below. 
    From the definition of electricity price profile and electricity demand during the baseline scenario, electricity demand 
    during the DR scenario is managed by shifting building loads toward the hours when electricity is cheaper.
    Specifically, electricity for each building drawn from the grid is minimized:
        min (Q_dr*(ele_grid_price/COP)), 
        with
        Q_dr = thermal load (decision variable).
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""
    """For each simulation day, calculation of the aggregate peak load. """

    Max         = []
    Max_dr      = []

    for n in range (n_days):
        Max.append(np.amax(np.array(Etot_bl[int(samples/n_days)*n:int(samples/n_days)*(n + 1)])))
        for i in range (int(samples/n_days)*n,int(samples/n_days)*(n + 1)):
            Max_dr.append(Max[n]) # maximum electricity consumption during the baseline scenario

    """Definition of the electricity demand of the cluster of buildings constrained to remain below a certain limit to reduce rebound effects.. """
    Eele_max = []

    for t in range(samples):
            Eele_max.append(float(Max_dr[t]*1000)*(float(f_limit)))
            
    """Definition of temperature profiles during the DR event to unlock building energy flexibility. """
    if Thermal_load == 'cooling':
        for i in range(n_SFH06_air):
            for t in range(samples):
                if hour[t] >= price_ranges["hour_day"] and hour[t] < price_ranges["hour_night"]:
                    upper["upper_SFH06_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for cooling to reduce electric consumption from the grid
                else:
                    lower["lower_SFH06_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling
        for i in range(n_SFH05_air):
            for t in range(samples):
                if hour[t] >= price_ranges["hour_day"] and hour[t] < price_ranges["hour_night"]:
                    upper["upper_SFH05_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for cooling to reduce electric consumption from the grid
                else:
                    lower["lower_SFH05_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling
        for i in range(n_SFH90_air):
            for t in range(samples):
                if hour[t] >= price_ranges["hour_day"] and hour[t] < price_ranges["hour_night"]:
                    upper["upper_SFH90_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for cooling to reduce electric consumption from the grid
                else:
                    lower["lower_SFH90_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling
        for i in range(n_SFH60_air):
            for t in range(samples):
                if hour[t] >= price_ranges["hour_day"] and hour[t] < price_ranges["hour_night"]:
                    upper["upper_SFH60_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for cooling to reduce electric consumption from the grid
                else:
                    lower["lower_SFH60_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling

    elif Thermal_load == 'heating':
        for i in range(n_SFH06_air):
            for t in range(samples):
                if hour[t] >= price_ranges["hour_day"] and hour[t] < price_ranges["hour_night"]:
                    lower["lower_SFH06_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to reduce electric consumption from the grid
                else:
                    upper["upper_SFH06_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for heating 
        for i in range(n_SFH05_air):
            for t in range(samples):
                if hour[t] >= price_ranges["hour_day"] and hour[t] < price_ranges["hour_night"]:
                    lower["lower_SFH05_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to reduce electric consumption from the grid
                else:
                    upper["upper_SFH05_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for heating 
        for i in range(n_SFH90_air):
            for t in range(samples):
                if hour[t] >= price_ranges["hour_day"] and hour[t] < price_ranges["hour_night"]:
                    lower["lower_SFH90_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to reduce electric consumption from the grid
                else:
                    upper["upper_SFH90_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for heating 
        for i in range(n_SFH60_air):
            for t in range(samples):
                if hour[t] >= price_ranges["hour_day"] and hour[t] < price_ranges["hour_night"]:
                    lower["lower_SFH60_{0}".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to reduce electric consumption from the grid
                else:
                    upper["upper_SFH60_{0}".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for heating 
        for i in range(n_SFH06_floor):
            for t in range(samples):
                if hour[t] >= price_ranges["hour_day"] and hour[t] < price_ranges["hour_night"]:
                    lower["lower_SFH06_{0}_floor".format(i+1)][t] -= float(th_tolerance_DR_min_air) # Reduction of indoor air temperature set point for cooling to reduce electric consumption from the grid
                else:
                    upper["upper_SFH06_{0}_floor".format(i+1)][t] += float(th_tolerance_DR_sup_air) # Increase in indoor air temperature set point for heating 
    
    """In the following, thermal loads provided as input to the state-space model formulation are calculated via linear programming. 
           Specifically, the thermal load of each building is minimized to maintain the indoor temperature set-point, 
           and the electricity demand of the cluster of buildings is constrained to remain below a certain 
           percentage reduction (Eele_max) for the duration of the peak-shaving event (duration_dr).
           Linear programming solves problems of the following form: 
               minimize c @ x
               such that:
               A_ub @ x <= b_ub
               A_eq @ x == b_eq
               lb <= x <= ub
               
               For more information see the documentation of the python library scipy.optimize.linprog"""
    
    """Definition of the number of constraints at the individual building and cluster level. """
    Nc_06_air   = (2)*n_SFH06_air          # 2 for indoor air temperature range 
    Nc_05_air   = (2)*n_SFH05_air          # 2 for indoor air temperature range 
    Nc_90_air   = (2)*n_SFH90_air          # 2 for indoor air temperature range 
    Nc_60_air   = (2)*n_SFH60_air          # 2 for indoor air temperature range 
    Nc_06_floor = (4)*n_SFH06_floor        # 2 for indoor air temperature range + 2 for temperature range of underfloor heating system
    Nc_bui     = Nc_06_air + Nc_05_air + Nc_90_air + Nc_60_air + Nc_06_floor
    
    Nc_cluster = 1                         # Electric demand at cluster level
    Nc_tot     = int(Nc_bui + Nc_cluster)  # Total number of constrains
    
    c      = [] # Coefficients of the linear objective function to be minimized (electric consumption)
    bounds = [] # Sequence of pairs (min, max) for each element in x, defining the minimum and maximum values of the decision variable
    
    """Setting the bounds and coefficients of 2006-today archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH06_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH06_air[r]))            # bounds of the decision variable (thermal power, Qth)
            c.append(float(ele_grid_price[r]/COP_SFH06_air[r]))          # having x = Qth and c = 1/COP minimizes the power consumption
    
    """Setting the bounds and coefficients of 1991-2005 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH05_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH05_air[r]))
            c.append(float(ele_grid_price[r]/COP_SFH05_air[r]))          # having x = Qth and c = 1/COP minimizes the electric consumption 
    
    """Setting the bounds and coefficients of 1976-1990 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH90_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH90_air[r]))
            c.append(float(ele_grid_price[r]/COP_SFH90_air[r]))          # having x = Qth and c = 1/COP minimizes the electric consumption 
    
    """Setting the bounds and coefficients of 1946-1960 archetypes with air cooling/heating systems. """
    for i in range(0,n_SFH60_air):
        for r in range(0,samples):
            bounds.append((0,Q_SFH60_air[r]))
            c.append(float(ele_grid_price[r]/COP_SFH60_air[r]))          # having x = Qth and c = 1/COP minimizes the electric consumption 
    
    """Setting the bounds and coefficients of 2006-today archetypes with underfloor heating system. """
    for i in range(0,n_SFH06_floor):
        for r in range(0,samples):
            bounds.append((0,Q_SFH06_floor[r])) 
            c.append(float(ele_grid_price[r]/COP_SFH06_floor[r]))    # having x = Qth and c = 1/COP minimizes the electric consumption 

    I       = np.identity(samples)
    
    """Setting matrix and vector size of inequality constraints. """
    Aub_opt = np.zeros([(Nc_tot)*samples,n_Bui*samples]) # Matrix of inequality constraints
    Bub_opt = np.zeros([(Nc_tot)*samples,1])             # Vector of inequality constraints
    
    """Compilation of the matrix and vector of inequality constraints. """
    for i in range (n_SFH06_air):
        for t in range (samples):
            """2006-today archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i)*2)*samples), ((i)*samples):((i + 1)*samples)] = np.array(Z["Z_SFH06_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i)*2)*samples), 0] = float(upper["upper_SFH06_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH06_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i)*2)+1)*samples), ((i)*samples):((i + 1)*samples)] = np.array(Z["Z_SFH06_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH06_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH06_{0}".format(i + 1)][t])
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot-Nc_cluster)*samples), ((i)*samples):((i + 1)*samples)] = I[t,:]*(1/COP_SFH06_air[t])
    
    for i in range (n_SFH05_air):
        for t in range (samples):
            """1991-2005 archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air)*2)*samples), ((i + n_SFH06_air)*samples):((i + n_SFH06_air + 1)*samples)] = np.array(Z["Z_SFH05_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air)*2)*samples), 0] = float(upper["upper_SFH05_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH05_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air)*2)+1)*samples), ((i + n_SFH06_air)*samples):((i + n_SFH06_air + 1)*samples)] = np.array(Z["Z_SFH05_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH05_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH05_{0}".format(i + 1)][t])
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot-Nc_cluster)*samples), ((i + n_SFH06_air)*samples):((i + n_SFH06_air + 1)*samples)] = I[t,:]*(1/COP_SFH05_air[t])

    for i in range (n_SFH90_air):
        for t in range (samples):
            """1976-1990 archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air)*2)*samples), ((i + n_SFH06_air + n_SFH05_air)*samples):((i + n_SFH06_air + n_SFH05_air + 1)*samples)] = np.array(Z["Z_SFH90_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air)*2)*samples), 0] = float(upper["upper_SFH90_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH90_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air)*2)+1)*samples), ((i + n_SFH06_air + n_SFH05_air)*samples):((i + n_SFH06_air + n_SFH05_air + 1)*samples)] = np.array(Z["Z_SFH90_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH90_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH90_{0}".format(i + 1)][t])
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot-Nc_cluster)*samples), ((i + n_SFH06_air + n_SFH05_air)*samples):((i + n_SFH06_air + n_SFH05_air + 1)*samples)] = I[t,:]*(1/COP_SFH90_air[t])

    for i in range (n_SFH60_air):
        for t in range (samples):
            """1946-1960 archetype with air systems indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = np.array(Z["Z_SFH60_{0}".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)*samples), 0] = float(upper["upper_SFH60_{0}".format(i + 1)][t]) - float(Mmax["Mmax_SFH60_{0}".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)+1)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = np.array(Z["Z_SFH60_{0}".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH60_{0}".format(i+ 1)][t][0]) - float(lower["lower_SFH60_{0}".format(i + 1)][t])
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot-Nc_cluster)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = I[t,:]*(1/COP_SFH60_air[t])

    for i in range (n_SFH06_floor):
        for t in range (samples):
            """2006-today archetype with underfloor heating system indoor air temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air+ 1)*samples)] = np.array(Z["Z_SFH06_{0}_floor".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2)*samples), 0] = float(upper["upper_SFH06_{0}_floor".format(i + 1)][t]) - float(Mmax["Mmax_SFH06_{0}_floor".format(i + 1)][t][0])
    
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+1)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air+ 1)*samples)] = np.array(Z["Z_SFH06_{0}_floor".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+1)*samples), 0] = float(Mmax["Mmax_SFH06_{0}_floor".format(i + 1)][t][0]) - float(lower["lower_SFH06_{0}_floor".format(i + 1)][t])
            
            """2006-today archetype with underfloor heating system floor temperature constraints. """
            Aub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2 + 2)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air+ 1)*samples)] = np.array(Z_floor["Z_SFH06_{0}_f".format(i + 1)])[t,:]
            Bub_opt[t+(((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air)*2 + 2)*samples), 0] = float(29) - float(Mmax_floor["Mmax_SFH06_{0}_f".format(i + 1)][t][0])
            
            Aub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+ 2 + 1)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)] = np.array(Z_floor["Z_SFH06_{0}_f".format(i + 1)])[t,:]*(-1)
            Bub_opt[t+((((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*2)+ 2 + 1)*samples), 0] = float(Mmax_floor["Mmax_SFH06_{0}_f".format(i+ 1)][t][0]) - float(15)
            
            """Constraint on power consumption at cluster level. """
            Aub_opt[t+((Nc_tot-Nc_cluster)*samples), ((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + n_SFH60_air + 1)*samples)] = I[t,:]*(1/COP_SFH06_floor[t])
    
    """Constraint on power consumption at cluster level. """
    for t in range (samples):
        Bub_opt[t+((Nc_bui)*samples), 0] = Eele_max[t]      # constrained electric consumption at cluster level
    
    Aub = Aub_opt.reshape((Nc_tot)*samples,n_Bui*samples)
    Bub = Bub_opt.reshape((Nc_tot)*samples,1)
    """Optimization result (Demand-Response scenario). """
    res = linprog(c=c, A_ub=Aub_opt, b_ub=Bub_opt, bounds=bounds, method='highs-ipm', options={'maxiter': 1000000})
    print("Cluster level:")
    print(res.message + "(Demand-Response scenario)")
    if res.success == False:
        print("\n              ")
        print("\nWARNING")
        print("\nThe thermal demand of the building (decision variable 1) depends on the state of the indoor air temperature Tair.")
        print("\nTo enable the resolution of the problem, it is recommended to change the tolerances set in the thermostat")
        print("\nThis changes the temperature range in which the indoor air temperature is constrained. ")
        print("\nAlternatively, we recommend checking the size of the heat pump for a change in the bounds of the decision variable. ")
    
    """Calculation of thermal loads (decision variable x1), 
            indoor air temperatures and building electricity consumption during the DR scenario. """
            
    Qres_dr = {}  # Values of the decision variable x1 (thermal demand in Wh)
    Tair_dr = {}  # Indoor air temeperature (°C) - Demand-Response scenario
    Q_dr    = {}     # Thermal demand (kWh) - Demand-Response scenario
    E_dr    = {}     # Electric demand (kWh) - Demand-Response scenario
    cost_dr = {}     # Electricity bill (Eur) - Demand-Response scenario
    
    for i in range(n_SFH06_air):
        """Values of the decision variables x1 (i.e., thermal load). """
        Qres_dr["Qres_dr_SFH06_{0}".format(i+1)] = (res.x[(i*samples):((i + 1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH06_{0}".format(i+1)] = list(M["M_SFH06_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}".format(i+1)], Qres_dr["Qres_dr_SFH06_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH06_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH06_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption drawn from the grid (DR scenario). """
        E_dr["E_dr_SFH06_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH06_{0}".format(i+1)][r]/COP_SFH06_air[r] for r in range(0, samples)]
        
        """Calculating the cost of electiricty drawn from the grid (DR scenario). """
        cost_dr["cost_dr_SFH06_{0}".format(i+1)] = [E_dr["E_dr_SFH06_{0}".format(i+1)][r]*ele_grid_price[r] for r in range(0, samples)]
        
    for i in range(n_SFH05_air):
        """Values of the decision variables x1 (i.e., thermal load). """
        Qres_dr["Qres_dr_SFH05_{0}".format(i+1)] = (res.x[((i+n_SFH06_air)*samples):((i+n_SFH06_air+1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH05_{0}".format(i+1)] = list(M["M_SFH05_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH05_{0}".format(i+1)], Qres_dr["Qres_dr_SFH05_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH05_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH05_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption drawn from the grid (DR scenario). """
        E_dr["E_dr_SFH05_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH05_{0}".format(i+1)][r]/COP_SFH05_air[r] for r in range(0, samples)]
        
        """Calculating the cost of electiricty drawn from the grid (DR scenario). """
        cost_dr["cost_dr_SFH05_{0}".format(i+1)] = [E_dr["E_dr_SFH05_{0}".format(i+1)][r]*ele_grid_price[r] for r in range(0, samples)]
        
    for i in range(n_SFH90_air):
        """Values of the decision variables x1 (i.e., thermal load). """
        Qres_dr["Qres_dr_SFH90_{0}".format(i+1)] = (res.x[((i+n_SFH06_air+n_SFH05_air)*samples):((i+n_SFH06_air+n_SFH05_air+1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH90_{0}".format(i+1)] = list(M["M_SFH90_{0}".format(i+1)][6, :].reshape(samples, 1) + np.dot(Z["Z_SFH90_{0}".format(i+1)], Qres_dr["Qres_dr_SFH90_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH90_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH90_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption drawn from the grid (DR scenario). """
        E_dr["E_dr_SFH90_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH90_{0}".format(i+1)][r]/COP_SFH90_air[r] for r in range(0, samples)]
        
        """Calculating the cost of electiricty drawn from the grid (DR scenario). """
        cost_dr["cost_dr_SFH90_{0}".format(i+1)] = [E_dr["E_dr_SFH90_{0}".format(i+1)][r]*ele_grid_price[r] for r in range(0, samples)]
        
    for i in range(n_SFH60_air):
        """Values of the decision variables x1 (i.e., thermal load). """
        Qres_dr["Qres_dr_SFH60_{0}".format(i+1)] = (res.x[((i + n_SFH06_air + n_SFH05_air + n_SFH90_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air + 1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH60_{0}".format(i+1)] = list(M["M_SFH60_{0}".format(i+1)][3, :].reshape(samples, 1) + np.dot(Z["Z_SFH60_{0}".format(i+1)], Qres_dr["Qres_dr_SFH60_{0}".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH60_{0}".format(i+1)]       = [Qres_dr["Qres_dr_SFH60_{0}".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption drawn from the grid (DR scenario). """
        E_dr["E_dr_SFH60_{0}".format(i+1)]       = [Q_dr["Q_dr_SFH60_{0}".format(i+1)][r]/COP_SFH60_air[r] for r in range(0, samples)]
        
        """Calculating the cost of electiricty drawn from the grid (DR scenario). """
        cost_dr["cost_dr_SFH60_{0}".format(i+1)] = [E_dr["E_dr_SFH60_{0}".format(i+1)][r]*ele_grid_price[r] for r in range(0, samples)]
        
    for i in range(n_SFH06_floor):
        """Values of the decision variables x1 (i.e., thermal load). """
        Qres_dr["Qres_dr_SFH06_{0}_floor".format(i+1)] = (res.x[((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air)*samples):((i + n_SFH06_air + n_SFH05_air + n_SFH90_air+ n_SFH60_air + 1)*samples)].reshape(samples,1))
        
        """Indoor air temperature calculation (DR scenario). """
        Tair_dr["Tair_dr_SFH06_{0}_floor".format(i+1)] = list(M["M_SFH06_{0}_floor".format(i+1)][7, :].reshape(samples, 1) + np.dot(Z["Z_SFH06_{0}_floor".format(i+1)], Qres_dr["Qres_dr_SFH06_{0}_floor".format(i+1)]))
        
        """W to kW thermal power conversion. """
        Q_dr["Q_dr_SFH06_{0}_floor".format(i+1)]       = [Qres_dr["Qres_dr_SFH06_{0}_floor".format(i+1)][r, 0]/1000 for r in range(0, samples)]
        
        """Calculating electric power consumption drawn from the grid (DR scenario). """
        E_dr["E_dr_SFH06_{0}_floor".format(i+1)]       = [Q_dr["Q_dr_SFH06_{0}_floor".format(i+1)][r]/COP_SFH06_floor[r] for r in range(0, samples)]
        
        """Calculating the cost of electiricty drawn from the grid (DR scenario). """
        cost_dr["cost_dr_SFH06_{0}_floor".format(i+1)] = [E_dr["E_dr_SFH06_{0}_floor".format(i+1)][r]*ele_grid_price[r] for r in range(0, samples)]
        

    """Conversion of dictionary to data frame for calculation of aggregate electric power consumption. """
    E_dr = pd.DataFrame.from_dict(E_dr)
    
    """Calculating electric consumption at the aggregate level during DR scenario. """
    Etot_dr = E_dr.sum(axis=1)
    Edr_sum = float(Etot_dr.sum())
    Etot_dr = [Etot_dr[i] for i in range(0,samples)]  
    
    """Conversion of dictionary to data frame for calculation of aggregate electricity bill. """
    cost_dr = pd.DataFrame.from_dict(cost_dr)
    
    """Calculating electricity bill at the aggregate level during DR scenario. """
    cost_tot_dr = cost_dr.sum(axis=1)
    cost_dr_sum = float(cost_tot_dr.sum())
    cost_tot_dr = [cost_tot_dr[i] for i in range(0,samples)]  
    
    
    import matplotlib.pyplot as plt
    time = [i*timestep for i in range(0, samples)] # time array
    
    """Do plot. """
    if do_plot_E_dr_cluster:
        """If True, plot cluster electric and PV consumption trends. """
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Time (hours)', fontsize=14)
        ax1.set_ylabel('$\.P_\mathrm{cluster}$ (kW$_\mathrm{el}$)', fontsize=14)
        ax1.plot(time, Etot_bl, 'k', label="BL_scenario", linewidth=1.0)
        ax1.plot(time, Etot_dr, 'r', label="DR_scenario", linewidth=1.0)
        ax1.set_xlim(0, duration-1)
        ax1.grid(alpha=0.3)
        plt.legend(frameon=False, fontsize=12)
        fig.set_size_inches(w=10, h=4)
        fig.tight_layout()
        plt.savefig("Eele_cluster", dpi=200)
        
        plt.show()
        plt.close('all')

    if do_plot_E_dr_single:
        """If True, plot single buildings electric and PV consumption trends. """
        for i in range (n_SFH06_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            ax1.plot(time,E_dr["E_dr_SFH06_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH06_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH06_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH06_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
            
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH05_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            ax1.plot(time,E_dr["E_dr_SFH05_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH05_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH05_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH05_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH90_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            ax1.plot(time,E_dr["E_dr_SFH90_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH90_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH90_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH90_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH60_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            ax1.plot(time,E_dr["E_dr_SFH60_{0}".format(i+1)],'r', linewidth=1, label = "E_dr_SFH60_{0}".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH60_{0}".format(i+1)],'k', linewidth=1, label = "E_bl_SFH60_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
    
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH06_floor):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('$\.P_\mathrm{building}$ (kW$_\mathrm{el}$)', fontsize=14)
            ax1.plot(time,E_dr["E_dr_SFH06_{0}_floor".format(i+1)],'r', linewidth=1, label = "E_dr_SFH06_{0}_floor".format(i+1))
            ax1.plot(time,E_bl["E_bl_SFH06_{0}_floor".format(i+1)],'k', linewidth=1, label = "E_bl_SFH06_{0}_floor".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
            
        plt.show()
        plt.close('all')

    if do_plot_Tair_dr_single:
        """If True, plot single buildings indoor air temperature trends. """
        for i in range (n_SFH06_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            ax1.plot(time,Tair_dr["Tair_dr_SFH06_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH06_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH06_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH06_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
            
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH05_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            ax1.plot(time,Tair_dr["Tair_dr_SFH05_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH05_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH05_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH05_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH90_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            ax1.plot(time,Tair_dr["Tair_dr_SFH90_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH90_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH90_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH90_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
        for i in range (n_SFH60_air):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            ax1.plot(time,Tair_dr["Tair_dr_SFH60_{0}".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH60_{0}".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH60_{0}".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH60_{0}".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
    
        plt.show()
        plt.close('all')
        
        for i in range (n_SFH06_floor):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Time (hours)')
            ax1.set_ylabel('Indoor air temperature (°C)', fontsize=14)
            ax1.plot(time,Tair_dr["Tair_dr_SFH06_{0}_floor".format(i+1)],'r--', linewidth=1, label = "Tair_dr_SFH06_{0}_floor".format(i+1))
            ax1.plot(time,Tair_bl["Tair_bl_SFH06_{0}_floor".format(i+1)],'k--', linewidth=1, label = "Tair_bl_SFH06_{0}_floor".format(i+1))

            ax1.set_xlim(0, duration-1)
            ax1.grid(alpha=0.3)
            plt.legend(frameon=False, fontsize=12)
            fig.set_size_inches(w=10, h=4)
            fig.tight_layout()
        
        plt.show()
        plt.close('all')
    
    """Data export to excel. """
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writerEele = pd.ExcelWriter(name_file + '.xlsx', engine='xlsxwriter')

    # Write each dataframe to a different worksheet.
    E_bl.to_excel(writerEele, sheet_name='Eele_BL')
    E_dr.to_excel(writerEele, sheet_name='Eele_DR')

    # Close the Pandas Excel writer and output the Excel file.
    writerEele.close()
    
    """Electricity consumption during BL and DR scenarios. """
    print("\nCluster electricity consumption of about " +
          str(round(Ebl_sum,2)) + " kWh during baseline (BL) scenario" +
          "\n and a cluster electricity bill of about " + str(round(cost_bl_sum,2)) + " Eur.")
    print("\n------------------------------------------------------------------------")
    
    ele_Diff   = ((Edr_sum - Ebl_sum)/Ebl_sum)*100
    
    if ele_Diff >= 0:
        print("\nCluster electricity consumption of about " +
              str(round(Edr_sum,2)) + " kWh during DR scenario with load-shifting strategy,"
              "\nwith a " + str(round(ele_Diff,2)) + "% increase compared to BL scenario.")
        print("\n------------------------------------------------------------------------")
    else:
        print("\nCluster electricity consumption of about " +
              str(round(Edr_sum,2)) + " kWh during DR scenario with load-shifting strategy,"
              "\nwith a " + str(round(ele_Diff*(-1),2)) + "% reduction compared to BL scenario.")
        print("\n------------------------------------------------------------------------")
    
    cost_Diff   = ((cost_dr_sum - cost_bl_sum)/cost_bl_sum)*100
    
    if cost_Diff >= 0:
        print("\nCluster electricity bill of about " +
              str(round(cost_dr_sum,2)) + " Eur during DR scenario with load-shifting strategy,"
              "\nwith a " + str(round(cost_Diff,2)) + "% increase compared to BL scenario.")
        print("\n------------------------------------------------------------------------")
    else:
        print("\nCluster electricity bill of about " +
              str(round(cost_dr_sum,2)) + " Eur during DR scenario with load-shifting strategy,"
              "\nwith a " + str(round(cost_Diff*(-1),2)) + "% reduction compared to BL scenario.")
        print("\n------------------------------------------------------------------------")
    
    
    return(Etot_dr)

