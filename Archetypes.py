"""
ClustEnergy OpTool version 1.0.0
Released on February 29, 2024
@author: DIISM UNIVPM
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
import numpy as np

from SingleFamilyHouses import SFH_06, SFH_05, SFH_90, SFH_60 
from Meteo import MeteoFile,Window_solar,Wall_solar
from Functions import optForm_build


#---------------------------------------------------------------------------------------------------------------------------------------
"""Check if the following python libraries needed to run the simulations are installed:
       numpy and math. """
#---------------------------------------------------------------------------------------------------------------------------------------

"""Building archetypes are defined below. 
    According to the preferences, matrices (Aub) and vectors (Bub) of the inequality constraints, coefficients (c) and maximum and minimum pair (bounds) 
    of the decision variable (c) are compiled, which are useful for calculating the thermal loads of the archetype during the baseline scenario. 
    
    In addition, for the calculation of thermal demand during the Demand-Response scenario, useful variables are output to describe the state of the 
    indoor air temperature node (M, Z and Mmax) and the upper (Upper) and lower (Lower) ranges of the temperature profiles set in the thermostat. """

GR  = 0.2                # ground reflectance
wall_absorption   = 0.3  # wall absorption

radiative_part  = 0.6    # radiative part of internal gains  
convective_part = 0.4    # convective part of internal gains  
#----------------------------------------------------------------------------------------------------
"""
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
2006-today Single Family House with air cooling/heating system
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
def SFH2006_air(Locality, day, month, duration, timestep, th_sp, th_tolerance_sup, th_tolerance_min, thermal_load, occupancy, internalGains, Q_hp):
    """Recall of variables useful for defining inputs. """
    T_est           =  MeteoFile(Locality,int(day),int(month),duration,timestep)[0]  # Outdoor temperature trends (°C)
    T_ground        =  MeteoFile(Locality,int(day),int(month),duration,timestep)[1]  # Ground temperature trends (°C)
    Diff_Whm2       =  MeteoFile(Locality,int(day),int(month),duration,timestep)[2]  # Diffuse solar radiation (Whm2)
    Dir_Whm2        =  MeteoFile(Locality,int(day),int(month),duration,timestep)[3]  # Direct solar radiation (Whm2)
    Glo_Whm2        =  MeteoFile(Locality,int(day),int(month),duration,timestep)[4]  # Global solar radiation (Whm2)
    sun_altitude    =  MeteoFile(Locality,int(day),int(month),duration,timestep)[5]  # Sun altitude
    sun_azimut      =  MeteoFile(Locality,int(day),int(month),duration,timestep)[6]  # Sun azimut
    sun_zenit       =  MeteoFile(Locality,int(day),int(month),duration,timestep)[7]  # Sun zenit
    latitude_deg    =  MeteoFile(Locality,int(day),int(month),duration,timestep)[8]  # Locality laitude
    longitude_deg   =  MeteoFile(Locality,int(day),int(month),duration,timestep)[9]  # Locality longitude
    start           =  MeteoFile(Locality,int(day),int(month),duration,timestep)[10] # Simulation start time
    UTC             =  MeteoFile(Locality,int(day),int(month),duration,timestep)[11] # Locality Coordinated Universal Time

    samples =  int(duration/timestep) # number of total steps
    
    """Calling 2006-today Single-Family House building class. """
    Archetype = SFH_06(ACH = 0.2, wall_absorption = wall_absorption, radiative_part = radiative_part, convective_part = convective_part)
    
    """Defining solar and internal gains. """
    southWindow  = Window_solar(azimut_surface=0.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window, area = Archetype.S_windows_area)   
    eastWindow  = Window_solar(azimut_surface=270.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window,  area = Archetype.E_windows_area) 
    westWindow  = Window_solar(azimut_surface=90.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window,  area = Archetype.W_windows_area)  
    northWindow  = Window_solar(azimut_surface=180.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window, area = Archetype.N_windows_area) 
    
    Archetype.Uvalues()
    southWall    = Wall_solar(azimut_surface=0.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0) , equivalent_transmittance = Archetype.Uwall, area = Archetype.south_area)
    eastWall     = Wall_solar(azimut_surface=270.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0),  equivalent_transmittance = Archetype.Uwall,area = Archetype.east_area)       
    westWall     = Wall_solar(azimut_surface=90.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0), equivalent_transmittance = Archetype.Uwall, area = Archetype.west_area)        
    northWall    = Wall_solar(azimut_surface=180.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0) , equivalent_transmittance = Archetype.Uwall , area = Archetype.north_area)    
    horRoof      = Wall_solar(azimut_surface=0.0, inclination_surface = 0.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0) , equivalent_transmittance = Archetype.Uroof , area = Archetype.roof_area)                            
    
    window_solar = []
    walls_solar  = []
    roof_solar   = []

    for h in range(0,len(sun_altitude)):     
        southWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        eastWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        westWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        northWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        southWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        eastWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        westWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        northWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        horRoof.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        """Solar gains to give in input to the state-space model formulation. """
        window_solar.append(southWindow.windows_solar_gains + eastWindow.windows_solar_gains + westWindow.windows_solar_gains + northWindow.windows_solar_gains) # Windows solar gains
        walls_solar.append(southWall.walls_solar_gains + eastWall.walls_solar_gains + westWall.walls_solar_gains + northWall.walls_solar_gains)                  # Walls solar gains
        roof_solar.append(horRoof.walls_solar_gains)                                                                                                             # Roof solar gains

    floor_area = Archetype.floor_area
    
    """Internal gains due to people, lighting, and equipment are calculated from the occupancy profile defined using the richardsonpy tool
        (see the python file "User_pattern.py"). Then the trends are given as input to the archetype function. 
        However, it is possible to consider them constant during the simulation period calculated according to the methodology described in UNI 11300-1. """
    
    # To consider internal gains constant, then activate the following lines of code
    """
    Af = floor_area*2   
    if Af<=120: 
       internalGains = 7.987 * Af - 0.0353 * Af^2  # Internal gains (UNI 11300-1 pr. 13) 
    else: 
       internalGains = 450
    """
    
    """Defining the upper (Upper) and lower (Lower) ranges of the temperature profiles set in the thermostat. """
    Upper = []
    Lower = []

    if thermal_load == 'cooling':
        for i in range (samples):
            if occupancy[i] > 0:
                Upper.append(float(th_sp[i])  + th_tolerance_sup)
            else:
                Upper.append(30) # meaning the shutdown of the cooling system (active cooling)
            Lower.append(float(th_sp[i]) - th_tolerance_min)
    elif thermal_load == 'heating':
        for i in range (samples):
            if occupancy[i] > 0:
                Lower.append(float(th_sp[i]) - th_tolerance_min)
            else:
                Lower.append(float(th_sp[i]) - th_tolerance_min)
                #Lower.append(10) # meaning the shutdown of the heating system (active heating)
            Upper.append(float(th_sp[i])  + th_tolerance_sup)
    
    """Configuring the state-space model representation for linear programming problem resolution. """
    if thermal_load == 'cooling':
        Archetype.stateSpace_air_cooling() # state-space representation for air cooling system
    elif thermal_load == 'heating':
        Archetype.stateSpace_air_heating() # state-space representation for air heating system
    
    U = np.empty([6,samples]) # input vector of the state-space model formulation
   
    for r in range(0,len(sun_altitude)):       
        U[0,r] = T_est[r]                  # Outdoor temperature (°C)
        U[1,r] = T_ground[r]               # Ground temperature (°C)
        U[2,r] = window_solar[r]           # Solar Gains Windows (W) 
        U[3,r] = walls_solar[r]            # Solar Gains-Wall (W) 
        U[4,r] = roof_solar[r]             # Solar Gains-Roof (W)  
        U[5,r] = internalGains[r]          # Internal gains (W)    

    Xo     = np.array([[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]]]) # initial state of indoor air temperature

    Archetype.sys_dis_forOp(Archetype.sys_tuple,timestep) # System discretization for linear resolution
    
    Ad = np.asarray(Archetype.sys_disc_A) # State matrix A of the discretized system written in SSM form
    Bd = np.asarray(Archetype.sys_disc_B) # Input matrix B of the discretized system written in SSM form
    Ed = np.asarray(Archetype.sys_disc_E) # State matrix E of the discretized system isolating the decision variable 

    Z     = optForm_build(samples, Ad, Xo, Bd, U, Ed, 6)[0] # Variable describing the state of the system containing the decision variable (e.g., thermal demand)
    M_max = optForm_build(samples, Ad, Xo, Bd, U, Ed, 6)[1] # Variable describing the state of the indoor air temperture
    M     = optForm_build(samples, Ad, Xo, Bd, U, Ed, 6)[2]    
    
    """Defining linear programming problem. """
    Nc         = 2  # Number of constraints (2 for indoor air temperature)

    c         = [] # Decision variable coefficient 
    bounds    = [] # Decision variable bounds

    for r in range(0,samples):
        bounds.append((0,Q_hp[r]))
        c.append(timestep)

    Aub_opt = np.zeros([(Nc)*samples,samples]) # Matrix of inequality constraints
    Bub_opt = np.zeros([(Nc)*samples,1])       # Vector of inequality constraints
    
    """Compiling matrix and vector of the inequality constraints. 
        The constraint imposed in the optimization refers only to the state of indoor air temperature (Tair). """
    """Temperature constraints on the state of indoor air temperature system. """
    for i in range(0,samples):
        # Tair <= Upper
        Aub_opt[i,0:(samples)] = Z[i,:]
        Bub_opt[i,0]           = float(Upper[i])-M_max[i]      
    for i in range(0,samples): 
        # Tair >= Lower
        Aub_opt[i+(1*samples),0:(samples)] = Z[i,:]*(-1)
        Bub_opt[i+(1*samples),0]           = M_max[i]-float(Lower[i])

    Aub = Aub_opt.reshape(samples*(Nc),samples)
    Bub = Bub_opt.reshape(samples*(Nc),1)  
       
    return(c,      # 0 - decision variable coefficients
           Aub,    # 1 - Matrix of inequality constraints
           Bub,    # 2 - Vector of inequlity constraints
           bounds, # 3 - Decision variable bounds
           M,      # 4 - Variable describing the state of the indoor air temperture
           Z,      # 5 - Variable describing the state of the system containing the thermal demand (decision variable)
           M_max,  # 6 - Variable describing the state of the indoor air temperture
           Upper,  # 7 - Upper range of indoor temperature trend
           Lower,  # 8 - Lower range of indoor temperature trend
           floor_area) # 9 - Floor area (m2)

"""
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
2006-today Single Family House with underfloor heating system
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
def SFH2006_floor(Locality, day, month, duration, timestep, th_sp, th_tolerance_sup, th_tolerance_min, occupancy, internalGains,Q_hp):     
    """Recall of variables useful for defining inputs. """
    T_est           =  MeteoFile(Locality,int(day),int(month),duration,timestep)[0]  # Outdoor temperature trends (°C)
    T_ground        =  MeteoFile(Locality,int(day),int(month),duration,timestep)[1]  # Ground temperature trends (°C)
    Diff_Whm2       =  MeteoFile(Locality,int(day),int(month),duration,timestep)[2]  # Diffuse solar radiation (Whm2)
    Dir_Whm2        =  MeteoFile(Locality,int(day),int(month),duration,timestep)[3]  # Direct solar radiation (Whm2)
    Glo_Whm2        =  MeteoFile(Locality,int(day),int(month),duration,timestep)[4]  # Global solar radiation (Whm2)
    sun_altitude    =  MeteoFile(Locality,int(day),int(month),duration,timestep)[5]  # Sun altitude
    sun_azimut      =  MeteoFile(Locality,int(day),int(month),duration,timestep)[6]  # Sun azimut
    sun_zenit       =  MeteoFile(Locality,int(day),int(month),duration,timestep)[7]  # Sun zenit
    latitude_deg    =  MeteoFile(Locality,int(day),int(month),duration,timestep)[8]  # Locality laitude
    longitude_deg   =  MeteoFile(Locality,int(day),int(month),duration,timestep)[9]  # Locality longitude
    start           =  MeteoFile(Locality,int(day),int(month),duration,timestep)[10] # Simulation start time
    UTC             =  MeteoFile(Locality,int(day),int(month),duration,timestep)[11] # Locality Coordinated Universal Time
    
    samples =  int(duration/timestep) # number of total steps
    
    """Calling 2006-today Single-Family House building class. """
    Archetype = SFH_06(ACH = 0.2, wall_absorption = wall_absorption, radiative_part = radiative_part, convective_part = convective_part)
    
    """Defining solar and internal gains. """
    southWindow  = Window_solar(azimut_surface=0.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window, area = Archetype.S_windows_area)   
    eastWindow  = Window_solar(azimut_surface=270.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window,  area = Archetype.E_windows_area) 
    westWindow  = Window_solar(azimut_surface=90.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window,  area = Archetype.W_windows_area)  
    northWindow  = Window_solar(azimut_surface=180.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window, area = Archetype.N_windows_area) 
    
    Archetype.Uvalues()
    southWall    = Wall_solar(azimut_surface=0.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0) , equivalent_transmittance = Archetype.Uwall, area = Archetype.south_area)
    eastWall     = Wall_solar(azimut_surface=270.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0),  equivalent_transmittance = Archetype.Uwall,area = Archetype.east_area)       
    westWall     = Wall_solar(azimut_surface=90.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0), equivalent_transmittance = Archetype.Uwall, area = Archetype.west_area)        
    northWall    = Wall_solar(azimut_surface=180.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0) , equivalent_transmittance = Archetype.Uwall , area = Archetype.north_area)    
    horRoof      = Wall_solar(azimut_surface=0.0, inclination_surface = 0.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0) , equivalent_transmittance = Archetype.Uroof , area = Archetype.roof_area)                            
    
    window_solar = []
    walls_solar  = []
    roof_solar   = []

    for h in range(0,len(sun_altitude)):     
        southWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        eastWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        westWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        northWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        southWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        eastWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        westWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        northWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        horRoof.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        """Solar gains to give in input to the state-space model formulation. """
        window_solar.append(southWindow.windows_solar_gains + eastWindow.windows_solar_gains + westWindow.windows_solar_gains + northWindow.windows_solar_gains)   # Solar gains
        walls_solar.append(southWall.walls_solar_gains + eastWall.walls_solar_gains + westWall.walls_solar_gains + northWall.walls_solar_gains)                    # Solar gains
        roof_solar.append(horRoof.walls_solar_gains)      
                                                                                                             
    floor_area = Archetype.floor_area
    
    """Internal gains due to people, lighting, and equipment are calculated from the occupancy profile defined using the richardsonpy tool
        (see the python file "User_pattern.py"). Then the trends are given as input to the archetype function. 
        However, it is possible to consider them constant during the simulation period calculated according to the methodology described in UNI 11300-1. """
    
    # To consider internal gains constant, then activate the following lines of code
    """
    Af = floor_area*2   
    if Af<=120: 
       internalGains = 7.987 * Af - 0.0353 * Af^2  # Internal gains (UNI 11300-1 pr. 13) 
    else: 
       internalGains = 450
    """
    
    """Defining the upper (Upper) and lower (Lower) ranges of the temperature profiles set in the thermostat. """
    Upper = []
    Lower = []

    for i in range (samples):
        if occupancy[i] > 0:
            Lower.append(float(th_sp[i]) - th_tolerance_min)
        else:
            Lower.append(float(th_sp[i]) - th_tolerance_min)
            #Lower.append(th_sp[i]) - th_tolerance_min - 2) # active heating
        Upper.append(float(th_sp[i])  + th_tolerance_sup)
    
    """Configuring the state-space model representation for linear programming problem resolution. """
    Archetype.stateSpace_floor_heating() # state-space representation for underfloor heating system
    
    U = np.empty([6,samples]) # input vector of the state-space model formulation
   
    for r in range(0,len(sun_altitude)):       
        U[0,r] = T_est[r]                  # Outdoor temperature (°C)
        U[1,r] = T_ground[r]               # Ground temperature (°C)
        U[2,r] = window_solar[r]           # Solar Gains Windows (W) 
        U[3,r] = walls_solar[r]            # Solar Gains-Wall (W) 
        U[4,r] = roof_solar[r]             # Solar Gains-Roof (W)  
        U[5,r] = internalGains[r]          # Internal gains (W)    

    Xo     = np.array([[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]]]) # initial state of indoor air temperature

    Archetype.sys_dis_forOp(Archetype.sys_tuple,timestep) # System discretization for linear resolution
    
    Ad = np.asarray(Archetype.sys_disc_A) # State matrix A of the discretized system written in SSM form
    Bd = np.asarray(Archetype.sys_disc_B) # Input matrix B of the discretized system written in SSM form
    Ed = np.asarray(Archetype.sys_disc_E) # State matrix E of the discretized system isolating the decision variable 

    Z     = optForm_build(samples, Ad, Xo, Bd, U, Ed, 7)[0] # Variable describing the state of the system containing the thermal demand (decision variable)
    M_max = optForm_build(samples, Ad, Xo, Bd, U, Ed, 7)[1] # Variable describing the state of the indoor air temperture
    M     = optForm_build(samples, Ad, Xo, Bd, U, Ed, 7)[2]    
    
    Z_floor     = optForm_build(samples, Ad, Xo, Bd, U, Ed, 6)[0] # Variable describing the state of the floor system containing the thermal demand (decision variable)
    M_max_floor = optForm_build(samples, Ad, Xo, Bd, U, Ed, 6)[1] # Variable describing the state of the floor temperture
    
    """Defining linear programming problem. """
    Nc         = 2 + 2  # Number of constraints (2 for indoor air temperature + 2 for floor temperature)

    c         = [] # Decision variable coefficient 
    bounds    = [] # Decision variable bounds

    for r in range(0,samples):
        bounds.append((0,Q_hp[r]))
        c.append(timestep)
    
    Aub_opt = np.zeros([(Nc)*samples,samples]) # Matrix of inequality constraints
    Bub_opt = np.zeros([(Nc)*samples,1])       # Vector of inequality constraints
    
    """Compiling matrix and vector of the inequality constraints. 
        The constraint imposed in the optimization refers to the state of indoor air temperature (Tair)
        and floor temperature (Tf). """
    """Temperature constraints on the state of indoor air temperature system. """
    for i in range(0,samples):
        # Tair <= Upper
        Aub_opt[i,0:(samples)] = Z[i,:]
        Bub_opt[i,0]           = Upper[i]-M_max[i]      
    for i in range(0,samples):
        # Tair >= Lower
        Aub_opt[i+(1*samples),0:(samples)] = Z[i,:]*(-1)
        Bub_opt[i+(1*samples),0]           = M_max[i]-Lower[i]
        
    """Temperature constraints on the state of radiant floor system. """
    for i in range(0,samples): 
        # Tair <= 29 
        Aub_opt[i+(2*samples),0:(samples)] = Z_floor[i,:]
        Bub_opt[i+(2*samples),0] = 29 - M_max_floor[i]
    for i in range(0,samples): 
        # Tair >= 15
        Aub_opt[i+(3*samples),0:(samples)] = Z_floor[i,:]*(-1)
        Bub_opt[i+(3*samples),0] = M_max_floor[i] - 15

    Aub = Aub_opt.reshape(samples*(Nc),samples)
    Bub = Bub_opt.reshape(samples*(Nc),1)     

    return(c,      # 0 - decision variable coefficients
           Aub,    # 1 - Matrix of inequality constraints
           Bub,    # 2 - Vector of inequlity constraints
           bounds, # 3 - Decision variable bounds
           M,      # 4 - Variable describing the state of the indoor air temperture
           Z,      # 5 - Variable describing the state of the system containing the thermal demand (decision variable)
           M_max,  # 6 - Variable describing the state of the indoor air temperture
           Upper,  # 7 - Upper range of indoor temperature trend
           Lower,  # 8 - Lower range of indoor temperature trend  # 8 - Lower range of indoor temperature trend
           Z_floor,     # 9 - Variable describing the state of the floor system containing the thermal demand (decision variable)
           M_max_floor, # 10 - Variable describing the state of the floor temperture
           floor_area)# 11 - Floor area (m2)

"""
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
1991-2005 Single Family House with air cooling/heating system
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
def SFH2005_air(Locality, day, month, duration, timestep, th_sp, th_tolerance_sup, th_tolerance_min, thermal_load, occupancy, internalGains, Q_hp):     
    """Recall of variables useful for defining inputs. """
    T_est           =  MeteoFile(Locality,int(day),int(month),duration,timestep)[0]  # Outdoor temperature trends (°C)
    T_ground        =  MeteoFile(Locality,int(day),int(month),duration,timestep)[1]  # Ground temperature trends (°C)
    Diff_Whm2       =  MeteoFile(Locality,int(day),int(month),duration,timestep)[2]  # Diffuse solar radiation (Whm2)
    Dir_Whm2        =  MeteoFile(Locality,int(day),int(month),duration,timestep)[3]  # Direct solar radiation (Whm2)
    Glo_Whm2        =  MeteoFile(Locality,int(day),int(month),duration,timestep)[4]  # Global solar radiation (Whm2)
    sun_altitude    =  MeteoFile(Locality,int(day),int(month),duration,timestep)[5]  # Sun altitude
    sun_azimut      =  MeteoFile(Locality,int(day),int(month),duration,timestep)[6]  # Sun azimut
    sun_zenit       =  MeteoFile(Locality,int(day),int(month),duration,timestep)[7]  # Sun zenit
    latitude_deg    =  MeteoFile(Locality,int(day),int(month),duration,timestep)[8]  # Locality laitude
    longitude_deg   =  MeteoFile(Locality,int(day),int(month),duration,timestep)[9]  # Locality longitude
    start           =  MeteoFile(Locality,int(day),int(month),duration,timestep)[10] # Simulation start time
    UTC             =  MeteoFile(Locality,int(day),int(month),duration,timestep)[11] # Locality Coordinated Universal Time
    
    samples =  int(duration/timestep) # number of total steps

    """Calling 1991-2005 Single-Family House building class. """
    Archetype = SFH_05(ACH = 0.2, wall_absorption = wall_absorption, radiative_part = radiative_part, convective_part = convective_part)

    """Defining solar and internal gains. """
    southWindow  = Window_solar(azimut_surface=0.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window, area = Archetype.S_windows_area)   
    eastWindow  = Window_solar(azimut_surface=270.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window,  area = Archetype.E_windows_area) 
    westWindow  = Window_solar(azimut_surface=90.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window,  area = Archetype.W_windows_area)  
    northWindow  = Window_solar(azimut_surface=180.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window, area = Archetype.N_windows_area) 
    
    Archetype.Uvalues()
    southWall    = Wall_solar(azimut_surface=0.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0) , equivalent_transmittance = Archetype.Uwall, area = Archetype.south_area)
    eastWall     = Wall_solar(azimut_surface=270.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0),  equivalent_transmittance = Archetype.Uwall,area = Archetype.east_area)       
    westWall     = Wall_solar(azimut_surface=90.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0), equivalent_transmittance = Archetype.Uwall, area = Archetype.west_area)        
    northWall    = Wall_solar(azimut_surface=180.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0) , equivalent_transmittance = Archetype.Uwall , area = Archetype.north_area)    
    horRoof      = Wall_solar(azimut_surface=0.0, inclination_surface = 0.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0) , equivalent_transmittance = Archetype.Uroof , area = Archetype.roof_area)                            
    
    window_solar = []
    walls_solar  = []
    roof_solar   = []

    for h in range(0,len(sun_altitude)):     
        southWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        eastWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        westWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        northWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        southWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        eastWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        westWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        northWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        horRoof.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        """Solar gains to give in input to the state-space model formulation. """
        window_solar.append(southWindow.windows_solar_gains + eastWindow.windows_solar_gains + westWindow.windows_solar_gains + northWindow.windows_solar_gains)   # Solar gains
        walls_solar.append(southWall.walls_solar_gains + eastWall.walls_solar_gains + westWall.walls_solar_gains + northWall.walls_solar_gains)                    # Solar gains
        roof_solar.append(horRoof.walls_solar_gains)      
    
    floor_area = Archetype.floor_area
    
    """Internal gains due to people, lighting, and equipment are calculated from the occupancy profile defined using the richardsonpy tool
        (see the python file "User_pattern.py"). Then the trends are given as input to the archetype function. 
        However, it is possible to consider them constant during the simulation period calculated according to the methodology described in UNI 11300-1. """
    
    # To consider internal gains constant, then activate the following lines of code
    """
    Af = floor_area*2   
    if Af<=120: 
       internalGains = 7.987 * Af - 0.0353 * Af^2  # Internal gains (UNI 11300-1 pr. 13) 
    else: 
       internalGains = 450
    """
    
    """Defining the upper (Upper) and lower (Lower) ranges of the temperature profiles set in the thermostat. """
    Upper = []
    Lower = []

    if thermal_load == 'cooling':
        for i in range (samples):
            if occupancy[i] > 0:
                Upper.append(float(th_sp[i])  + th_tolerance_sup)
            else:
                Upper.append(30) # meaning the shutdown of the cooling system (active cooling)
            Lower.append(float(th_sp[i]) - th_tolerance_min)
    elif thermal_load == 'heating':
        for i in range (samples):
            if occupancy[i] > 0:
                Lower.append(float(th_sp[i]) - th_tolerance_min)
            else:
                Lower.append(float(th_sp[i]) - th_tolerance_min)
                #Lower.append(10) # meaning the shutdown of the heating system (active heating)
            Upper.append(float(th_sp[i])  + th_tolerance_sup)
    
    """Configuring the state-space model representation for linear programming problem resolution. """
    if thermal_load == 'cooling':
        Archetype.stateSpace_air_cooling() # state-space representation for air cooling system
    elif thermal_load == 'heating':
        Archetype.stateSpace_air_heating() # state-space representation for air heating system
    
    U = np.empty([6,samples]) # input vector of the state-space model formulation
   
    for r in range(0,len(sun_altitude)):       
        U[0,r] = T_est[r]                  # Outdoor temperature (°C)
        U[1,r] = T_ground[r]               # Ground temperature (°C)
        U[2,r] = window_solar[r]           # Solar Gains Windows (W) 
        U[3,r] = walls_solar[r]            # Solar Gains-Wall (W) 
        U[4,r] = roof_solar[r]             # Solar Gains-Roof (W)  
        U[5,r] = internalGains[r]          # Internal gains (W)    

    Xo     = np.array([[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]]]) # initial state of indoor air temperature

    Archetype.sys_dis_forOp(Archetype.sys_tuple,timestep) # System discretization for linear resolution
    
    Ad = np.asarray(Archetype.sys_disc_A) # State matrix A of the discretized system written in SSM form
    Bd = np.asarray(Archetype.sys_disc_B) # Input matrix B of the discretized system written in SSM form
    Ed = np.asarray(Archetype.sys_disc_E) # State matrix E of the discretized system isolating the decision variable 
    
    Z     = optForm_build(samples, Ad, Xo, Bd, U, Ed, 6)[0] # Variable describing the state of the system containing the thermal demand (decision variable)
    M_max = optForm_build(samples, Ad, Xo, Bd, U, Ed, 6)[1] # Variable describing the state of the indoor air temperture
    M     = optForm_build(samples, Ad, Xo, Bd, U, Ed, 6)[2]    
    
    """Defining linear programming problem. """
    Nc         = 2  # Number of constraints (2 for indoor air temperature)

    c         = [] # Decision variable coefficient 
    bounds    = [] # Decision variable bounds

    for r in range(0,samples):
        bounds.append((0,Q_hp[r]))
        c.append(timestep)

    Aub_opt = np.zeros([(Nc)*samples,samples]) # Matrix of inequality constraints
    Bub_opt = np.zeros([(Nc)*samples,1])       # Vector of inequality constraints
    
    """Compiling matrix and vector of the inequality constraints. 
        The constraint imposed in the optimization refers only to the state of indoor air temperature (Tair). """
    """Temperature constraints on the state of indoor air temperature system. """
    for i in range(0,samples):
        # Tair <= Upper
        Aub_opt[i,0:(samples)] = Z[i,:]
        Bub_opt[i,0]           = float(Upper[i])-M_max[i]      
    for i in range(0,samples): 
        # Tair >= Lower
        Aub_opt[i+(1*samples),0:(samples)] = Z[i,:]*(-1)
        Bub_opt[i+(1*samples),0]           = M_max[i]-float(Lower[i])

    Aub = Aub_opt.reshape(samples*(Nc),samples)
    Bub = Bub_opt.reshape(samples*(Nc),1)     

    return(c,      # 0 - decision variable coefficients
           Aub,    # 1 - Matrix of inequality constraints
           Bub,    # 2 - Vector of inequlity constraints
           bounds, # 3 - Decision variable bounds
           M,      # 4 - Variable describing the state of the indoor air temperture
           Z,      # 5 - Variable describing the state of the system containing the thermal demand (decision variable)
           M_max,  # 6 - Variable describing the state of the indoor air temperture
           Upper,  # 7 - Upper range of indoor temperature trend
           Lower,  # 8 - Lower range of indoor temperature trend
           floor_area) # 9 - Floor area (m2)

"""
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
1976-1990 Single Family House with air cooling/heating system
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
def SFH1990_air(Locality, day, month, duration, timestep, th_sp, th_tolerance_sup, th_tolerance_min, thermal_load, occupancy, internalGains, Q_hp):     
    """Recall of variables useful for defining inputs. """
    T_est           =  MeteoFile(Locality,int(day),int(month),duration,timestep)[0]  # Outdoor temperature trends (°C)
    T_ground        =  MeteoFile(Locality,int(day),int(month),duration,timestep)[1]  # Ground temperature trends (°C)
    Diff_Whm2       =  MeteoFile(Locality,int(day),int(month),duration,timestep)[2]  # Diffuse solar radiation (Whm2)
    Dir_Whm2        =  MeteoFile(Locality,int(day),int(month),duration,timestep)[3]  # Direct solar radiation (Whm2)
    Glo_Whm2        =  MeteoFile(Locality,int(day),int(month),duration,timestep)[4]  # Global solar radiation (Whm2)
    sun_altitude    =  MeteoFile(Locality,int(day),int(month),duration,timestep)[5]  # Sun altitude
    sun_azimut      =  MeteoFile(Locality,int(day),int(month),duration,timestep)[6]  # Sun azimut
    sun_zenit       =  MeteoFile(Locality,int(day),int(month),duration,timestep)[7]  # Sun zenit
    latitude_deg    =  MeteoFile(Locality,int(day),int(month),duration,timestep)[8]  # Locality laitude
    longitude_deg   =  MeteoFile(Locality,int(day),int(month),duration,timestep)[9]  # Locality longitude
    start           =  MeteoFile(Locality,int(day),int(month),duration,timestep)[10] # Simulation start time
    UTC             =  MeteoFile(Locality,int(day),int(month),duration,timestep)[11] # Locality Coordinated Universal Time
    
    samples =  int(duration/timestep) # number of total steps
    
    """Calling 2006-today Single-Family House building class. """
    Archetype = SFH_90(ACH = 0.2, wall_absorption = wall_absorption, radiative_part = radiative_part, convective_part = convective_part)
    
    """Defining solar and internal gains. """
    southWindow  = Window_solar(azimut_surface=0.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window, area = Archetype.S_windows_area)   
    eastWindow  = Window_solar(azimut_surface=270.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window,  area = Archetype.E_windows_area) 
    westWindow  = Window_solar(azimut_surface=90.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window,  area = Archetype.W_windows_area)  
    northWindow  = Window_solar(azimut_surface=180.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window, area = Archetype.N_windows_area) 

    Archetype.Uvalues()
    southWall    = Wall_solar(azimut_surface=0.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0) , equivalent_transmittance = Archetype.Uwall, area = Archetype.south_area)
    eastWall     = Wall_solar(azimut_surface=270.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0),  equivalent_transmittance = Archetype.Uwall,area = Archetype.east_area)       
    westWall     = Wall_solar(azimut_surface=90.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0), equivalent_transmittance = Archetype.Uwall, area = Archetype.west_area)        
    northWall    = Wall_solar(azimut_surface=180.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0) , equivalent_transmittance = Archetype.Uwall , area = Archetype.north_area)    
    horRoof      = Wall_solar(azimut_surface=0.0, inclination_surface = 0.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0) , equivalent_transmittance = Archetype.Uroof , area = Archetype.roof_area)
    
    window_solar = []
    walls_solar  = []
    roof_solar   = []

    for h in range(0,len(sun_altitude)):     
        southWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        eastWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        westWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        northWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        southWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        eastWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        westWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        northWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        horRoof.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        """Solar gains to give in input to the state-space model formulation. """
        window_solar.append(southWindow.windows_solar_gains + eastWindow.windows_solar_gains + westWindow.windows_solar_gains + northWindow.windows_solar_gains)   # Solar gains
        walls_solar.append(southWall.walls_solar_gains + eastWall.walls_solar_gains + westWall.walls_solar_gains + northWall.walls_solar_gains)                    # Solar gains
        roof_solar.append(horRoof.walls_solar_gains)      

    floor_area = Archetype.floor_area

    """Internal gains due to people, lighting, and equipment are calculated from the occupancy profile defined using the richardsonpy tool
        (see the python file "User_pattern.py"). Then the trends are given as input to the archetype function. 
        However, it is possible to consider them constant during the simulation period calculated according to the methodology described in UNI 11300-1. """
    
    # To consider internal gains constant, then activate the following lines of code
    """
    Af = floor_area*2   
    if Af<=120: 
       internalGains = 7.987 * Af - 0.0353 * Af^2  # Internal gains (UNI 11300-1 pr. 13) 
    else: 
       internalGains = 450
    """
    
    """Defining the upper (Upper) and lower (Lower) ranges of the temperature profiles set in the thermostat. """
    Upper = []
    Lower = []

    if thermal_load == 'cooling':
        for i in range (samples):
            if occupancy[i] > 0:
                Upper.append(float(th_sp[i])  + th_tolerance_sup)
            else:
                Upper.append(30) # meaning the shutdown of the cooling system (active cooling)
            Lower.append(float(th_sp[i]) - th_tolerance_min)
    elif thermal_load == 'heating':
        for i in range (samples):
            if occupancy[i] > 0:
                Lower.append(float(th_sp[i]) - th_tolerance_min)
            else:
                Lower.append(float(th_sp[i]) - th_tolerance_min)
                #Lower.append(10) # meaning the shutdown of the heating system (active heating)
            Upper.append(float(th_sp[i])  + th_tolerance_sup)
    
    """Configuring the state-space model representation for linear programming problem resolution. """
    if thermal_load == 'cooling':
        Archetype.stateSpace_air_cooling() # state-space representation for air cooling system
    elif thermal_load == 'heating':
        Archetype.stateSpace_air_heating() # state-space representation for air heating system
    
    U = np.empty([6,samples]) # input vector of the state-space model formulation
   
    for r in range(0,len(sun_altitude)):       
        U[0,r] = T_est[r]                  # Outdoor temperature (°C)
        U[1,r] = T_ground[r]               # Ground temperature (°C)
        U[2,r] = window_solar[r]           # Solar Gains Windows (W) 
        U[3,r] = walls_solar[r]            # Solar Gains-Wall (W) 
        U[4,r] = roof_solar[r]             # Solar Gains-Roof (W)  
        U[5,r] = internalGains[r]          # Internal gains (W)    

    Xo     = np.array([[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]]]) # initial state of indoor air temperature

    Archetype.sys_dis_forOp(Archetype.sys_tuple,timestep) # System discretization for linear resolution
    
    Ad = np.asarray(Archetype.sys_disc_A) # State matrix A of the discretized system written in SSM form
    Bd = np.asarray(Archetype.sys_disc_B) # Input matrix B of the discretized system written in SSM form
    Ed = np.asarray(Archetype.sys_disc_E) # State matrix E of the discretized system isolating the decision variable 

    Z     = optForm_build(samples, Ad, Xo, Bd, U, Ed, 6)[0] # Variable describing the state of the system containing the thermal demand (decision variable)
    M_max = optForm_build(samples, Ad, Xo, Bd, U, Ed, 6)[1] # Variable describing the state of the indoor air temperture
    M     = optForm_build(samples, Ad, Xo, Bd, U, Ed, 6)[2]    
    
    """Defining linear programming problem. """
    Nc         = 2  # Number of constraints (2 for indoor air temperature)

    c         = [] # Decision variable coefficient 
    bounds    = [] # Decision variable bounds

    for r in range(0,samples):
        bounds.append((0,Q_hp[r]))
        c.append(timestep)

    Aub_opt = np.zeros([(Nc)*samples,samples]) # Matrix of inequality constraints
    Bub_opt = np.zeros([(Nc)*samples,1])       # Vector of inequality constraints
    
    """Compiling matrix and vector of the inequality constraints. 
        The constraint imposed in the optimization refers only to the state of indoor air temperature (Tair). """
    """Temperature constraints on the state of indoor air temperature system. """
    for i in range(0,samples):
        # Tair <= Upper
        Aub_opt[i,0:(samples)] = Z[i,:]
        Bub_opt[i,0]           = float(Upper[i])-M_max[i]      
    for i in range(0,samples): 
        # Tair >= Lower
        Aub_opt[i+(1*samples),0:(samples)] = Z[i,:]*(-1)
        Bub_opt[i+(1*samples),0]           = M_max[i]-float(Lower[i])

    Aub = Aub_opt.reshape(samples*(Nc),samples)
    Bub = Bub_opt.reshape(samples*(Nc),1)     

    return(c,      # 0 - decision variable coefficients
           Aub,    # 1 - Matrix of inequality constraints
           Bub,    # 2 - Vector of inequlity constraints
           bounds, # 3 - Decision variable bounds
           M,      # 4 - Variable describing the state of the indoor air temperture
           Z,      # 5 - Variable describing the state of the system containing the thermal demand (decision variable)
           M_max,  # 6 - Variable describing the state of the indoor air temperture
           Upper,  # 7 - Upper range of indoor temperature trend
           Lower,  # 8 - Lower range of indoor temperature trend
           floor_area) # 9 - Floor area (m2)

"""
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
1946-1960 Single Family House with air cooling/heating system
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
def SFH1960_air(Locality, day, month, duration, timestep, th_sp, th_tolerance_sup, th_tolerance_min, thermal_load, occupancy, internalGains, Q_hp):     
    """Recall of variables useful for defining inputs. """
    T_est           =  MeteoFile(Locality,int(day),int(month),duration,timestep)[0]  # Outdoor temperature trends (°C)
    T_ground        =  MeteoFile(Locality,int(day),int(month),duration,timestep)[1]  # Ground temperature trends (°C)
    Diff_Whm2       =  MeteoFile(Locality,int(day),int(month),duration,timestep)[2]  # Diffuse solar radiation (Whm2)
    Dir_Whm2        =  MeteoFile(Locality,int(day),int(month),duration,timestep)[3]  # Direct solar radiation (Whm2)
    Glo_Whm2        =  MeteoFile(Locality,int(day),int(month),duration,timestep)[4]  # Global solar radiation (Whm2)
    sun_altitude    =  MeteoFile(Locality,int(day),int(month),duration,timestep)[5]  # Sun altitude
    sun_azimut      =  MeteoFile(Locality,int(day),int(month),duration,timestep)[6]  # Sun azimut
    sun_zenit       =  MeteoFile(Locality,int(day),int(month),duration,timestep)[7]  # Sun zenit
    latitude_deg    =  MeteoFile(Locality,int(day),int(month),duration,timestep)[8]  # Locality laitude
    longitude_deg   =  MeteoFile(Locality,int(day),int(month),duration,timestep)[9]  # Locality longitude
    start           =  MeteoFile(Locality,int(day),int(month),duration,timestep)[10] # Simulation start time
    UTC             =  MeteoFile(Locality,int(day),int(month),duration,timestep)[11] # Locality Coordinated Universal Time
    
    samples =  int(duration/timestep) # number of total steps
    
    """Calling 2006-today Single-Family House building class. """
    Archetype = SFH_60(ACH = 0.2,wall_absorption = wall_absorption, radiative_part = radiative_part, convective_part = convective_part)

    """Defining solar and internal gains. """
    southWindow  = Window_solar(azimut_surface=0.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window, area = Archetype.S_windows_area)   
    eastWindow  = Window_solar(azimut_surface=270.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window,  area = Archetype.E_windows_area) 
    westWindow  = Window_solar(azimut_surface=90.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window,  area = Archetype.W_windows_area)  
    northWindow  = Window_solar(azimut_surface=180.0, inclination_surface = 90.0, glass_solar_transmittance = Archetype.g_window, area = Archetype.N_windows_area) 
    
    Archetype.Uvalues()
    southWall    = Wall_solar(azimut_surface=0.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0) , equivalent_transmittance = Archetype.Uwall, area = Archetype.south_area)
    eastWall     = Wall_solar(azimut_surface=270.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0),  equivalent_transmittance = Archetype.Uwall,area = Archetype.east_area)       
    westWall     = Wall_solar(azimut_surface=90.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0), equivalent_transmittance = Archetype.Uwall, area = Archetype.west_area)        
    northWall    = Wall_solar(azimut_surface=180.0, inclination_surface = 90.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0) , equivalent_transmittance = Archetype.Uwall , area = Archetype.north_area)    
    horRoof      = Wall_solar(azimut_surface=0.0, inclination_surface = 0.0, absorption= wall_absorption, surface_thermal_resistance = (1/25.0) , equivalent_transmittance = Archetype.Uroof , area = Archetype.roof_area)                            
    
    window_solar = []
    walls_solar  = []
    roof_solar   = []

    for h in range(0,len(sun_altitude)):     
        southWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        eastWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        westWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        northWindow.window_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        southWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        eastWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        westWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        northWall.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        horRoof.wall_solar_gains(sun_altitude = sun_altitude[h], sun_azimut = sun_azimut[h], sun_zenit = sun_zenit[h], dirnorrad_Whm2 = Dir_Whm2[h], difhorrad_Whm2 = Diff_Whm2[h], time = (start+(h*timestep)), latitude_deg = latitude_deg , longitude_deg = longitude_deg, UTC = UTC, glohorrad_Whm2 = Glo_Whm2[h], ground_reflectance = GR, timestep = timestep)
        """Solar gains to give in input to the state-space model formulation. """
        window_solar.append(southWindow.windows_solar_gains + eastWindow.windows_solar_gains + westWindow.windows_solar_gains + northWindow.windows_solar_gains)   # Solar gains
        walls_solar.append(southWall.walls_solar_gains + eastWall.walls_solar_gains + westWall.walls_solar_gains + northWall.walls_solar_gains)                    # Solar gains
        roof_solar.append(horRoof.walls_solar_gains)

    floor_area = Archetype.floor_area

    """Internal gains due to people, lighting, and equipment are calculated from the occupancy profile defined using the richardsonpy tool
        (see the python file "User_pattern.py"). Then the trends are given as input to the archetype function. 
        However, it is possible to consider them constant during the simulation period calculated according to the methodology described in UNI 11300-1. """
    
    # To consider internal gains constant, then activate the following lines of code
    """
    Af = floor_area*2   
    if Af<=120: 
       internalGains = 7.987 * Af - 0.0353 * Af^2  # Internal gains (UNI 11300-1 pr. 13) 
    else: 
       internalGains = 450
    """
        
    """Defining the upper (Upper) and lower (Lower) ranges of the temperature profiles set in the thermostat. """
    Upper = []
    Lower = []

    if thermal_load == 'cooling':
        for i in range (samples):
            if occupancy[i] > 0:
                Upper.append(float(th_sp[i])  + th_tolerance_sup)
            else:
                Upper.append(30) # meaning the shutdown of the cooling system (active cooling)
            Lower.append(float(th_sp[i]) - th_tolerance_min)
    elif thermal_load == 'heating':
        for i in range (samples):
            if occupancy[i] > 0:
                Lower.append(float(th_sp[i]) - th_tolerance_min)
            else:
                Lower.append(float(th_sp[i]) - th_tolerance_min)
                #Lower.append(10) # meaning the shutdown of the heating system (active heating)
            Upper.append(float(th_sp[i])  + th_tolerance_sup)
    
    """Configuring the state-space model representation for linear programming problem resolution. """
    if thermal_load == 'cooling':
        Archetype.stateSpace_air_cooling() # state-space representation for air cooling system
    elif thermal_load == 'heating':
        Archetype.stateSpace_air_heating() # state-space representation for air heating system
    
    U = np.empty([6,samples]) # input vector of the state-space model formulation
   
    for r in range(0,len(sun_altitude)):       
        U[0,r] = T_est[r]                  # Outdoor temperature (°C)
        U[1,r] = T_ground[r]               # Ground temperature (°C)
        U[2,r] = window_solar[r]           # Solar Gains Windows (W) 
        U[3,r] = walls_solar[r]            # Solar Gains-Wall (W) 
        U[4,r] = roof_solar[r]             # Solar Gains-Roof (W)  
        U[5,r] = internalGains[r]          # Internal gains (W)    

    Xo     = np.array([[th_sp[0]],[th_sp[0]],[th_sp[0]],[th_sp[0]]]) # initial state of indoor air temperature

    Archetype.sys_dis_forOp(Archetype.sys_tuple,timestep) # System discretization for linear resolution
    
    Ad = np.asarray(Archetype.sys_disc_A) # State matrix A of the discretized system written in SSM form
    Bd = np.asarray(Archetype.sys_disc_B) # Input matrix B of the discretized system written in SSM form
    Ed = np.asarray(Archetype.sys_disc_E) # State matrix E of the discretized system isolating the decision variable 

    Z     = optForm_build(samples, Ad, Xo, Bd, U, Ed, 3)[0] # Variable describing the state of the system containing the thermal demand (decision variable)
    M_max = optForm_build(samples, Ad, Xo, Bd, U, Ed, 3)[1] # Variable describing the state of the indoor air temperture
    M     = optForm_build(samples, Ad, Xo, Bd, U, Ed, 3)[2]    
    
    """Defining linear programming problem. """
    Nc         = 2  # Number of constraints (2 for indoor air temperature)

    c         = [] # Decision variable coefficient 
    bounds    = [] # Decision variable bounds

    for r in range(0,samples):
        bounds.append((0,Q_hp[r]))
        c.append(timestep)

    Aub_opt = np.zeros([(Nc)*samples,samples]) # Matrix of inequality constraints
    Bub_opt = np.zeros([(Nc)*samples,1])       # Vector of inequality constraints
    
    """Compiling matrix and vector of the inequality constraints. 
        The constraint imposed in the optimization refers only to the state of indoor air temperature (Tair). """
    """Temperature constraints on the state of indoor air temperature system. """
    for i in range(0,samples):
        # Tair <= Upper
        Aub_opt[i,0:(samples)] = Z[i,:]
        Bub_opt[i,0]           = float(Upper[i])-M_max[i]      
    for i in range(0,samples): 
        # Tair >= Lower
        Aub_opt[i+(1*samples),0:(samples)] = Z[i,:]*(-1)
        Bub_opt[i+(1*samples),0]           = M_max[i]-float(Lower[i])

    Aub = Aub_opt.reshape(samples*(Nc),samples)
    Bub = Bub_opt.reshape(samples*(Nc),1)     
       
    return(c,      # 0 - decision variable coefficients
           Aub,    # 1 - Matrix of inequality constraints
           Bub,    # 2 - Vector of inequlity constraints
           bounds, # 3 - Decision variable bounds
           M,      # 4 - Variable describing the state of the indoor air temperture
           Z,      # 5 - Variable describing the state of the system containing the thermal demand (decision variable)
           M_max,  # 6 - Variable describing the state of the indoor air temperture
           Upper,  # 7 - Upper range of indoor temperature trend
           Lower,  # 8 - Lower range of indoor temperature trend
           floor_area) # 9 - Floor area (m2)
