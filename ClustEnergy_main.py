"""
ClustEnergy OpTool
Released on January, 2024
@author: DIISM UNIVPM
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
#---------------------------------------------------------------------------------------------------------------------------------------
"""Check if the following python libraries needed to run the simulations are installed:
        pandas, numpy, scipy, pvlib, matplotlib, os, sys, psychrolib, math and datetime.
        
        🡲 To start, define the following simulation settings🡰 """
#---------------------------------------------------------------------------------------------------------------------------------------
""" SINULATION SETTINGS """
#---------------------------------------------------------------------------------------------------------------------------------------
"""Locality (epw file): Reference location (simulation of the outdoor environment).
        Available locations: "Ancona" (IT), "Messina" (IT), "Milano" (IT), "Roma" (IT) and "Denver" (USA).
        For additions, download from https://climate.onebuilding.org/ and add location within "MeteoF" fuction in MeteoFile.py"""

Locality            = "Roma"
#---------------------------------------------------------------------------------------------------------------------------------------
"""Thermal_load: Choose whether to consider "cooling" or "heating".
        By choosing "cooling" the thermal loads given as input to the space state modeling formulation will be with negative sign. 
        In addition, the performance characteristics of a heat pump and thermostat temperature set-point for cooling will be considered.
        In contrast, by choosing "heating" the thermal loads given as input will be with positive sign. 
        In addition, the performance characteristics of a heat pump and thermostat temperature set-points for heating will be considered."""

Thermal_load        = "heating"
#---------------------------------------------------------------------------------------------------------------------------------------
"""month: Simulation start month from 1 (January) to 12 (December)."""

month               = 2
#---------------------------------------------------------------------------------------------------------------------------------------
"""day: Simulation start day."""

day                 = 16
#---------------------------------------------------------------------------------------------------------------------------------------
"""n_days: Simulation days (max days of simulation = 365 - simulation start day of the year)."""

n_days              = 1 # 97° giorno = 7 aprile (fine riscaldamento), 319° giorno = 15 novembre (inizio riscaldamento)
#---------------------------------------------------------------------------------------------------------------------------------------
"""mins: Time step in minutes."""

mins                = 15         # minutes
#---------------------------------------------------------------------------------------------------------------------------------------
"""Selection of the number of buildings composing the cluster (example of max 20 buildings).
        Please note that the code is currently written to consider a maximum of 20 buildings, which is useful for performing simulations
        on a representative cluster without excessive run times. Considering more than 20 buildings requires the addition of additional 
        occupancy profiles, relative internal gains and temperature set-points beyond those predefined in 'User_pattern.py.' """
"""n_SFH06_air: Number of Single-Family Houses with class year of construction 2006-present equipped with air emission system.
        Thermal loads are given as input to the air node.
        U-values (W m-2 K-1) for archetype SFH06 according to Tabula Project (Corrado et al., 2014)
        SFH        	Ext_walls	Floor	Roof	Windows
        2006-today 	0.34     	0.33 	0.28	2.20   """
n_SFH06_air         = 1

"""n_SFH05_air: Number of Single-Family Houses with class year of construction 1991-2005 equipped with air cooling/heating system.
        Thermal loads are given as input to the air node.
        U-values (W m-2 K-1) for archetype SFH05 according to Tabula Project (Corrado et al., 2014)
        SFH        	Ext_walls	Floor	Roof	Windows
        1991-2005  	0.59     	0.63 	0.57	2.40   """
n_SFH05_air         = 1

"""n_SFH90_air: Number of Single-Family Houses with class year of construction 1976-1990 equipped with air cooling/heating system.
        Thermal loads are given as input to the air node.
        U-values (W m-2 K-1) for archetype SFH90 according to Tabula Project (Corrado et al., 2014)
        SFH        	Ext_walls	Floor	Roof	Windows
        1976-1990  	0.76     	0.76 	1.14	2.80   """
n_SFH90_air         = 1

"""n_SFH60_air: Number of Single-Family Houses with class year of construction 1946-1960 equipped with air cooling/heating system.
        Thermal loads are given as input to the air node.
        U-values (W m-2 K-1) for archetype SFH60 according to Tabula Project (Corrado et al., 2014)
        SFH        	Ext_walls	Floor	Roof	Windows
        1946-1960  	1.48     	2.00 	2.20	4.90   """
n_SFH60_air         = 1

"""n_SFH06_floor: Number of Single-Family Houses with class year of construction 2006-present equipped with underfloor heating system.
        Thermal loads are given as input to the floor node.
        U-values (W m-2 K-1) for archetype SFH06 according to Tabula Project (Corrado et al., 2014)
        SFH        	Ext_walls	Floor	Roof	Windows
        2006-today 	0.34     	0.33 	0.28	2.20   """
n_SFH06_floor       = 1 # only heating

ACH = 0.5                # Air changes per hour

GR  = 0.2                # ground reflectance
wall_absorption   = 0.3  # wall absorption

#---------------------------------------------------------------------------------------------------------------------------------------
"""Internal gains."""
int_gains_profile = "fixed" # or variable (occupancy dependent)

radiative_part  = 0.6    # radiative part of internal gains  
convective_part = 0.4    # convective part of internal gains  
#---------------------------------------------------------------------------------------------------------------------------------------
"""Set the heat pumps characteristics for space cooling or heating.
        With "pre_defined_HP" in the following "HP_data_type" it is possible to choose between some commercial heat pumps provided in the library.
        While, with "user_defined_HP" it is possibile to consider normalized performance characteristic curves provided by the user via external Excel file 
        (see the documentation for further information). """

HP_data_type = "pre_defined_HP" 

"""Set the supply temperatures, rated capacities and COP of heat pumps for space cooling.
        These input operating points refer to a commercial air-water heat pump. 
        Once defined, as the outdoor air temperature changes, the Coefficient of Performance (COP) 
        and the maximum thermal power (Qmax) will be output, based on the normalized performance characteristic curves provided 
        by the manufacturer at capacity ratio (CR) = 1 (if pre_defined_HP data) or provided by the user (if user_provided_HP data). """

Supply_air    = 40.0 # °C - Supply temperature for air systems
Supply_floor  = 35.0 # °C - Supply temperature for underfloor heating system

"""rated_cap_n_air/floor (with n for each archetype of SFHs): Rated heating or cooling capacity of heat pumps.
        If pre_defined_HP data select between the following heat pumps available in the library in reference to the data provided by the manufacturer (Viessmann, 2020).
        
        HEATING (with reference to supply temperature of 35°C and outdoor air temperature of 7°C)
        - rated_capacity: 8.64 and rated_COP: 4.19
        - rated_capacity: 14.00 and rated_COP: 4.08
        - rated_capacity: 18.00 and rated_COP: 3.90
        
        COOLING (with reference to supply temperature of 18°C and outdoor air temperature of 35°C)
        - rated_capacity: 10.00 and rated_COP: 2.70
        - rated_capacity: 11.50 and rated_COP: 2.94
        - rated_capacity: 14.22 and rated_COP: 3.43
        
        REFERENCES:
        Viessmann, 2020. “Dati integrativi pompe di calore VITOCAL 250-S per il calcolo delle prestazioni energetiche degli edifici, secondo UNI/TS 11300 parte 4”.
        """
        
"""Air system."""
# SFH 2006-today
rated_cap_06_air     = 8.64  # kW
rated_COP_06_air     = 4.19

# SFH 1991-2005
rated_cap_05_air     = 14.00 # kW
rated_COP_05_air     = 4.08

# SFH 1976-1990
rated_cap_90_air     = 14.00 # kW
rated_COP_90_air     = 4.08

# SFH 1946-1960
rated_cap_60_air     = 18.00 # kW
rated_COP_60_air     = 3.90

"""Underfloor heating system."""
# SFH 2006-today
rated_cap_06_floor   = 14.00  # kW
rated_COP_06_floor   = 4.08

#---------------------------------------------------------------------------------------------------------------------------------------
"""Thermostat settings."""
"""tsp_profile: Choose whether to consider a "pre_defined" daily set-point temperature profile or to set a "user_defined". 
        If "pre_defined", see in the documentation how to define new temperature setpoint profiles."""
tsp_profile               = "pre_defined"

"""user_def_tsp: if "user_defined" profile set the thermostat set-point temperature (e.g., 20°C for heating or 26°C for cooling)."""
user_def_tsp              = 26.0              # °C
#---------------------------------------------------------------------------------------------------------------------------------------
"""Thermostat tolerances (baseline scenario)."""
"Set thermostat tolerances during baseline (BL) and Demand-Response (DR) scenarios."
"""th_tolerance_BL_i: Upper/lower thermostat tolerance during baseline (BL) scenario for air cooling/heating and underfloor heating systems.
        As a result, the indoor air temperature during the BL scenario will range between 
        upper_BL = tsp_profile + th_tolerance_BL_up and lower_BL = tsp_profile - th_tolerance_BL_low."""
th_tolerance_BL_up_air    = 2.0    # °C - BL upper tolerance for air cooling/heating systems
th_tolerance_BL_low_air   = 0.0    # °C - BL lower tolerance for air cooling/heating systems

th_tolerance_BL_up_floor  = 3.0    # °C - BL upper tolerance for underfloor heating systems
th_tolerance_BL_low_floor = 2.0    # °C - BL lower tolerance for underfloor heating systems


#-----------------------------------------------------DOMESTIC HOT WATER----------------------------------------------------------------
"""To consider domestic hot water production, type "True" in DHW_on variable and set the following parameters."""
DHW_on         = True
T_inlet        = 10        # °C - temperature of the water entering the tank

tank_temp      = 60        # °C - water set-point temperature
tank_tol       = 5.0       # °C - temperature tolerance
tank_volume    = 178       # m3 - water volume of the tank
h_tank         = 1.555     # height - m
l_tank         = 0.586     # length - m
U_tank         = 0.23      # W/(m2*K)

# HP for DHW production
rated_COP_tank = 3.1
rated_cap_tank = 1.3

Q_draw         = "variable"   # DHW draw profile (fixed or variable (occupancy dependent))
                              # To create new occupancy dependent DHW draw profile use the module "DHWdraw_profile_calc"
#-----------------------------------------------DEMAND RESPONSE STRATEGY----------------------------------------------------------------
#
"""Demand_response: Choose the type of Demand Response scenario for managing aggregate loads.
        Available management strategies:
            - "base_load" (baseline scenario without comparing to a DR scenario);
            - "peak_shaving" (reduction of peak loads); 
            - "pv_centralized" (electric demand shifting under centralized PV generation, e.g. shared energy resources);
            - "pv_distributed" (electric demand shifting under distributed PV generation, e.g. without shared energy resources).
            - "load_shifting" (electric demand shifting under electricity price signal).
              "pv_central_price" (electric demand shifting under electricity price signal and centralized PV genration)
            """

Demand_response     = "pv_central_price"

"""Define the limit of aggregate electricity consumption relative to the baseline scenario in order to limit rebound effects. 
        Accordingly, Eele_max = Eele_BL*f_limit. """
f_limit = 1000 # consumption limit factor compared with the baseline profile (low values limit rebound effects)
#---------------------------------------------------------------------------------------------------------------------------------------
"""Peak shaving settings."""
"""f_red: Peak load reduction factor during peak_shaving scenario (from 0 (no reduction) to 1 (100% reduction))."""
f_red               = 0.2       # dimensionless

"""duration_dr: Duration of the peak load reduction event (to be considered from the time step in which the maximum electrical demand occurs)."""
duration_dr         = 1         # hours
#---------------------------------------------------------------------------------------------------------------------------------------
"""Photovoltaic generation settings."""
"""pv_size: Rated power of the photovoltaic systems available to archetypes "calculated" according to the Italian regulation Dlgs. 199/2021 or "user_defined"."""
pv_size             = "calculated"

"""panel_dim                If "calculated", specify the rated power in W of the PV panel."""
panel_dim           = 396.0 # Watt - PV panel

"""pv_rp_n_air: If "user_defined", specify the rated power of th PV system in kW for archetypes with air cooling/heating system."""
pv_rp_06_air        = 3.0   # kW - PV system 
pv_rp_05_air        = 3.0   # kW - PV system 
pv_rp_90_air        = 3.0   # kW - PV system 
pv_rp_60_air        = 3.0   # kW - PV system 

"""pv_rp_n_floor: If "user_defined", specify the rated power of the PV system in kW for archetype with underfloor heating system."""
pv_rp_06_floor      = 3.0   # kW - PV system 

#---------------------------------------------------------------------------------------------------------------------------------------
"""Electricity price settings."""
"""A Time Of Use rate is considered as a price signal during a load_shifting strategy.
        Specifically, below in price_ranges can be defined the time bands and the electricity price (Eur/kWh)."""
# A day-night (bi-hourly) electricity tariff is set below.
price_ranges = {"hour_day": 8.0, "price_day": 0.0286,      # Define the start time and electricity price (Eur/kWh) during day hours
                "hour_night": 23.0, "price_night": 0.011} # Define the start time and electricity price (Eur/kWh) during night hours

#---------------------------------------------------------------------------------------------------------------------------------------
"""Thermostat tolerances (Demand-Response scenario)."""
"""th_tolerance_DR_i: Upper/lower thermostat tolerances during Demand-Response (DR) scenario to unlock energy flexibility of buildings.
        As a result, the indoor air temperature during the BL scenario will range between 
        upper_DR = upper_BL + th_tolerance_DR_up and lower_DR = lower_BL - th_tolerance_DR_low.
        During the peak load reduction scenario for space cooling, increase the upper_BL by th_tolerance_DR_up.
        Similarly, during the peak load reduction scenario for space heating, reduce the lower_BL by th_tolerance_DR_low.
        While, during load shifting scenario for space cooling under PV generation, reduce lower_BL by th_tolerance_DR_low.
        Similarly, during load shifting scenario for space heating under PV generation, increase upper_BL by th_tolerance_DR_up."""
th_tolerance_DR_up_air    = 0.5    # °C - DR upper tolerance for air cooling/heating systems
th_tolerance_DR_low_air   = 0.5    # °C - DR lower tolerance for air cooling/heating systems

th_tolerance_DR_up_floor  = 0.0    # °C - DR upper tolerance for underfloor heating systems
th_tolerance_DR_low_floor = 0.0    # °C - DR lower tolerance for underfloor heating systems
#---------------------------------------------------------------------------------------------------------------------------------------
"""Plot settings."""
"""Decide whether to plot aggregated electric demand trends (True) or not (False), for comparison between BL and DR scenarios."""
do_plot_Ebl_cluster = True
do_plot_Edr_cluster = True

"""Decide whether to plot buildings electrical demand trends (True) or not (False), for comparison between BL and DR scenarios."""
do_plot_E_single  =  True

"""Decide whether to plot buildings indoor air temperatures trends (True) or not (False), for comparison between BL and DR scenarios."""
do_plot_Tair_single  =  True

"""Choose excel file name for exporting data regarding electricity demand of individual buildings."""
name_file = "TOU_shift"

#---------------------------------------------------------------------------------------------------------------------------------------
"""Import libraries"""
from datetime import datetime
from MeteoFile import CheckMeteoFile
import strategies

""" SINULATION START """
start_calc_sim = datetime.now()

CheckMeteoFile(Locality)                           # Check the locality
duration = n_days*24                               # Calculating simulation duration - hours
timestep     = mins/60                             # Timestep conversion - hours
samples = int(duration/timestep)                   # Calculation of total time steps
time = [i*timestep for i in range(0, samples)]     # Time array
# --------------------------------------------------------------------------------------------------------------------------------------
if Demand_response == "base_load":
    """If Demand_response = "base_load", call 'base_load' function from module 'strategies.py'."""
    getattr(strategies,Demand_response)(Locality, day, month, duration, timestep, Thermal_load, 
                                n_SFH06_air, n_SFH05_air, n_SFH90_air, n_SFH60_air,n_SFH06_floor, 
                                ACH, GR, wall_absorption, int_gains_profile, radiative_part, convective_part,
                                HP_data_type, Supply_air, rated_cap_06_air, rated_cap_05_air, rated_cap_90_air, rated_cap_60_air, rated_COP_06_air, rated_COP_05_air, rated_COP_90_air, rated_COP_60_air,
                                Supply_floor, rated_cap_06_floor, rated_COP_06_floor,
                                tsp_profile, user_def_tsp, th_tolerance_BL_up_air, th_tolerance_BL_low_air,th_tolerance_BL_up_floor, th_tolerance_BL_low_floor,int_gains_profile,
                                DHW_on,T_inlet,tank_temp,tank_tol,tank_volume,h_tank,l_tank,U_tank,rated_COP_tank,rated_cap_tank,Q_draw,
                                do_plot_E_single, do_plot_Ebl_cluster, do_plot_Tair_single, name_file)
    
elif Demand_response == "peak_shaving":
    """If Demand_response = "peak_shaving", call 'peak_shaving' function from module 'strategies.py'."""
    getattr(strategies,Demand_response)(Locality, day, month, duration, timestep, Thermal_load, 
                                n_SFH06_air, n_SFH05_air, n_SFH90_air, n_SFH60_air,n_SFH06_floor, 
                                ACH, GR, wall_absorption, int_gains_profile, radiative_part, convective_part,
                                HP_data_type, Supply_air, rated_cap_06_air, rated_cap_05_air, rated_cap_90_air, rated_cap_60_air, rated_COP_06_air, rated_COP_05_air, rated_COP_90_air, rated_COP_60_air,
                                Supply_floor, rated_cap_06_floor, rated_COP_06_floor,
                                tsp_profile, user_def_tsp, th_tolerance_BL_up_air, th_tolerance_BL_low_air, th_tolerance_BL_up_floor, th_tolerance_BL_low_floor,
                                DHW_on,T_inlet,tank_temp,tank_tol,tank_volume,h_tank,l_tank,U_tank,rated_COP_tank,rated_cap_tank,Q_draw,
                                th_tolerance_DR_up_air, th_tolerance_DR_low_air, th_tolerance_DR_up_floor, th_tolerance_DR_low_floor, 
                                f_red, duration_dr, f_limit,
                                do_plot_E_single, do_plot_Edr_cluster, do_plot_Tair_single, name_file)
    
elif Demand_response == "pv_centralized":
    """If Demand_response = "pv_centralized", call 'pv_centralized' function from module 'strategies.py'."""
    getattr(strategies,Demand_response)(Locality, day, month, duration, timestep, Thermal_load, 
                                n_SFH06_air, n_SFH05_air, n_SFH90_air, n_SFH60_air,n_SFH06_floor, 
                                ACH, GR, wall_absorption, int_gains_profile, radiative_part, convective_part,
                                HP_data_type, Supply_air, rated_cap_06_air, rated_cap_05_air, rated_cap_90_air, rated_cap_60_air, rated_COP_06_air, rated_COP_05_air, rated_COP_90_air, rated_COP_60_air,
                                Supply_floor, rated_cap_06_floor, rated_COP_06_floor,
                                tsp_profile, user_def_tsp, th_tolerance_BL_up_air, th_tolerance_BL_low_air, th_tolerance_BL_up_floor, th_tolerance_BL_low_floor,
                                DHW_on,T_inlet,tank_temp,tank_tol,tank_volume,h_tank,l_tank,U_tank,rated_COP_tank,rated_cap_tank,Q_draw,
                                th_tolerance_DR_up_air, th_tolerance_DR_low_air, th_tolerance_DR_up_floor, th_tolerance_DR_low_floor,
                                pv_size, panel_dim, pv_rp_06_air, pv_rp_05_air, pv_rp_90_air, pv_rp_60_air, pv_rp_06_floor, f_limit,
                                do_plot_E_single, do_plot_Edr_cluster,do_plot_Tair_single, name_file)
    
elif Demand_response == "pv_distributed":
    """If Demand_response = "pv_distributed", call 'pv_distributed' function from module 'strategies.py'."""
    getattr(strategies,Demand_response)(Locality, day, month, duration, timestep, Thermal_load, 
                                n_SFH06_air, n_SFH05_air, n_SFH90_air, n_SFH60_air,n_SFH06_floor, 
                                ACH, GR, wall_absorption, int_gains_profile, radiative_part, convective_part,
                                HP_data_type, Supply_air, rated_cap_06_air, rated_cap_05_air, rated_cap_90_air, rated_cap_60_air, rated_COP_06_air, rated_COP_05_air, rated_COP_90_air, rated_COP_60_air,
                                Supply_floor, rated_cap_06_floor, rated_COP_06_floor,
                                tsp_profile, user_def_tsp, th_tolerance_BL_up_air, th_tolerance_BL_low_air, th_tolerance_BL_up_floor, th_tolerance_BL_low_floor,
                                DHW_on,T_inlet,tank_temp,tank_tol,tank_volume,h_tank,l_tank,U_tank,rated_COP_tank,rated_cap_tank,Q_draw,
                                th_tolerance_DR_up_air, th_tolerance_DR_low_air, th_tolerance_DR_up_floor, th_tolerance_DR_low_floor,
                                pv_size, panel_dim, pv_rp_06_air, pv_rp_05_air, pv_rp_90_air, pv_rp_60_air, pv_rp_06_floor, f_limit,
                                do_plot_E_single,do_plot_Edr_cluster,do_plot_Tair_single, name_file)

elif Demand_response == "load_shifting":
    """If Demand_response = "pv_distributed", call 'pv_distributed' function from module 'strategies.py'."""
    var=getattr(strategies,Demand_response)(Locality, day, month, duration, timestep, Thermal_load, 
                                n_SFH06_air, n_SFH05_air, n_SFH90_air, n_SFH60_air,n_SFH06_floor, 
                                ACH, GR, wall_absorption, int_gains_profile, radiative_part, convective_part,
                                HP_data_type, Supply_air, rated_cap_06_air, rated_cap_05_air, rated_cap_90_air, rated_cap_60_air, rated_COP_06_air, rated_COP_05_air, rated_COP_90_air, rated_COP_60_air,
                                Supply_floor, rated_cap_06_floor, rated_COP_06_floor,
                                tsp_profile, user_def_tsp, th_tolerance_BL_up_air, th_tolerance_BL_low_air, th_tolerance_BL_up_floor, th_tolerance_BL_low_floor,
                                th_tolerance_DR_up_air, th_tolerance_DR_low_air, th_tolerance_DR_up_floor, th_tolerance_DR_low_floor,
                                price_ranges, f_limit,
                                DHW_on,T_inlet,tank_temp,tank_tol,tank_volume,h_tank,l_tank,U_tank,rated_COP_tank,rated_cap_tank,Q_draw,
                                do_plot_E_single, do_plot_Edr_cluster,do_plot_Tair_single, name_file)

elif Demand_response == "pv_central_price":
    """If Demand_response = "pv_distributed", call 'pv_distributed' function from module 'strategies.py'."""
    var=getattr(strategies,Demand_response)(Locality, day, month, duration, timestep, Thermal_load, 
                                n_SFH06_air, n_SFH05_air, n_SFH90_air, n_SFH60_air,n_SFH06_floor, 
                                ACH, GR, wall_absorption, int_gains_profile, radiative_part, convective_part,
                                HP_data_type, Supply_air, rated_cap_06_air, rated_cap_05_air, rated_cap_90_air, rated_cap_60_air, rated_COP_06_air, rated_COP_05_air, rated_COP_90_air, rated_COP_60_air,
                                Supply_floor, rated_cap_06_floor, rated_COP_06_floor,
                                tsp_profile, user_def_tsp, th_tolerance_BL_up_air, th_tolerance_BL_low_air, th_tolerance_BL_up_floor, th_tolerance_BL_low_floor,
                                th_tolerance_DR_up_air, th_tolerance_DR_low_air, th_tolerance_DR_up_floor, th_tolerance_DR_low_floor,
                                price_ranges, f_limit,
                                pv_size, panel_dim, pv_rp_06_air, pv_rp_05_air, pv_rp_90_air, pv_rp_60_air, pv_rp_06_floor, 
                                DHW_on,T_inlet,tank_temp,tank_tol,tank_volume,h_tank,l_tank,U_tank,rated_COP_tank,rated_cap_tank,Q_draw,
                                do_plot_E_single, do_plot_Edr_cluster,do_plot_Tair_single, name_file)
else:
    print('Demand Response strategy not available')


stop_calc_sim = datetime.now()-start_calc_sim
print("program run time:",stop_calc_sim)

""" SIMULATION STOP """