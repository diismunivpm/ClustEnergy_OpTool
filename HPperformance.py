"""
ClustEnergy OpTool version 1.0.0
Released on February 29, 2024
@author: DIISM UNIVPM
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""

"""The capacity and coefficient of performance of heat pumps are calculated here as the outdoor temperature varies.
    With reference to commercial heat pumps, the model is based on catalog data for capacity (total capacity in heating/cooling mode)
    and coefficient of performance (COP) provided by the manufacturer (Viessmann, 2020).
    
    REFERENCES:
    Viessmann, 2020. “Dati integrativi pompe di calore VITOCAL 250-S per il calcolo delle prestazioni energetiche degli edifici, secondo UNI/TS 11300 parte 4”.
    """
from scipy.interpolate import interp2d
import numpy as np
import os 
import pandas as pd


def HP_cooling(Toutdoor,Tsupply,rated_cap,rated_COP,data_type): 

    """Heat pumps operating points with normalized data at a 7 °C outdoor temperature and 35 °C supply temperature. """
    if data_type == "pre_defined_HP":
        if rated_cap == 10.00 and rated_COP == 2.70:
            HP_perf = {"supply_temperature": [18,7], 
                      "dry_bulb_temperature": [35,30,27,25,20],
                      "Q_fraction": np.array([(1.0000,           # Capacity fraction at T_supply = 18 °C, Dry_bulb = 35 °C and partial load = 1
                                               1.0680,           # Capacity fraction at T_supply = 18 °C, Dry_bulb = 30 °C and partial load = 1
                                               1.1140,           # Capacity fraction at T_supply = 18 °C, Dry_bulb = 27 °C and partial load = 1
                                               1.1250,           # Capacity fraction at T_supply = 18 °C, Dry_bulb = 25 °C and partial load = 1
                                               1.1390),          # Capacity fraction at T_supply = 18 °C, Dry_bulb = 20 °C and partial load = 1
                                              (0.8000,           # Capacity fraction at T_supply = 7 °C, Dry_bulb = 35 °C and partial load = 1
                                               0.8770,           # Capacity fraction at T_supply = 7 °C, Dry_bulb = 30 °C and partial load = 1
                                               0.9290,           # Capacity fraction at T_supply = 7 °C, Dry_bulb = 27 °C and partial load = 1
                                               0.9680,           # Capacity fraction at T_supply = 7 °C, Dry_bulb = 25 °C and partial load = 1
                                               1.1390)])         # Capacity fraction at T_supply = 7 °C, Dry_bulb = 20 °C and partial load = 1
                                      , 
                      "COP_fraction": np.array([(1.0000,         # COP fraction at T_supply = 18 °C, Dry_bulb = 35 °C and partial load = 1
                                                 1.1926,         # COP fraction at T_supply = 18 °C, Dry_bulb = 30 °C and partial load = 1
                                                 1.3111,         # COP fraction at T_supply = 18 °C, Dry_bulb = 27 °C and partial load = 1
                                                 1.3889,         # COP fraction at T_supply = 18 °C, Dry_bulb = 25 °C and partial load = 1
                                                 1.5926),        # COP fraction at T_supply = 18 °C, Dry_bulb = 20 °C and partial load = 1
                                                (0.7778,         # COP fraction at T_supply = 7 °C, Dry_bulb = 35 °C and partial load = 1
                                                 0.9259,         # COP fraction at T_supply = 7 °C, Dry_bulb = 30 °C and partial load = 1
                                                 1.0407,         # COP fraction at T_supply = 7 °C, Dry_bulb = 27 °C and partial load = 1
                                                 1.1148,         # COP fraction at T_supply = 7 °C, Dry_bulb = 25 °C and partial load = 1
                                                 1.3000)])}      # COP fraction at T_supply = 7 °C, Dry_bulb = 20 °C and partial load = 1
        elif rated_cap == 11.50 and rated_COP == 2.94:
            HP_perf = {"supply_temperature": [18,7], 
                      "dry_bulb_temperature": [35,30,27,25,20],
                      "Q_fraction": np.array([(1.0000,           # Capacity fraction at T_supply = 18 °C, Dry_bulb = 35 °C and partial load = 1
                                               1.0661,           # Capacity fraction at T_supply = 18 °C, Dry_bulb = 30 °C and partial load = 1
                                               1.1209,           # Capacity fraction at T_supply = 18 °C, Dry_bulb = 27 °C and partial load = 1
                                               1.1583,           # Capacity fraction at T_supply = 18 °C, Dry_bulb = 25 °C and partial load = 1
                                               1.1800),          # Capacity fraction at T_supply = 18 °C, Dry_bulb = 20 °C and partial load = 1
                                              (0.7391,           # Capacity fraction at T_supply = 7 °C, Dry_bulb = 35 °C and partial load = 1
                                               0.7409,           # Capacity fraction at T_supply = 7 °C, Dry_bulb = 30 °C and partial load = 1
                                               0.7957,           # Capacity fraction at T_supply = 7 °C, Dry_bulb = 27 °C and partial load = 1
                                               0.7939,           # Capacity fraction at T_supply = 7 °C, Dry_bulb = 25 °C and partial load = 1
                                               0.8383)])         # Capacity fraction at T_supply = 7 °C, Dry_bulb = 20 °C and partial load = 1
                                      , 
                      "COP_fraction": np.array([(1.0000,         # COP fraction at T_supply = 18 °C, Dry_bulb = 35 °C and partial load = 1
                                                 1.2109,         # COP fraction at T_supply = 18 °C, Dry_bulb = 30 °C and partial load = 1
                                                 1.3639,         # COP fraction at T_supply = 18 °C, Dry_bulb = 27 °C and partial load = 1
                                                 1.4796,         # COP fraction at T_supply = 18 °C, Dry_bulb = 25 °C and partial load = 1
                                                 1.6429),        # COP fraction at T_supply = 18 °C, Dry_bulb = 20 °C and partial load = 1
                                                (0.8061,         # COP fraction at T_supply = 7 °C, Dry_bulb = 35 °C and partial load = 1
                                                 0.9864,         # COP fraction at T_supply = 7 °C, Dry_bulb = 30 °C and partial load = 1
                                                 1.1122,         # COP fraction at T_supply = 7 °C, Dry_bulb = 27 °C and partial load = 1
                                                 1.2041,         # COP fraction at T_supply = 7 °C, Dry_bulb = 25 °C and partial load = 1
                                                 1.3367)])}      # COP fraction at T_supply = 7 °C, Dry_bulb = 20 °C and partial load = 1
        elif rated_cap == 14.22 and rated_COP == 3.43:
            HP_perf = {"supply_temperature": [18,7], 
                      "dry_bulb_temperature": [35,30,27,25,20],
                      "Q_fraction": np.array([(1.0000,           # Capacity fraction at T_supply = 18 °C, Dry_bulb = 35 °C and partial load = 1
                                               1.0162,           # Capacity fraction at T_supply = 18 °C, Dry_bulb = 30 °C and partial load = 1
                                               1.0316,           # Capacity fraction at T_supply = 18 °C, Dry_bulb = 27 °C and partial load = 1
                                               1.0556,           # Capacity fraction at T_supply = 18 °C, Dry_bulb = 25 °C and partial load = 1
                                               1.0949),          # Capacity fraction at T_supply = 18 °C, Dry_bulb = 20 °C and partial load = 1
                                              (0.8340,           # Capacity fraction at T_supply = 7 °C, Dry_bulb = 35 °C and partial load = 1
                                               0.8608,           # Capacity fraction at T_supply = 7 °C, Dry_bulb = 30 °C and partial load = 1
                                               0.9248,           # Capacity fraction at T_supply = 7 °C, Dry_bulb = 27 °C and partial load = 1
                                               0.9613,           # Capacity fraction at T_supply = 7 °C, Dry_bulb = 25 °C and partial load = 1
                                               0.9979)])         # Capacity fraction at T_supply = 7 °C, Dry_bulb = 20 °C and partial load = 1
                                      , 
                      "COP_fraction": np.array([(1.0000,         # COP fraction at T_supply = 18 °C, Dry_bulb = 35 °C and partial load = 1
                                                 1.2128,         # COP fraction at T_supply = 18 °C, Dry_bulb = 30 °C and partial load = 1
                                                 1.3061,         # COP fraction at T_supply = 18 °C, Dry_bulb = 27 °C and partial load = 1
                                                 1.3878,         # COP fraction at T_supply = 18 °C, Dry_bulb = 25 °C and partial load = 1
                                                 1.5773),        # COP fraction at T_supply = 18 °C, Dry_bulb = 20 °C and partial load = 1
                                                (0.5423,         # COP fraction at T_supply = 7 °C, Dry_bulb = 35 °C and partial load = 1
                                                 0.5948,         # COP fraction at T_supply = 7 °C, Dry_bulb = 30 °C and partial load = 1
                                                 0.7143,         # COP fraction at T_supply = 7 °C, Dry_bulb = 27 °C and partial load = 1
                                                 0.7405,         # COP fraction at T_supply = 7 °C, Dry_bulb = 25 °C and partial load = 1
                                                 0.7930)])}      # COP fraction at T_supply = 7 °C, Dry_bulb = 20 °C and partial load = 1
        else:
            print("\nHeat pump of " + str(round(rated_cap,2)) + "kW rated capacity and" + str(round(rated_COP,2)) + " rated COP is not available in the library."
                  + "\nCheck available sizes or follow the instructions in the documentation to consider a user provided heat pump.")
    elif data_type == "user_provided_HP":
        # User provided heat pump data in Excel format
        
        """Reading the excel file. """
        currentPath     = os.path.dirname("__file__")                
        userPath     = os.path.join(currentPath, 'userProfiles')
        HPper_FileName = "HP_performance_cooling.xlsx"
        HPper_File     = os.path.join(userPath, HPper_FileName)
        
        HPper_data = dict(pd.read_excel(HPper_File)) # excel data to dictionary
        
        """Reading data. """
        supply       = []
        Test         = []
        Qth_fraction = []
        Cop_fraction = []
        
        for value in HPper_data["supply_temperature"]:
            if value > 0:
                supply.append(value)
                
        for value in HPper_data["dry_bulb_temperature"]:
            Test.append(value)
            
        for i in range(len(supply)):
            for value in HPper_data["Q_fraction_at_{0}".format(supply[i])]:
                Qth_fraction.append(value)
            for value in HPper_data["COP_fraction_at_{0}".format(supply[i])]:
                Cop_fraction.append(value)
        
        """Collecting data to be given as input in interpolation functions. """
        HP_perf      = {}
        HP_perf["supply_temperature"] = supply
        HP_perf["dry_bulb_temperature"] = Test
        HP_perf["Q_fraction"] = np.array([Qth_fraction[i*(len(Test)):(i+1)*(len(Test))] for i in range(len(supply))])
        HP_perf["COP_fraction"] = np.array([Cop_fraction[i*(len(Test)):(i+1)*(len(Test))] for i in range(len(supply))])
        
    """Creating the interpolation functions."""
    f_Qmax = interp2d(x=HP_perf["dry_bulb_temperature"], y=HP_perf["supply_temperature"], z=HP_perf["Q_fraction"])
    f_COP = interp2d(x=HP_perf["dry_bulb_temperature"], y=HP_perf["supply_temperature"], z=HP_perf["COP_fraction"])
    
    """Calculating Qmax and COP at the desired outdoor and supply temperatures. """
    Qmax   = float(f_Qmax(Toutdoor,Tsupply)*rated_cap*1000)
    COP   = float(f_COP(Toutdoor,Tsupply)*rated_COP)
    
    return (COP,   # 0 - Cefficient Of Performance
            Qmax)  # 1 - Total capacity (W)


def HP_heating(Toutdoor,Tsupply,rated_cap,rated_COP,data_type): 

    """Heat pumps operating points with normalized data at a 7 °C outdoor temperature and 35 °C supply temperature. """
    if data_type == "pre_defined_HP":
        if rated_cap == 8.64 and rated_COP == 4.19:
            HP_perf = {"supply_temperature": [35,45,55], 
                      "dry_bulb_temperature": [-7,-2,2,7,12],
                      "Q_fraction": np.array([(0.7639, # Capacity fraction at T_supply = 35 °C, Dry_bulb = -7 °C and partial load = 1
                                               0.8287,           # Capacity fraction at T_supply = 35 °C, Dry_bulb = -2 °C and partial load = 1
                                               0.8912,           # Capacity fraction at T_supply = 35 °C, Dry_bulb = 2 °C and partial load = 1
                                               1.0000,           # Capacity fraction at T_supply = 35 °C, Dry_bulb = 7 °C and partial load = 1
                                               1.1551),          # Capacity fraction at T_supply = 35 °C, Dry_bulb = 12 °C and partial load = 1
                                              (0.6840,          # Capacity fraction at T_supply = 45 °C, Dry_bulb = -7 °C and partial load = 1
                                               0.7419,           # Capacity fraction at T_supply = 45 °C, Dry_bulb = -2 °C and partial load = 1
                                               0.7975,           # Capacity fraction at T_supply = 45 °C, Dry_bulb = 2 °C and partial load = 1
                                               0.9838,           # Capacity fraction at T_supply = 45 °C, Dry_bulb = 7 °C and partial load = 1
                                               1.0336),          # Capacity fraction at T_supply = 45 °C, Dry_bulb = 12 °C and partial load = 1
                                              (0.6574,          # Capacity fraction at T_supply = 55 °C, Dry_bulb = -7 °C and partial load = 1
                                               0.7141,          # Capacity fraction at T_supply = 55 °C, Dry_bulb = -2 °C and partial load = 1
                                               0.7674,          # Capacity fraction at T_supply = 55 °C, Dry_bulb = 2 °C and partial load = 1
                                               0.9468,          # Capacity fraction at T_supply = 55 °C, Dry_bulb = 7 °C and partial load = 1
                                               0.9942)])        # Capacity fraction at T_supply = 55 °C, Dry_bulb = 12 °C and partial load = 1
                                      , 
                      "COP_fraction": np.array([(0.5943,       # COP fraction at T_supply = 35 °C, Dry_bulb = -7 °C and partial load = 1
                                                 0.7494,          # COP fraction at T_supply = 35 °C, Dry_bulb = -2 °C and partial load = 1
                                                 0.8138,          # COP fraction at T_supply = 35 °C, Dry_bulb = 2 °C and partial load = 1
                                                 1.0000,          # COP fraction at T_supply = 35 °C, Dry_bulb = 7 °C and partial load = 1
                                                 1.0358),         # COP fraction at T_supply = 35 °C, Dry_bulb = 12 °C and partial load = 1
                                                (0.4606,         # COP fraction at T_supply = 45 °C, Dry_bulb = -7 °C and partial load = 1
                                                 0.5107,          # COP fraction at T_supply = 45 °C, Dry_bulb = -2 °C and partial load = 1
                                                 0.5537,          # COP fraction at T_supply = 45 °C, Dry_bulb = 2 °C and partial load = 1
                                                 0.6683,          # COP fraction at T_supply = 45 °C, Dry_bulb = 7 °C and partial load = 1
                                                 0.7064),         # COP fraction at T_supply = 45 °C, Dry_bulb = 12 °C and partial load = 1
                                                (0.3652,         # COP fraction at T_supply = 55 °C, Dry_bulb = -7 °C and partial load = 1
                                                 0.4033,          # COP fraction at T_supply = 55 °C, Dry_bulb = -2 °C and partial load = 1
                                                 0.4368,          # COP fraction at T_supply = 55 °C, Dry_bulb = 2 °C and partial load = 1
                                                 0.5322,          # COP fraction at T_supply = 55 °C, Dry_bulb = 7 °C and partial load = 1
                                                 0.5585)])}       # COP fraction at T_supply = 55 °C, Dry_bulb = 12 °C and partial load = 1
        elif rated_cap == 14.00 and rated_COP == 4.08:
            HP_perf = {"supply_temperature": [35,45,55], 
                      "dry_bulb_temperature": [-7,-2,2,7,12],
                      "Q_fraction": np.array([(0.6429, # Capacity fraction at T_supply = 35 °C, Dry_bulb = -7 °C and partial load = 1
                                               0.6679,           # Capacity fraction at T_supply = 35 °C, Dry_bulb = -2 °C and partial load = 1
                                               0.7071,           # Capacity fraction at T_supply = 35 °C, Dry_bulb = 2 °C and partial load = 1
                                               1.0000,           # Capacity fraction at T_supply = 35 °C, Dry_bulb = 7 °C and partial load = 1
                                               1.1307),          # Capacity fraction at T_supply = 35 °C, Dry_bulb = 12 °C and partial load = 1
                                              (0.5786,          # Capacity fraction at T_supply = 45 °C, Dry_bulb = -7 °C and partial load = 1
                                               0.6050,           # Capacity fraction at T_supply = 45 °C, Dry_bulb = -2 °C and partial load = 1
                                               0.6743,           # Capacity fraction at T_supply = 45 °C, Dry_bulb = 2 °C and partial load = 1
                                               0.9071,           # Capacity fraction at T_supply = 45 °C, Dry_bulb = 7 °C and partial load = 1
                                               0.9779),          # Capacity fraction at T_supply = 45 °C, Dry_bulb = 12 °C and partial load = 1
                                              (0.5286,          # Capacity fraction at T_supply = 55 °C, Dry_bulb = -7 °C and partial load = 1
                                               0.5564,          # Capacity fraction at T_supply = 55 °C, Dry_bulb = -2 °C and partial load = 1
                                               0.6236,          # Capacity fraction at T_supply = 55 °C, Dry_bulb = 2 °C and partial load = 1
                                               0.8571,          # Capacity fraction at T_supply = 55 °C, Dry_bulb = 7 °C and partial load = 1
                                               0.9214)])        # Capacity fraction at T_supply = 55 °C, Dry_bulb = 12 °C and partial load = 1
                                      , 
                      "COP_fraction": np.array([(0.6348,       # COP fraction at T_supply = 35 °C, Dry_bulb = -7 °C and partial load = 1
                                                 0.7279,          # COP fraction at T_supply = 35 °C, Dry_bulb = -2 °C and partial load = 1
                                                 0.8358,          # COP fraction at T_supply = 35 °C, Dry_bulb = 2 °C and partial load = 1
                                                 1.0000,          # COP fraction at T_supply = 35 °C, Dry_bulb = 7 °C and partial load = 1
                                                 1.1569),         # COP fraction at T_supply = 35 °C, Dry_bulb = 12 °C and partial load = 1
                                                (0.4412,         # COP fraction at T_supply = 45 °C, Dry_bulb = -7 °C and partial load = 1
                                                 0.5074,          # COP fraction at T_supply = 45 °C, Dry_bulb = -2 °C and partial load = 1
                                                 0.5735,          # COP fraction at T_supply = 45 °C, Dry_bulb = 2 °C and partial load = 1
                                                 0.7941,          # COP fraction at T_supply = 45 °C, Dry_bulb = 7 °C and partial load = 1
                                                 0.8603),         # COP fraction at T_supply = 45 °C, Dry_bulb = 12 °C and partial load = 1
                                                (0.3480,         # COP fraction at T_supply = 55 °C, Dry_bulb = -7 °C and partial load = 1
                                                 0.4167,          # COP fraction at T_supply = 55 °C, Dry_bulb = -2 °C and partial load = 1
                                                 0.4632,          # COP fraction at T_supply = 55 °C, Dry_bulb = 2 °C and partial load = 1
                                                 0.6029,          # COP fraction at T_supply = 55 °C, Dry_bulb = 7 °C and partial load = 1
                                                 0.7181)])}       # COP fraction at T_supply = 55 °C, Dry_bulb = 12 °C and partial load = 1
        elif rated_cap == 18.00 and rated_COP == 3.90:
            HP_perf = {"supply_temperature": [35,45,55], 
                      "dry_bulb_temperature": [-7,-2,2,7,12],
                      "Q_fraction": np.array([(0.7389, # Capacity fraction at T_supply = 35 °C, Dry_bulb = -7 °C and partial load = 1
                                               0.7844,           # Capacity fraction at T_supply = 35 °C, Dry_bulb = -2 °C and partial load = 1
                                               0.8489,           # Capacity fraction at T_supply = 35 °C, Dry_bulb = 2 °C and partial load = 1
                                               1.0000,           # Capacity fraction at T_supply = 35 °C, Dry_bulb = 7 °C and partial load = 1
                                               1.2056),          # Capacity fraction at T_supply = 35 °C, Dry_bulb = 12 °C and partial load = 1
                                              (0.6278,          # Capacity fraction at T_supply = 45 °C, Dry_bulb = -7 °C and partial load = 1
                                               0.6700,           # Capacity fraction at T_supply = 45 °C, Dry_bulb = -2 °C and partial load = 1
                                               0.7244,           # Capacity fraction at T_supply = 45 °C, Dry_bulb = 2 °C and partial load = 1
                                               0.9206,           # Capacity fraction at T_supply = 45 °C, Dry_bulb = 7 °C and partial load = 1
                                               1.0244),          # Capacity fraction at T_supply = 45 °C, Dry_bulb = 12 °C and partial load = 1
                                              (0.5833,          # Capacity fraction at T_supply = 55 °C, Dry_bulb = -7 °C and partial load = 1
                                               0.6222,          # Capacity fraction at T_supply = 55 °C, Dry_bulb = -2 °C and partial load = 1
                                               0.6733,          # Capacity fraction at T_supply = 55 °C, Dry_bulb = 2 °C and partial load = 1
                                               0.8550,          # Capacity fraction at T_supply = 55 °C, Dry_bulb = 7 °C and partial load = 1
                                               0.9517)])        # Capacity fraction at T_supply = 55 °C, Dry_bulb = 12 °C and partial load = 1
                                      , 
                      "COP_fraction": np.array([(0.6641,       # COP fraction at T_supply = 35 °C, Dry_bulb = -7 °C and partial load = 1
                                                 0.7821,          # COP fraction at T_supply = 35 °C, Dry_bulb = -2 °C and partial load = 1
                                                 0.8974,          # COP fraction at T_supply = 35 °C, Dry_bulb = 2 °C and partial load = 1
                                                 1.0000,          # COP fraction at T_supply = 35 °C, Dry_bulb = 7 °C and partial load = 1
                                                 1.1231),         # COP fraction at T_supply = 35 °C, Dry_bulb = 12 °C and partial load = 1
                                                (0.5462,         # COP fraction at T_supply = 45 °C, Dry_bulb = -7 °C and partial load = 1
                                                 0.6436,          # COP fraction at T_supply = 45 °C, Dry_bulb = -2 °C and partial load = 1
                                                 0.7128,          # COP fraction at T_supply = 45 °C, Dry_bulb = 2 °C and partial load = 1
                                                 0.7821,          # COP fraction at T_supply = 45 °C, Dry_bulb = 7 °C and partial load = 1
                                                 0.9231),         # COP fraction at T_supply = 45 °C, Dry_bulb = 12 °C and partial load = 1
                                                (0.5026,         # COP fraction at T_supply = 55 °C, Dry_bulb = -7 °C and partial load = 1
                                                 0.5077,          # COP fraction at T_supply = 55 °C, Dry_bulb = -2 °C and partial load = 1
                                                 0.5333,          # COP fraction at T_supply = 55 °C, Dry_bulb = 2 °C and partial load = 1
                                                 0.6179,          # COP fraction at T_supply = 55 °C, Dry_bulb = 7 °C and partial load = 1
                                                 0.7231)])}       # COP fraction at T_supply = 55 °C, Dry_bulb = 12 °C and partial load = 1
        else:
            print("\nHeat pump of " + str(round(rated_cap,2)) + "kW rated capacity and" + str(round(rated_COP,2)) + " rated COP is not available in the library."
                  + "\nCheck available sizes or follow the instructions in the documentation to consider a user provided heat pump.")
    elif data_type == "user_provided_HP":
        # User provided heat pump data in Excel format
        
        """Reading the excel file. """
        currentPath     = os.path.dirname("__file__")                
        userPath     = os.path.join(currentPath, 'userProfiles')
        HPper_FileName = "HP_performance_heating.xlsx"
        HPper_File     = os.path.join(userPath, HPper_FileName)
        
        HPper_data = dict(pd.read_excel(HPper_File)) # excel data to dictionary
        
        """Reading data. """
        supply       = []
        Test         = []
        Qth_fraction = []
        Cop_fraction = []
        
        for value in HPper_data["supply_temperature"]:
            if value > 0:
                supply.append(value)
                
        for value in HPper_data["dry_bulb_temperature"]:
            Test.append(value)
            
        for i in range(len(supply)):
            for value in HPper_data["Q_fraction_at_{0}".format(supply[i])]:
                Qth_fraction.append(value)
            for value in HPper_data["COP_fraction_at_{0}".format(supply[i])]:
                Cop_fraction.append(value)
        
        """Collecting data to be given as input in interpolation functions. """
        HP_perf      = {}
        HP_perf["supply_temperature"] = supply
        HP_perf["dry_bulb_temperature"] = Test
        HP_perf["Q_fraction"] = np.array([Qth_fraction[i*(len(Test)):(i+1)*(len(Test))] for i in range(len(supply))])
        HP_perf["COP_fraction"] = np.array([Cop_fraction[i*(len(Test)):(i+1)*(len(Test))] for i in range(len(supply))])
        
    """Creating the interpolation functions."""
    f_Qmax = interp2d(x=HP_perf["dry_bulb_temperature"], y=HP_perf["supply_temperature"], z=HP_perf["Q_fraction"])
    f_COP = interp2d(x=HP_perf["dry_bulb_temperature"], y=HP_perf["supply_temperature"], z=HP_perf["COP_fraction"])
    
    """Calculating Qmax and COP at the desired outdoor and supply temperatures. """
    Qmax   = float(f_Qmax(Toutdoor,Tsupply)*rated_cap*1000)
    COP   = float(f_COP(Toutdoor,Tsupply)*rated_COP)
    
    return (COP,   # 0 - Cefficient Of Performance
            Qmax)  # 1 - Total capacity (W)


