# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 13:42:55 2020

@author: DIISM
"""
import numpy as np
import pandas as pd
from datetime import datetime
start_calc_sim = datetime.now()
"""
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
File in which thermostat set-point profiles can be calculated according to a normal distribution.
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
month    = 1
day      = 1
n_days   = 1
mins     = 30       # in minutes
timestep = mins/60  # in hours

duration_dr_h = 1     # in hours
duration = n_days*24  # in hours

# thermostat:
mu, sigma     = 26.0, 0.15 # mean and standard deviation of temperature for normal (Gaussian) distribution
Tprofile_type = 'hourly'   # temperature profile type: hourly or daily

samples       = int(duration/timestep)
time          = [i*timestep for i in range(0, samples)]
timesteps     = int(24 / timestep)  # Number of timesteps per day

tsp   = {}
tsp_SFH   = {}

if Tprofile_type == 'hourly':
    
    th_sp_Bui1        = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui2        = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui3        = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui4        = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui5        = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui6        = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui7        = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui8        = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui9        = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui10       = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui11       = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui12       = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui13       = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui14       = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui15       = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui16       = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui17       = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui18       = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui19       = np.random.normal(mu, sigma, duration) # °C
    th_sp_Bui20       = np.random.normal(mu, sigma, duration) # °C

    th_sp_SFH06    = []
    th_sp_SFH05_1  = []
    th_sp_SFH05_2  = []
    th_sp_SFH90_1  = []
    th_sp_SFH90_2  = []
    th_sp_SFH90_3  = []
    th_sp_SFH90_4  = []
    th_sp_SFH90_5  = []
    th_sp_SFH90_6  = []
    th_sp_SFH60_1  = []
    th_sp_SFH60_2  = []
    th_sp_SFH60_3  = []
    th_sp_SFH60_4  = []
    th_sp_SFH60_5  = []
    th_sp_SFH60_6  = []
    th_sp_SFH60_7  = []
    th_sp_SFH60_8  = []
    th_sp_SFH60_9  = []
    th_sp_SFH60_10 = []
    th_sp_SFH60_11 = []
    
    for n in range (duration):
        th_sp_SFH06.extend((round((th_sp_Bui1[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH05_1.extend((round((th_sp_Bui2[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH05_2.extend((round((th_sp_Bui3[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH90_1.extend((round((th_sp_Bui4[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH90_2.extend((round((th_sp_Bui5[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH90_3.extend((round((th_sp_Bui6[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH90_4.extend((round((th_sp_Bui7[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH90_5.extend((round((th_sp_Bui8[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH90_6.extend((round((th_sp_Bui9[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH60_1.extend((round((th_sp_Bui10[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH60_2.extend((round((th_sp_Bui11[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH60_3.extend((round((th_sp_Bui12[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH60_4.extend((round((th_sp_Bui13[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH60_5.extend((round((th_sp_Bui14[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH60_6.extend((round((th_sp_Bui15[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH60_7.extend((round((th_sp_Bui16[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH60_8.extend((round((th_sp_Bui17[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH60_9.extend((round((th_sp_Bui18[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH60_10.extend((round((th_sp_Bui19[n]),2)) for i in range (int(60/mins)))
        th_sp_SFH60_11.extend((round((th_sp_Bui20[n]),2)) for i in range (int(60/mins)))


if Tprofile_type == 'daily':
    th_sp_Bui1        = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui2        = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui3        = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui4        = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui5        = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui6        = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui7        = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui8        = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui9        = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui10       = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui11       = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui12       = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui13       = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui14       = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui15       = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui16       = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui17       = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui18       = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui19       = np.random.normal(mu, sigma, n_days) # °C
    th_sp_Bui20       = np.random.normal(mu, sigma, n_days) # °C

    th_sp_SFH06    = []
    th_sp_SFH05_1  = []
    th_sp_SFH05_2  = []
    th_sp_SFH90_1  = []
    th_sp_SFH90_2  = []
    th_sp_SFH90_3  = []
    th_sp_SFH90_4  = []
    th_sp_SFH90_5  = []
    th_sp_SFH90_6  = []
    th_sp_SFH60_1  = []
    th_sp_SFH60_2  = []
    th_sp_SFH60_3  = []
    th_sp_SFH60_4  = []
    th_sp_SFH60_5  = []
    th_sp_SFH60_6  = []
    th_sp_SFH60_7  = []
    th_sp_SFH60_8  = []
    th_sp_SFH60_9  = []
    th_sp_SFH60_10 = []
    th_sp_SFH60_11 = []
    
    for n in range (n_days):
        th_sp_SFH06.extend((round((th_sp_Bui1[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH05_1.extend((round((th_sp_Bui2[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH05_2.extend((round((th_sp_Bui3[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH90_1.extend((round((th_sp_Bui4[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH90_2.extend((round((th_sp_Bui5[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH90_3.extend((round((th_sp_Bui6[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH90_4.extend((round((th_sp_Bui7[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH90_5.extend((round((th_sp_Bui8[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH90_6.extend((round((th_sp_Bui9[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH60_1.extend((round((th_sp_Bui10[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH60_2.extend((round((th_sp_Bui11[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH60_3.extend((round((th_sp_Bui12[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH60_4.extend((round((th_sp_Bui13[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH60_5.extend((round((th_sp_Bui14[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH60_6.extend((round((th_sp_Bui15[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH60_7.extend((round((th_sp_Bui16[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH60_8.extend((round((th_sp_Bui17[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH60_9.extend((round((th_sp_Bui18[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH60_10.extend((round((th_sp_Bui19[n]),2)) for i in range (int(samples/n_days)))
        th_sp_SFH60_11.extend((round((th_sp_Bui20[n]),2)) for i in range (int(samples/n_days)))


time = (time[0:samples])
Tsp_Bui1 = (th_sp_SFH06[0:samples])
Tsp_Bui2 = (th_sp_SFH05_1[0:samples])
Tsp_Bui3 = (th_sp_SFH05_2[0:samples])
Tsp_Bui4 = (th_sp_SFH90_1[0:samples])
Tsp_Bui5 = (th_sp_SFH90_2[0:samples])
Tsp_Bui6 = (th_sp_SFH90_3[0:samples])
Tsp_Bui7 = (th_sp_SFH90_4[0:samples])
Tsp_Bui8 = (th_sp_SFH90_5[0:samples])
Tsp_Bui9 = (th_sp_SFH90_6[0:samples])
Tsp_Bui10 = (th_sp_SFH60_1[0:samples])
Tsp_Bui11 = (th_sp_SFH60_2[0:samples])
Tsp_Bui12 = (th_sp_SFH60_3[0:samples])
Tsp_Bui13 = (th_sp_SFH60_4[0:samples])
Tsp_Bui14 = (th_sp_SFH60_5[0:samples])
Tsp_Bui15 = (th_sp_SFH60_6[0:samples])
Tsp_Bui16 = (th_sp_SFH60_7[0:samples])
Tsp_Bui17 = (th_sp_SFH60_8[0:samples])
Tsp_Bui18 = (th_sp_SFH60_9[0:samples])
Tsp_Bui19 = (th_sp_SFH60_10[0:samples])
Tsp_Bui20 = (th_sp_SFH60_11[0:samples])


dfTsp   = pd.DataFrame(list(zip(time,Tsp_Bui1,Tsp_Bui2,Tsp_Bui3,Tsp_Bui4,Tsp_Bui5,Tsp_Bui6,Tsp_Bui7,
                                Tsp_Bui8,Tsp_Bui9,Tsp_Bui10,Tsp_Bui11,Tsp_Bui12,Tsp_Bui13,Tsp_Bui14,Tsp_Bui15,
                                Tsp_Bui16,Tsp_Bui17,Tsp_Bui18,Tsp_Bui19,Tsp_Bui20)), 
                    columns= ['time (hr)','Tsp_Bui1','Tsp_Bui2','Tsp_Bui3','Tsp_Bui4','Tsp_Bui5','Tsp_Bui6','Tsp_Bui7',
                                                    'Tsp_Bui8','Tsp_Bui9','Tsp_Bui10','Tsp_Bui11','Tsp_Bui12','Tsp_Bui13','Tsp_Bui14','Tsp_Bui15',
                                                    'Tsp_Bui16','Tsp_Bui17','Tsp_Bui18','Tsp_Bui19','Tsp_Bui20'], dtype = float)


dfTsp.to_csv('./T_sp.csv', index = False)

print(tsp)