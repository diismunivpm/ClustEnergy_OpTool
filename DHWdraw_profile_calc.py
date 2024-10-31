# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 13:42:55 2020

@author: DIISM
"""
import numpy as np
import pandas as pd
from datetime import datetime
from User_pattern import occProfile

start_calc_sim = datetime.now()
"""
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
File in which DHW profiles can be calculated according to according to occupancy and maximum daily water consumption (e.g., 200 l/day).

During active occupation, a random amount of water consumption is assigned based on probabilities derived from the profiles 
specified in EN 15316-3-1. These profiles link each activity type (such as handwashing, showering, cleaning, etc.) with a 
corresponding amount of water, estimated by the time needed to complete the activity.
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
month    = 1
day      = 1
n_days   = 365
mins     = 15       # in minutes
timestep = mins/60  # in hours

n_Bui = 20            # number of buildings
duration = n_days*24  # in hours

samp_day = int(24/timestep)
samples       = int(duration/timestep)
time          = [i*timestep for i in range(0, samples)]
timesteps     = int(24 / timestep)  # Number of timesteps per day

"""Definition of water quatities to carry out certain household activities and their probabilities."""
water_amounts    = [6.02, 40.13, 3.01, 6.02, 3.01]     # for small utilization, shower, floor cleaning, dish washing and household cleaning
probabilities    = [0.696, 0.087, 0.043, 0.087, 0.087] # defining probabilities for each water draw (the total must be 1)
max_daily_water  = 206                                 # liters - maximum allowed daily consumption of DHW

occupancy   = {}
"""Reading the occupancy patterns."""
for i in range(n_Bui):
    occupancy["occ_{0}".format(i+1)] = occProfile(n_Bui,day,month,duration,timestep)[i]

draw = {}

# Generate water values for each building based on occupancy and maximum daily consumption
for i in range(n_Bui):
    water_for_building = []  # List to store water values
    for d in range(n_days):
        # Calculate the maximum allowed daily consumption
        max_daily_water = max_daily_water
        
        # Initialize the total daily consumption to zero
        daily_consumption = 0
        
        for t in range(timesteps):
            if occupancy["occ_{0}".format(i+1)][d*timesteps + t] > 0:
                # Calculate the maximum water consumption for this time interval
                max_interval_water = max_daily_water - daily_consumption
                
                if max_interval_water > 0:
                    # Randomly generate water consumption for this time interval
                    water_draw = np.random.choice(water_amounts, p=probabilities, size=1)
                    
                    # Limit water consumption to the maximum allowed interval
                    water_draw = min(water_draw, max_interval_water)
                    
                    # Update the total daily consumption
                    daily_consumption += water_draw
                
                # Add the drawn value to the list
                water_for_building.append(float(water_draw))
            else:
                water_for_building.append(0)  # Add 0 if occupancy is zero
                
    draw["draw_building_{0}".format(i+1)] = [water_for_building[n] for n in range(samples)]  # Add the list of water values to the corresponding building



draw = pd.DataFrame.from_dict(draw)

# export
draw.to_csv('./DHWdraw_profile' + str(mins) +'mins_step.csv', index = False)
