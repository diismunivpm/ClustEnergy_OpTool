"""
ClustEnergy OpTool
Released on January, 2024
@author: DIISM UNIVPM
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
import os
import sys
import numpy as np
import pandas as pd
import math

"""Solar gain calculation:
        (1) position of the sun
        (2) opaque components gains
        (3) transparent components gains. """

class Weather(object): 
   def __init__(self, epwfile_path):
        # Set EPW Labels and import epw file
        epw_labels = ['year', 'month', 'day', 'hour', 'minute', 'datasource', 'drybulb_C', 'dewpoint_C', 'relhum_percent',
                      'atmos_Pa', 'exthorrad_Whm2', 'extdirrad_Whm2', 'horirsky_Whm2', 'glohorrad_Whm2',
                      'dirnorrad_Whm2', 'difhorrad_Whm2', 'glohorillum_lux', 'dirnorillum_lux', 'difhorillum_lux',
                      'zenlum_lux', 'winddir_deg', 'windspd_ms', 'totskycvr_tenths', 'opaqskycvr_tenths', 'visibility_km',
                      'ceiling_hgt_m', 'presweathobs', 'presweathcodes', 'precip_wtr_mm', 'aerosol_opt_thousandths',
                      'snowdepth_cm', 'days_last_snow', 'Albedo', 'liq_precip_depth_mm', 'liq_precip_rate_Hour']
         
        # Import EPW file
        self.weather_data = pd.read_csv(epwfile_path, skiprows=8, header=None, names=epw_labels).drop('datasource', axis=1)        
   
   """Extrapolate weather data from epw file."""
   @property
   def hour(self):        
        return self.weather_data['hour'].tolist()            # time of the day (hr)
   
   @property
   def drybulb_C(self):                                      # Dry bulb temperature (°C)
        return self.weather_data['drybulb_C'].tolist()
    
   @property
   def dewpoint_C(self):                                     # Dew point temperature (°C)
        return self.weather_data['dewpoint_C'].tolist()
    
   @property
   def relhum_percent(self):                                 # Relative humidity (%)
        return self.weather_data['relhum_percent'].tolist()
    
   @property
   def atmos_Pa(self):                                       # Atmospheric pressure (Pa)
        return self.weather_data['atmos_Pa'].tolist()
    
   @property
   def dirnorrad_Whm2(self):                                 # Direct normal radiation (Wh/m2)
        return self.weather_data['dirnorrad_Whm2'].tolist()
    
   @property
   def difhorrad_Whm2(self):                                 # Diffuse horizontal radiation(Wh/m2)
        return self.weather_data['difhorrad_Whm2'].tolist()
    
   @property
   def glohorrad_Whm2(self):                                 # Global horizontal radiation (Wh/m2)
        return self.weather_data['glohorrad_Whm2'].tolist()

   @property
   def extdirrad_Whm2(self):                                 # Extraterrestrial horizontal radiation (Wh/m2)
        return self.weather_data['extdirrad_Whm2'].tolist()
        
   def calc_sun_position(self, latitude_deg, longitude_deg, dirnorrad_Whm2, difhorrad_Whm2, glohorrad_Whm2, time, UTC, timestep):
        """The calc_sun_position calculates the position of the sun during the simulation period. """
        
        latitude = latitude_deg
        longitude = longitude_deg
        latitude_rad = math.radians(latitude) # latitude in radians
        
        # time conversion
        day_of_year = int (time/(24))  
        time_day    = time -timestep/2-(day_of_year)*24
        time_day_th    = time -(day_of_year)*24 #for thermostat
        
        # Fractional year:
        B = 2*math.pi*(day_of_year-1)/365
        
        # Equation of time (in minutes) - eq. of Spencer (1971)
        E = (24/(2*math.pi))*(0.000075-0.001868*math.cos(B)-0.032077*math.sin(B)-0.014615*math.cos(2*B)-0.04089*math.sin(2*B))

        # Solar time is the time used in all of the sun-angle relationships
        Solar_time = time_day + E + (UTC*15 - longitude)/15
        
        # Declination
        Declination = 23.45*math.sin(math.radians((360/365)*(284+day_of_year)))  # Solar declination (ASHRAE,2007): equation of Cooper (1969)
        Declination_rad = math.radians(Declination)
        
        # Solar hour angle (degrees)
        w = 15 * (Solar_time - 12)
        w_rad = math.radians(w)        
        
        """Calculating the sun zenith. """
        cos_zenith  = math.sin(Declination_rad)*math.sin(latitude_rad)+math.cos(latitude_rad)*math.cos(Declination_rad)*math.cos(w_rad)           
        sin_zenith  = math.sin(math.acos(cos_zenith))
        
        if cos_zenith>0:
            zenith = math.degrees(math.acos(cos_zenith))
        else:
            zenith = 90 
        
        """Calculating the sun azimuth. """
        sin_azimuth = math.cos(Declination_rad)*math.sin(w_rad)/sin_zenith
        cos_azimuth = (math.cos(Declination_rad)*math.cos(w_rad)*math.sin(latitude_rad)-math.sin(Declination_rad)*math.cos(latitude_rad))/sin_zenith
        
        if zenith < 90:
            if sin_azimuth < 0:
                azimuth = -math.degrees(math.acos(cos_azimuth))
            else:
                azimuth = math.degrees(math.acos(cos_azimuth))
        else:
            azimuth = -90
        
        """Calculating the solar altitude angle (function only of time of day and declination). """
        altitude = math.degrees(math.asin(math.cos(latitude_rad) * math.cos(Declination_rad)* math.cos(w_rad) + math.sin(latitude_rad)* math.sin(Declination_rad)))

        return (altitude,     # 0 - sun altitude
                azimuth,      # 1 - sun azimuth
                zenith,       # 2 - sun zenith
                time_day_th)  # 3 - time of the day
     
class Window_solar(object):
         
   def __init__(self, azimut_surface, inclination_surface, glass_solar_transmittance, area):
       
        self.azimut_surface = azimut_surface 
        self.inclination_surface = inclination_surface   
        self.inclination_rad = math.radians(inclination_surface)
        self.azimuth_surf_rad = math.radians(azimut_surface)
        self.glass_solar_transmittance = glass_solar_transmittance              # (g_gl_n) - glass solar transmittance
        self.area = area                                                        # (A_w) - window area

   def window_solar_gains(self, sun_altitude, sun_azimut, sun_zenit, dirnorrad_Whm2, difhorrad_Whm2, time, latitude_deg, longitude_deg, UTC, glohorrad_Whm2, ground_reflectance, timestep):
        """
        The window_solar_gains function calculates solar gains through transparent surfaces (i.e., windows). 
        """
        solar_zenith_rad  = math.radians(sun_zenit)   # solar zenith in radians
        solar_azimuth_rad = math.radians(sun_azimut)  # solar azimuth in radians
        
        # cosine of incidence angle
        cos_teta = (math.cos(solar_zenith_rad)*math.cos(self.inclination_rad)+math.sin(solar_zenith_rad)*math.cos(solar_azimuth_rad-self.azimuth_surf_rad)*math.sin(self.inclination_rad))

        # Rb = beam radiation factor
        if cos_teta > 0:
            Rb = cos_teta
        else:      
            Rb = 0
        
        # teta = incidence angle
        teta = math.acos(cos_teta)
        
        # Rd = diffuse radiation factor
        Rd = (1 + math.cos(math.radians(self.inclination_surface)))/2        

        # Rr = reflected radiation
        Rr = ((1 - math.cos(math.radians(self.inclination_surface)))/2)*ground_reflectance
        
        self.teta = teta
        self.Rb = Rb
        self.direct_radiation_Whm2 =  dirnorrad_Whm2 * Rb                                                                                             # Direct radiation (Wh/m2)
        self.ground_reflected_Whm2 = glohorrad_Whm2 * Rr                                                                                              # Ground reflected radiation (Wh/m2)
        self.sky_diffuse_Whm2 = difhorrad_Whm2 * Rd                                                                                                   # Sky diffuse radiation (Wh/m2)
        self.total_diffuse_Whm2 = self.ground_reflected_Whm2 + self.sky_diffuse_Whm2                                                                  # Total diffuse radiation (Wh/m2)
        self.solar_radiation_Whm2 = self.direct_radiation_Whm2 + self.total_diffuse_Whm2                                                              # Total solar radiation (Wh/m2)
        self.incident_solar = (self.solar_radiation_Whm2) * self.area                                                                                 # Total incident solar radiation (Wh)
        self.windows_solar_gains = self.glass_solar_transmittance*self.area*(0.84*self.total_diffuse_Whm2 + self.direct_radiation_Whm2*math.sqrt(Rb)) # Calculating windows solar gains (Wh)
        
class Wall_solar(object):
    
    def __init__(self, azimut_surface, inclination_surface, absorption, surface_thermal_resistance, equivalent_transmittance, area):
       
        self.azimut_surface = azimut_surface 
        self.inclination_surface = inclination_surface
        self.inclination_rad = math.radians(inclination_surface)
        self.azimuth_surf_rad = math.radians(azimut_surface)
        self.absorption = absorption                                            # (alfa) - wall absorption
        self.surface_thermal_resistance = surface_thermal_resistance            # (R_se) - surface thermal resistance
        self.equivalent_transmittance = equivalent_transmittance                # (U_c_eq) - equivalent transmittance
        self.area = area                                                        # (A_c) - surface area
        
    def wall_solar_gains(self, sun_altitude, sun_azimut, sun_zenit, dirnorrad_Whm2, difhorrad_Whm2, time, latitude_deg, longitude_deg, UTC , glohorrad_Whm2, ground_reflectance, timestep):
        """
        The wall_solar_gains function calculates solar gains through opaque surfaces (i.e., walls and roof). 
        """
        
        solar_zenith_rad  = math.radians(sun_zenit)   # solar zenith in radians
        solar_azimuth_rad = math.radians(sun_azimut)  # solar azimuth in radians
        
        # cosine of incidence angle:
        cos_teta = (math.cos(solar_zenith_rad)*math.cos(self.inclination_rad)+math.sin(solar_zenith_rad)*math.cos(solar_azimuth_rad-self.azimuth_surf_rad)*math.sin(self.inclination_rad))

        # Rb = beam radiation factor
        if cos_teta > 0:
            Rb = cos_teta
        else:      
            Rb = 0
        
        # teta = incidence angle
        teta = ((math.acos(cos_teta)))
        
        # Rd = diffuse radiation factor
        Rd = (1 + math.cos(math.radians(self.inclination_surface)))/2        

        # Rr = reflected radiation
        Rr = ((1 - math.cos(math.radians(self.inclination_surface)))/2)*ground_reflectance
        
        self.teta = teta
        self.Rb = Rb
        self.direct_radiation_Whm2 =  dirnorrad_Whm2 * Rb                                                                                # Direct radiation (Wh/m2)
        self.ground_reflected_Whm2 = glohorrad_Whm2 * Rr                                                                                 # Ground reflected radiation (Wh/m2)
        self.sky_diffuse_Whm2 = difhorrad_Whm2 * Rd                                                                                      # Sky diffuse radiation (Wh/m2)
        self.total_diffuse_Whm2 = self.ground_reflected_Whm2 + self.sky_diffuse_Whm2                                                     # Total diffuse radiation (Wh/m2)
        self.solar_radiation_Whm2 = self.direct_radiation_Whm2 + self.total_diffuse_Whm2                                                 # Total solar radiation (Wh/m2)
        self.incident_solar = (self.solar_radiation_Whm2) * self.area                                                                    # Total incident solar radiation (Wh)
        self.walls_solar_gains = self.incident_solar * self.absorption * self.surface_thermal_resistance * self.equivalent_transmittance # Calculating walls solar gains (Wh)

def CheckMeteoFile(loc): 
    """The CheckMeteoFile function checks whether the chosen location is available.
            In case a new weather file is added, we recommend adding the location in the function. """
    if loc != 'Milano' and loc != 'Roma' and loc !='Ancona' and loc != 'Messina' and loc != 'Denver': 
       sys.exit('Weather file not available')

def MeteoFile(loc,day,month,duration,timestep): 
    """The MeteoF function opens the epw file of the chosen location. 
            It is possible to add further locations by following the steps: 
                •	insert the climate file in epw format in the "epwFile" folder, 
                •	add the climate file in the check function "CheckMeteoFile" in "MeteoFile.py" with the following code line in the if statement
                    and loc != 'NEW_LOCATION' 
                •	to open the file, in the "MeteoF" function in "MeteoFile.py" insert the following lines of code
                    if loc == 'Milano': 
                        psychrolib.SetUnitSystem(psychrolib.SI)  
                        currentPath     = os.path.dirname("__file__")    
                        weatherPath     = os.path.join(currentPath, 'epwFile')
                        weatherFileName = "NEW_LOCATION.epw"
                        weatherFile     = os.path.join(weatherPath, weatherFileName)
                        epw             = Weather(weatherFile) 
                        latitude_deg    = LAT(NEW_LOCATION)
                        longitude_deg   = LONG(NEW_LOCATION)
                        UTC             = UTC(NEW_LOCATION)"""
    if loc == 'Milano': 
        currentPath     = os.path.dirname("__file__")                
        weatherPath     = os.path.join(currentPath, 'epwFile')
        weatherFileName = "ITA_Milano-Linate.epw"
        weatherFile     = os.path.join(weatherPath, weatherFileName)
        epw             = Weather(weatherFile) 
        latitude_deg    = 45.43
        longitude_deg   = 9.28
        UTC             = 1.0
    if loc == 'Roma': 
        currentPath     = os.path.dirname("__file__")                
        weatherPath     = os.path.join(currentPath, 'epwFile')
        weatherFileName = "ITA_Roma-Ciampino.epw"
        weatherFile     = os.path.join(weatherPath, weatherFileName)
        epw             = Weather(weatherFile)   
        latitude_deg    = 41.80
        longitude_deg   = 18.58
        UTC             = 1.0
    if loc == 'Ancona': 
        currentPath     = os.path.dirname("__file__")                
        weatherPath     = os.path.join(currentPath, 'epwFile')
        weatherFileName = "ITA_Ancona.epw"
        weatherFile     = os.path.join(weatherPath, weatherFileName)
        epw             = Weather(weatherFile)      
        latitude_deg    = 43.62
        longitude_deg   = 13.52
        UTC             = 1.0
    if loc == 'Messina': 
        currentPath     = os.path.dirname("__file__")                
        weatherPath     = os.path.join(currentPath, 'epwFile')
        weatherFileName = "ITA_Messina.epw"
        weatherFile     = os.path.join(weatherPath, weatherFileName)
        epw             = Weather(weatherFile)            
        latitude_deg    = 38.20
        longitude_deg   = 15.55
        UTC             = 1.0
    if loc == 'Denver': # Denver International Airport
        currentPath     = os.path.dirname("__file__")                
        weatherPath     = os.path.join(currentPath, 'epwFile')
        weatherFileName = "USA_Denver.epw"
        weatherFile     = os.path.join(weatherPath, weatherFileName)
        epw             = Weather(weatherFile) 
        latitude_deg    = 39.83
        longitude_deg   = -104.65
        UTC             = -7.0
    
    """Simulation period calculation. """
    """Calculation of the start time of the chosen month. """
    if month == 1:   # January
        day_s = 0
    if month == 2:   # February
        day_s = 744
    if month == 3:   # March
        day_s = 1416
    if month == 4:   # April
        day_s = 2160
    if month == 5:   # May
        day_s = 2880
    if month == 6:   # June
        day_s = 3624
    if month == 7:   # July
        day_s = 4344
    if month == 8:   # August
        day_s = 5088
    if month == 9:   # September
        day_s = 5832
    if month == 10:  # October
        day_s = 6552
    if month == 11:  # November
        day_s = 7296
    if month == 12:  # December
        day_s = 8016
    
    start        = int((day_s + (day-1)*24))             # start time of simulation
    samples      = int(duration/timestep)                # number of steps
    stop         = int(start + duration)                 # stop time of simulation
    samples_h    = int((stop-start))                     # duration in hours
    duration_sim = range(0,samples)                      # total duration of simulation
    hours        = np.linspace(start, stop, samples_h)   # time array in hours
    time         = np.linspace(start, stop, samples)     # time array
    
    """Extrapolation of data from epw file. """
    time_day        = []
    T_outdoor       = []
    Difhorrad_Whm2  = []
    Dirnorrad_Whm2  = []
    Glohorrad_Whm2  = []
    
    for i in range(start,stop): 
        time_day.append(epw.hour[i-1])
        T_outdoor.append(epw.drybulb_C[i-1])           # outdoor temperature
        Difhorrad_Whm2.append(epw.difhorrad_Whm2[i-1]) # diffuse solar radiation
        Dirnorrad_Whm2.append(epw.dirnorrad_Whm2[i-1]) # direct solar radiation
        Glohorrad_Whm2.append(epw.glohorrad_Whm2[i-1]) # global solar radiation
    
    """Data interpolation. """
    hour      = np.interp(time, hours, time_day).tolist()
    T_est     = np.interp(time, hours, T_outdoor).tolist()
    Diff_Whm2 = np.interp(time, hours, Difhorrad_Whm2).tolist()
    Dir_Whm2  = np.interp(time, hours, Dirnorrad_Whm2).tolist()
    Glo_Whm2  = np.interp(time, hours, Glohorrad_Whm2).tolist()
    
    """Calculating the sun altitude, azimuth and zenith via the calc_sun_position function. """
    sun_altitude    = []
    sun_azimut      = []
    sun_zenit       = []
    
    for i in duration_sim:
        time_d = (start+(i*timestep)) 
        sun_altitude.append(epw.calc_sun_position(latitude_deg, longitude_deg ,Dir_Whm2[i],Diff_Whm2[i], Glo_Whm2[i], time_d, UTC, timestep)[0])
        sun_azimut.append(epw.calc_sun_position(latitude_deg, longitude_deg, Dir_Whm2[i],Diff_Whm2[i], Glo_Whm2[i], time_d, UTC, timestep)[1])
        sun_zenit.append(epw.calc_sun_position(latitude_deg, longitude_deg, Dir_Whm2[i],Diff_Whm2[i], Glo_Whm2[i], time_d, UTC, timestep)[2])
    
    """Ground temperature. """
    T_ground        = []
    for i in range(0,len(T_est)):
        T_ground.append(T_est[i]) # Ground temperature equal to outdoor air temperature, according to ANSI/ASHRAE Standard 140 (i.e., BESTEST)
    
    """Outputs. """
    return(T_est,           # 0 - outdoor temperature trend
           T_ground,        # 1 - ground temperature
           Diff_Whm2,       # 2 - diffuse solar radiation
           Dir_Whm2,        # 3 - direct solar radiation
           Glo_Whm2,        # 4 - global solar radiation
           sun_altitude,    # 5 - altitude of the sun
           sun_azimut,      # 6 - azimuth of the sun
           sun_zenit,       # 7 - zenith of the sun
           latitude_deg,    # 8 - latitude of the location
           longitude_deg,   # 9 - longitude of the location 
           start,           # 10 - simulation start time
           UTC,             # 11 - Coordinated Universal Time of the location
           hour)            # 12 - time of the day

"""Calculating solar gain factors. """
def SF_double_glazing(roof_area, wall_wide, wall_long, wall_high, S_wind_area, E_wind_area, W_wind_area, N_wind_area, wall_absorption, U_window, 
                      ext_heat_trans_coef, int_heat_trans_coef, window_reflectance, ext_window_absorption, int_window_absorption): 

    windows_area = S_wind_area + E_wind_area + W_wind_area + N_wind_area
    
    south_area    = wall_wide*wall_high - S_wind_area
    north_area    = wall_wide*wall_high - N_wind_area
    east_area     = wall_long*wall_high - E_wind_area
    west_area     = wall_long*wall_high - W_wind_area
    
    floor_area        = wall_long*wall_wide                       # Floor Area (m2)

    total_opaque_area = (wall_wide*wall_high)*2+ (wall_long*wall_high)*2 - (windows_area) # Total opaque Area (m2)

    total_surf_area   = floor_area + roof_area + total_opaque_area
    
    R_window = 1/U_window
    R_int    = 1/int_heat_trans_coef
    R_ext    = 1/ext_heat_trans_coef
    
    N_int    = R_int/R_window # inward conducted fraction of cavity reflected absorbed solar radiation for inner pane
    N_ext    = R_ext/R_window # inward conducted fraction of cavity reflected absorbed solar radiation for outer pane
    
    Fraction_window_S   = S_wind_area/(wall_wide*wall_high)
    Fraction_wall_S = 1 - Fraction_window_S
    
    Fraction_window_N   = N_wind_area/(wall_wide*wall_high)
    Fraction_wall_N = 1 - Fraction_window_N
    
    Fraction_window_E   = E_wind_area/(wall_long*wall_high)
    Fraction_wall_E = 1 - Fraction_window_E
    
    Fraction_window_W   = W_wind_area/(wall_long*wall_high)
    Fraction_wall_W = 1 - Fraction_window_W
    
    """
    ---------------------------------------------------------------------------------------------------------------
    B1 describes the first bounce of incident shortwave radiation, assuming all of it initially hits the floor
    ---------------------------------------------------------------------------------------------------------------
    """
    B1_floor   = wall_absorption
    B1_south   = 0.000
    B1_east    = 0.000
    B1_west    = 0.000
    B1_north   = 0.000
    B1_ceiling = 0.000

    """
    ---------------------------------------------------------------------------------------------------------------
    B2 describes the second bounce such that shortwave radiation diffusely reflected by the floor is distributed 
    over other surfaces in proportion to their view factor-absorptance product
    ---------------------------------------------------------------------------------------------------------------
    """

    y = wall_long
    x = wall_wide
    z = wall_high
    Y = y/x
    Z = z/x

    FF_floor_south = (1/(math.pi*Y))*(Y*math.atan(1/Y) + Z*math.atan(1/Z) - math.sqrt(Z**2 + Y**2)*math.atan(1/(math.sqrt(Z**2 + Y**2))) + 
                                      1/4*math.log(((1+Y**2)*(1+Z**2)/(1+Y**2+Z**2))*(((Y**2)*(1+Y**2+Z**2)/((1+Y**2)*(Y**2+Z**2)))**(Y**2))*
                                      (((Z**2)*(1+Z**2+Y**2)/((1+Z**2)*(Z**2+Y**2)))**(Z**2))))

    FF_floor_north = FF_floor_south

    y = wall_wide
    x = wall_long
    z = wall_high
    Y = y/x
    Z = z/x

    FF_floor_east = (1/(math.pi*Y))*(Y*math.atan(1/Y) + Z*math.atan(1/Z) - math.sqrt(Z**2 + Y**2)*math.atan(1/(math.sqrt(Z**2 + Y**2))) + 
                                      1/4*math.log(((1+Y**2)*(1+Z**2)/(1+Y**2+Z**2))*(((Y**2)*(1+Y**2+Z**2)/((1+Y**2)*(Y**2+Z**2)))**(Y**2))*
                                      (((Z**2)*(1+Z**2+Y**2)/((1+Z**2)*(Z**2+Y**2)))**(Z**2))))

    FF_floor_west = FF_floor_east

    y = wall_wide
    x = wall_long
    D = wall_high
    X = x/D
    Y = y/D

    FF_floor_ceiling = (2/(math.pi*X*Y))*((math.log((1+X**2)*(1+Y**2)/(1+X**2+Y**2)))**(1/2) + X*math.sqrt(1+Y**2)*math.atan(X/math.sqrt(1+Y**2)) +
                                          Y*math.sqrt(1+X**2)*math.atan(Y/math.sqrt(1+X**2)) - X*math.atan(X) - Y*math.atan(Y))

    B2_floor   = 0
    B2_south   = (1 - wall_absorption)*(FF_floor_south*Fraction_wall_S)*wall_absorption
    B2_north   = (1 - wall_absorption)*(FF_floor_north*Fraction_wall_N)*wall_absorption
    B2_east    = (1 - wall_absorption)*(FF_floor_east*Fraction_wall_E)*wall_absorption
    B2_west    = (1 - wall_absorption)*(FF_floor_west*Fraction_wall_W)*wall_absorption
    B2_ceiling = (1 - wall_absorption)*FF_floor_ceiling*wall_absorption
    
    if S_wind_area > 0:
        B2_window_S = (1 - wall_absorption)*(FF_floor_south*Fraction_window_S)*(1-(window_reflectance + N_int*int_window_absorption + N_ext*ext_window_absorption))
    else:
        B2_window_S = 0
    
    if N_wind_area > 0:
        B2_window_N = (1 - wall_absorption)*(FF_floor_north*Fraction_window_N)*(1-(window_reflectance + N_int*int_window_absorption + N_ext*ext_window_absorption))
    else:
        B2_window_N = 0
    
    if  E_wind_area > 0:
        B2_window_E = (1 - wall_absorption)*(FF_floor_east*Fraction_window_E)*(1-(window_reflectance + N_int*int_window_absorption + N_ext*ext_window_absorption))
    else:
        B2_window_E = 0
        
    if W_wind_area > 0:
        B2_window_W = (1 - wall_absorption)*(FF_floor_south*Fraction_window_W)*(1-(window_reflectance + N_int*int_window_absorption + N_ext*ext_window_absorption))
    else:
        B2_window_W = 0

    """
    ---------------------------------------------------------------------------------------------------------------
    B3 describes the third bounce such that the remaining nonabsorbed shortwave radiation is distributed over each
    surf in proportion to its area-absorptance product
    ---------------------------------------------------------------------------------------------------------------
    """
    Sum_B2 = (B2_floor + B2_south + B2_north + B2_east + B2_west + B2_ceiling + 
              B2_window_S +  B2_window_N + B2_window_E + B2_window_W)

    B3_floor   = (1 - wall_absorption - Sum_B2)*(floor_area/total_surf_area)*wall_absorption
    B3_south   = (1 - wall_absorption - Sum_B2)*(south_area/total_surf_area)*wall_absorption
    B3_east    = (1 - wall_absorption - Sum_B2)*(east_area/total_surf_area)*wall_absorption
    B3_west    = (1 - wall_absorption - Sum_B2)*(west_area/total_surf_area)*wall_absorption
    B3_north   = (1 - wall_absorption - Sum_B2)*(north_area/total_surf_area)*wall_absorption
    B3_ceiling = (1 - wall_absorption - Sum_B2)*(roof_area/total_surf_area)*wall_absorption
    
    if S_wind_area > 0:
        B3_window_S  = (1 - wall_absorption - Sum_B2)*(S_wind_area/total_surf_area)*(1-(window_reflectance + N_int*int_window_absorption + N_ext*ext_window_absorption))
    else:
        B3_window_S = 0
    
    if N_wind_area > 0:
        B3_window_N  = (1 - wall_absorption - Sum_B2)*(N_wind_area/total_surf_area)*(1-(window_reflectance + N_int*int_window_absorption + N_ext*ext_window_absorption))
    else:
        B3_window_N = 0
    
    if  E_wind_area > 0:
        B3_window_E  = (1 - wall_absorption - Sum_B2)*(E_wind_area/total_surf_area)*(1-(window_reflectance + N_int*int_window_absorption + N_ext*ext_window_absorption))
    else:
        B3_window_E = 0
        
    if W_wind_area > 0:
        B3_window_W  = (1 - wall_absorption - Sum_B2)*(W_wind_area/total_surf_area)*(1-(window_reflectance + N_int*int_window_absorption + N_ext*ext_window_absorption))
    else:
        B3_window_W = 0

    """
    ---------------------------------------------------------------------------------------------------------------
    BR describes the distribution of all remaining bounces based on distribution fractions from calculations for 
    B3n above
    ---------------------------------------------------------------------------------------------------------------
    """
    Sum_B3 = (B3_floor + B3_south + B3_north + B3_east + B3_west + B3_ceiling + 
              B3_window_S + B3_window_N + B3_window_E + B3_window_W)

    BR_floor   = (1 - wall_absorption - Sum_B2 - Sum_B3)*(B3_floor/Sum_B3)
    BR_south   = (1 - wall_absorption - Sum_B2 - Sum_B3)*(B3_south/Sum_B3)
    BR_east    = (1 - wall_absorption - Sum_B2 - Sum_B3)*(B3_east/Sum_B3)
    BR_west    = (1 - wall_absorption - Sum_B2 - Sum_B3)*(B3_west/Sum_B3)
    BR_north   = (1 - wall_absorption - Sum_B2 - Sum_B3)*(B3_north/Sum_B3)
    BR_ceiling = (1 - wall_absorption - Sum_B2 - Sum_B3)*(B3_ceiling/Sum_B3)


    SF_floor   = round(B1_floor + B2_floor + B3_floor + BR_floor,3)
    SF_south   = round(B1_south + B2_south + B3_south + BR_south,3)
    SF_east    = round(B1_east + B2_east + B3_east + BR_east,3)
    SF_west    = round(B1_west + B2_west + B3_west + BR_west,3)
    SF_north   = round(B1_north + B2_north + B3_north + BR_north,3)
    SF_ceiling = round(B1_ceiling + B2_ceiling + B3_ceiling + BR_ceiling,3)
        
    return(SF_south, SF_east, SF_west, SF_north, SF_ceiling, SF_floor)
    
def SF_single_glazing(roof_area, wall_wide, wall_long, wall_high, S_wind_area, E_wind_area, W_wind_area, N_wind_area, wall_absorption, U_window, 
                      int_heat_trans_coef, window_reflectance, int_window_absorption): 

    windows_area = S_wind_area + E_wind_area + W_wind_area + N_wind_area
    
    south_area    = wall_wide*wall_high - S_wind_area
    north_area    = wall_wide*wall_high - N_wind_area
    east_area     = wall_long*wall_high - E_wind_area
    west_area     = wall_long*wall_high - W_wind_area
    
    floor_area        = wall_long*wall_wide                       # Floor Area (m2)

    total_opaque_area = (wall_wide*wall_high)*2+ (wall_long*wall_high)*2 - (windows_area) # Total opaque Area (m2)

    total_surf_area   = floor_area + roof_area + total_opaque_area
    
    R_window = 1/U_window
    R_int    = 1/int_heat_trans_coef

    
    N_int    = R_int/R_window # inward conducted fraction of cavity reflected absorbed solar radiation for inner pane
    
    Fraction_window_S   = S_wind_area/(wall_wide*wall_high)
    Fraction_wall_S = 1 - Fraction_window_S
    
    Fraction_window_N   = N_wind_area/(wall_wide*wall_high)
    Fraction_wall_N = 1 - Fraction_window_N
    
    Fraction_window_E   = E_wind_area/(wall_long*wall_high)
    Fraction_wall_E = 1 - Fraction_window_E
    
    Fraction_window_W   = W_wind_area/(wall_long*wall_high)
    Fraction_wall_W = 1 - Fraction_window_W
    
    """
    ---------------------------------------------------------------------------------------------------------------
    B1 describes the first bounce of incident shortwave radiation, assuming all of it initially hits the floor
    ---------------------------------------------------------------------------------------------------------------
    """
    B1_floor   = wall_absorption
    B1_south   = 0.000
    B1_east    = 0.000
    B1_west    = 0.000
    B1_north   = 0.000
    B1_ceiling = 0.000

    """
    ---------------------------------------------------------------------------------------------------------------
    B2 describes the second bounce such that shortwave radiation diffusely reflected by the floor is distributed 
    over other surfaces in proportion to their view factor-absorptance product
    ---------------------------------------------------------------------------------------------------------------
    """

    y = wall_long
    x = wall_wide
    z = wall_high
    Y = y/x
    Z = z/x

    FF_floor_south = (1/(math.pi*Y))*(Y*math.atan(1/Y) + Z*math.atan(1/Z) - math.sqrt(Z**2 + Y**2)*math.atan(1/(math.sqrt(Z**2 + Y**2))) + 
                                      1/4*math.log(((1+Y**2)*(1+Z**2)/(1+Y**2+Z**2))*(((Y**2)*(1+Y**2+Z**2)/((1+Y**2)*(Y**2+Z**2)))**(Y**2))*
                                      (((Z**2)*(1+Z**2+Y**2)/((1+Z**2)*(Z**2+Y**2)))**(Z**2))))

    FF_floor_north = FF_floor_south

    y = wall_wide
    x = wall_long
    z = wall_high
    Y = y/x
    Z = z/x

    FF_floor_east = (1/(math.pi*Y))*(Y*math.atan(1/Y) + Z*math.atan(1/Z) - math.sqrt(Z**2 + Y**2)*math.atan(1/(math.sqrt(Z**2 + Y**2))) + 
                                      1/4*math.log(((1+Y**2)*(1+Z**2)/(1+Y**2+Z**2))*(((Y**2)*(1+Y**2+Z**2)/((1+Y**2)*(Y**2+Z**2)))**(Y**2))*
                                      (((Z**2)*(1+Z**2+Y**2)/((1+Z**2)*(Z**2+Y**2)))**(Z**2))))

    FF_floor_west = FF_floor_east

    y = wall_wide
    x = wall_long
    D = wall_high
    X = x/D
    Y = y/D

    FF_floor_ceiling = (2/(math.pi*X*Y))*((math.log((1+X**2)*(1+Y**2)/(1+X**2+Y**2)))**(1/2) + X*math.sqrt(1+Y**2)*math.atan(X/math.sqrt(1+Y**2)) +
                                          Y*math.sqrt(1+X**2)*math.atan(Y/math.sqrt(1+X**2)) - X*math.atan(X) - Y*math.atan(Y))

    B2_floor   = 0
    B2_south   = (1 - wall_absorption)*(FF_floor_south*Fraction_wall_S)*wall_absorption
    B2_north   = (1 - wall_absorption)*(FF_floor_north*Fraction_wall_N)*wall_absorption
    B2_east    = (1 - wall_absorption)*(FF_floor_east*Fraction_wall_E)*wall_absorption
    B2_west    = (1 - wall_absorption)*(FF_floor_west*Fraction_wall_W)*wall_absorption
    B2_ceiling = (1 - wall_absorption)*FF_floor_ceiling*wall_absorption
    
    if S_wind_area > 0:
        B2_window_S = (1 - wall_absorption)*(FF_floor_south*Fraction_window_S)*(1-(window_reflectance + N_int*int_window_absorption))
    else:
        B2_window_S = 0
    
    if N_wind_area > 0:
        B2_window_N = (1 - wall_absorption)*(FF_floor_north*Fraction_window_N)*(1-(window_reflectance + N_int*int_window_absorption))
    else:
        B2_window_N = 0
    
    if  E_wind_area > 0:
        B2_window_E = (1 - wall_absorption)*(FF_floor_east*Fraction_window_E)*(1-(window_reflectance + N_int*int_window_absorption))
    else:
        B2_window_E = 0
        
    if W_wind_area > 0:
        B2_window_W = (1 - wall_absorption)*(FF_floor_south*Fraction_window_W)*(1-(window_reflectance + N_int*int_window_absorption))
    else:
        B2_window_W = 0

    """
    ---------------------------------------------------------------------------------------------------------------
    B3 describes the third bounce such that the remaining nonabsorbed shortwave radiation is distributed over each
    surf in proportion to its area-absorptance product
    ---------------------------------------------------------------------------------------------------------------
    """
    Sum_B2 = (B2_floor + B2_south + B2_north + B2_east + B2_west + B2_ceiling + 
              B2_window_S +  B2_window_N + B2_window_E + B2_window_W)

    B3_floor   = (1 - wall_absorption - Sum_B2)*(floor_area/total_surf_area)*wall_absorption
    B3_south   = (1 - wall_absorption - Sum_B2)*(south_area/total_surf_area)*wall_absorption
    B3_east    = (1 - wall_absorption - Sum_B2)*(east_area/total_surf_area)*wall_absorption
    B3_west    = (1 - wall_absorption - Sum_B2)*(west_area/total_surf_area)*wall_absorption
    B3_north   = (1 - wall_absorption - Sum_B2)*(north_area/total_surf_area)*wall_absorption
    B3_ceiling = (1 - wall_absorption - Sum_B2)*(roof_area/total_surf_area)*wall_absorption
    
    if S_wind_area > 0:
        B3_window_S  = (1 - wall_absorption - Sum_B2)*(S_wind_area/total_surf_area)*(1-(window_reflectance + N_int*int_window_absorption))
    else:
        B3_window_S = 0
    
    if N_wind_area > 0:
        B3_window_N  = (1 - wall_absorption - Sum_B2)*(N_wind_area/total_surf_area)*(1-(window_reflectance + N_int*int_window_absorption))
    else:
        B3_window_N = 0
    
    if  E_wind_area > 0:
        B3_window_E  = (1 - wall_absorption - Sum_B2)*(E_wind_area/total_surf_area)*(1-(window_reflectance + N_int*int_window_absorption))
    else:
        B3_window_E = 0
        
    if W_wind_area > 0:
        B3_window_W  = (1 - wall_absorption - Sum_B2)*(W_wind_area/total_surf_area)*(1-(window_reflectance + N_int*int_window_absorption))
    else:
        B3_window_W = 0

    """
    ---------------------------------------------------------------------------------------------------------------
    BR describes the distribution of all remaining bounces based on distribution fractions from calculations for 
    B3n above
    ---------------------------------------------------------------------------------------------------------------
    """
    Sum_B3 = (B3_floor + B3_south + B3_north + B3_east + B3_west + B3_ceiling + 
              B3_window_S + B3_window_N + B3_window_E + B3_window_W)

    BR_floor   = (1 - wall_absorption - Sum_B2 - Sum_B3)*(B3_floor/Sum_B3)
    BR_south   = (1 - wall_absorption - Sum_B2 - Sum_B3)*(B3_south/Sum_B3)
    BR_east    = (1 - wall_absorption - Sum_B2 - Sum_B3)*(B3_east/Sum_B3)
    BR_west    = (1 - wall_absorption - Sum_B2 - Sum_B3)*(B3_west/Sum_B3)
    BR_north   = (1 - wall_absorption - Sum_B2 - Sum_B3)*(B3_north/Sum_B3)
    BR_ceiling = (1 - wall_absorption - Sum_B2 - Sum_B3)*(B3_ceiling/Sum_B3)


    SF_floor   = round(B1_floor + B2_floor + B3_floor + BR_floor,3)
    SF_south   = round(B1_south + B2_south + B3_south + BR_south,3)
    SF_east    = round(B1_east + B2_east + B3_east + BR_east,3)
    SF_west    = round(B1_west + B2_west + B3_west + BR_west,3)
    SF_north   = round(B1_north + B2_north + B3_north + BR_north,3)
    SF_ceiling = round(B1_ceiling + B2_ceiling + B3_ceiling + BR_ceiling,3)
        
    return(SF_south, SF_east, SF_west, SF_north, SF_ceiling, SF_floor)

"""
---------------------------------------------------------------------------------------------------------------
    The photovoltaic production profile is defined below by implementing the python pvlib library lines code. 
    -----------------------------------------------------------------------------------------------------------
    pvlib.py
    
    License
        BSD 3-Clause License
    Copyright
        Copyright (c) 2023 pvlib python Contributors
        Copyright (c) 2014 PVLIB python Development Team
        Copyright (c) 2013 Sandia National Laboratories
    
    Note from the author
        - The pvlib python Contributors comprises all authors of content of the pvlib python project. 
          A complete list of contributors can be found in the documentation.
        - The PVLIB python Development Team is the collection of developers of the pvlib python project. 
          Members of this team are included in the pvlib python Contributors.
        - Sandia National Laboratories originally developed pvlib python based on code in PVLib MATLAB.
---------------------------------------------------------------------------------------------------------------
"""
import pvlib
from pvlib import iotools, location
from pvlib.irradiance import get_total_irradiance
from pvlib.pvarray import pvefficiency_adr


def PV_load(rated_power,start,duration,timestep,Locality): 
    
    samples      = int(duration/timestep)
    stop         = int(start + duration)
    samples_h    = int((stop-start))    
    time_array         = np.linspace(start, stop, samples)   
    hours        = np.linspace(start, stop, samples_h)

    # Read a TMY3 file containing weather data and select needed columns
    #
    if Locality == 'Milano': 
        currentPath     = os.path.dirname("__file__")                
        weatherPath     = os.path.join(currentPath, 'epwFile')
        weatherFileName = "ITA_Milano-Linate.epw"
        weatherFile     = os.path.join(weatherPath, weatherFileName)
    if Locality == 'Roma': 
        currentPath     = os.path.dirname("__file__")                
        weatherPath     = os.path.join(currentPath, 'epwFile')
        weatherFileName = "ITA_Roma-Ciampino.epw"
        weatherFile     = os.path.join(weatherPath, weatherFileName)
    if Locality == 'Ancona': 
        currentPath     = os.path.dirname("__file__")                
        weatherPath     = os.path.join(currentPath, 'epwFile')
        weatherFileName = "ITA_Ancona.epw"
        weatherFile     = os.path.join(weatherPath, weatherFileName)
    if Locality == 'Messina': 
        currentPath     = os.path.dirname("__file__")                
        weatherPath     = os.path.join(currentPath, 'epwFile')
        weatherFileName = "ITA_Messina.epw"
        weatherFile     = os.path.join(weatherPath, weatherFileName)
    if Locality == 'Denver': # Denver International Airport
        currentPath     = os.path.dirname("__file__")                
        weatherPath     = os.path.join(currentPath, 'epwFile')
        weatherFileName = "USA_Denver.epw"
        weatherFile     = os.path.join(weatherPath, weatherFileName)
    
    epw, metadata = iotools.read_epw(weatherFile, coerce_year=None)

    df = pd.DataFrame({'ghi': epw['ghi'], 'dhi': epw['dhi'], 'dni': epw['dni'],
                       'temp_air': epw['temp_air'], 'wind_speed': epw['wind_speed'],
                       })

    # Shift timestamps to middle of hour and then calculate sun positions
    df.index = df.index - pd.Timedelta(minutes=30)

    #loc = location.Location.from_tmy(metadata)
    loc = location.Location.from_epw(metadata)
    solpos = loc.get_solarposition(df.index)

    # Determine  total irradiance on a fixed-tilt array
    TILT = metadata['latitude']
    ORIENT = 180

    total_irrad = get_total_irradiance(TILT, ORIENT,
                                       solpos.apparent_zenith, solpos.azimuth,
                                       df.dni, df.ghi, df.dhi)

    df['poa_global'] = total_irrad.poa_global

    # Estimate the expected operating temperature of the PV modules
    df['temp_pv'] = pvlib.temperature.faiman(df.poa_global, df.temp_air,
                                             df.wind_speed)

    # Now we're ready to calculate PV array DC output power based
    # on POA irradiance and PV module operating temperature.
    # Among the models available in pvlib-python to do this are:

    #  - PVWatts
    #  - SAPM
    #  - single-diode model variations

    # And now also the ADR PV efficiency model
    
    # Simulation is done in two steps:
    #  - first calculate efficiency using the ADR model,
    #  - then convert (scale up) efficiency to power.

    # ADR model parameters from the other example:
    adr_params = {'k_a': 0.999,
                  'k_d': -5.852,
                  'tc_d': 0.0194,
                  'k_rs': 0.06963,
                  'k_rsh': 0.2104
                  }

    df['eta_rel'] = pvefficiency_adr(df['poa_global'], df['temp_pv'], **adr_params)

    # Set the desired array size:
    P_STC = rated_power*1000.   # (W)

    # and the irradiance level needed to achieve this output:
    G_STC = 1000.   # (W/m2)

    
    df['p_mp'] = P_STC * df['eta_rel'] * (df['poa_global'] / G_STC)

    PV_ele     = np.interp(time_array, hours, (df['p_mp'][(start):(start+duration)])).tolist()

    return(PV_ele)
