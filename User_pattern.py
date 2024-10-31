"""
ClustEnergy OpTool
Released on January, 2024
@author: DIISM UNIVPM
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
#---------------------------------------------------------------------------------------------------------------------------------------
"""Occupancy profiles (occProfile), relative internal gains (intGains) and 
        setpoint temperature profiles for cooling (th_cooling) and heating (th_heating) spaces are "pre_defined" in this file.
        
        For modeling the occupational profiles, the richardsonpy tool was used. Whereas, temperature profiles were randomly 
        assigned according to a normal distribution (Tsetpoint_calc.py file in the folder). 
        However, these are modeling suggestions, which can be modified according to the users' needs in the ClustEnergy OpTool.
        
        Specifically, to model occupancy profiles, richardsonpy (RWTH-EBC, 2017) which is the Python implementation of the tool developed 
        by Richardson et al. (Richardson et al., 2008) was used. This tool generates stochastic occupancy profiles based on input parameters
        such as time step resolution and the number of occupants in a dwelling. Then, the generated occupancy model is used to create utility
        profiles for appliances and lights, facilitating the calculation of electricity demand. The household electricity demand profile 
        is established by incorporating time-use data from occupants (Richardson et al., 2010). Accordingly, the relative internal gains 
        of individual households, attributed to people, appliances, and lighting, are determined from the generated profiles. 
        This is done by following the methodology described in (RWTH-EBC, 2022). 
        For indoor air temperature profiles, a daily set-point temperature (Tsp) is associated with each dwelling. 
        This set-point is defined through a probabilistic normal distribution. Accordingly, by specifying the mean and standard deviation
        of the distribution, a daily set-point value is randomly assigned to each building.
        
        Occupancy and indoor gain profiles were generated considering a scenario with 4 occupants. 
        To define the temperature profiles, the design indoor air temperature for cooling (UNI/TS 11300-1, 2014) of 26.00 °C was considered 
        as the mean value and a standard deviation of 0.15 °C for normal distribution. In contrast, the indoor air design temperature for heating 
        (UNI/TS 11300-1, 2014) of 20.00 °C and a standard deviation of 0.15 °C for another normal distribution was considered for heating.
        
        REFERENCES
        RWTH-EBC, 2017. richarsonpy. RWTH-EBC. Accessed February 2023. https://github.com/RWTH-EBC/richardsonpy
        Richardson, I., Thomson, M. and Infield, D., 2008. A high-resolution domestic building occupancy model for energy demand simulations. Energy and Buildings. 40, 1560-1566.
        Richardson, I., Thomson, M., Infield, D. and Clifford, C., 2010. Domestic electricity use: A high-resolution energy demand model. Energy and Buildings. 42, 1878-1887.
        RWTH-EBC, 2022. districtgenerator. RWTH-EBC. Accessed February 2023. https://github.com/RWTH-EBC/districtgenerator
        UNI/TS 11300-1, 2014. Prestazioni energetiche degli edifici - Parte 1: Determinazione del fabbisogno di energia termica dell'edificio per la climatizzazione estiva ed invernale. UNI. Accessed February 2023
"""
#---------------------------------------------------------------------------------------------------------------------------------------
import os
import pandas as pd


class int_gains_nBui20(object): 
    def __init__(self, int_gains_FileName):
         int_gains_labels = ['time', 'n_Bui1', 'n_Bui2', 'n_Bui3', 'n_Bui4', 'n_Bui5', 'n_Bui6', 'n_Bui7', 'n_Bui8', 'n_Bui9',
                         'n_Bui10', 'n_Bui11', 'n_Bui12', 'n_Bui13', 'n_Bui14', 'n_Bui15', 'n_Bui16', 'n_Bui17', 'n_Bui18',
                         'n_Bui19', 'n_Bui20']
          
         # Import csv file
         self.int_gains_data = pd.read_csv(int_gains_FileName, skiprows=1, header=None, names=int_gains_labels, low_memory=(False))
         
    @property
    def gains_Bui1(self):        
         return self.int_gains_data['n_Bui1'].tolist()
     
    @property
    def gains_Bui2(self):        
          return self.int_gains_data['n_Bui2'].tolist()
    
    @property
    def gains_Bui3(self):        
          return self.int_gains_data['n_Bui3'].tolist()
      
    @property
    def gains_Bui4(self):        
          return self.int_gains_data['n_Bui4'].tolist()
      
    @property
    def gains_Bui5(self):        
         return self.int_gains_data['n_Bui5'].tolist()
     
    @property
    def gains_Bui6(self):        
          return self.int_gains_data['n_Bui6'].tolist()
    
    @property
    def gains_Bui7(self):        
          return self.int_gains_data['n_Bui7'].tolist()
      
    @property
    def gains_Bui8(self):        
          return self.int_gains_data['n_Bui8'].tolist()
      
    @property
    def gains_Bui9(self):        
         return self.int_gains_data['n_Bui9'].tolist()
     
    @property
    def gains_Bui10(self):        
          return self.int_gains_data['n_Bui10'].tolist()
    
    @property
    def gains_Bui11(self):        
          return self.int_gains_data['n_Bui11'].tolist()
      
    @property
    def gains_Bui12(self):        
          return self.int_gains_data['n_Bui12'].tolist()
      
    @property
    def gains_Bui13(self):        
         return self.int_gains_data['n_Bui13'].tolist()
     
    @property
    def gains_Bui14(self):        
          return self.int_gains_data['n_Bui14'].tolist()
    
    @property
    def gains_Bui15(self):        
          return self.int_gains_data['n_Bui15'].tolist()
      
    @property
    def gains_Bui16(self):        
          return self.int_gains_data['n_Bui16'].tolist()
      
    @property
    def gains_Bui17(self):        
          return self.int_gains_data['n_Bui17'].tolist()
      
    @property
    def gains_Bui18(self):        
         return self.int_gains_data['n_Bui18'].tolist()
     
    @property
    def gains_Bui19(self):        
          return self.int_gains_data['n_Bui19'].tolist()
    
    @property
    def gains_Bui20(self):        
          return self.int_gains_data['n_Bui20'].tolist()

def intGains(n_Bui,day,month,duration,timestep): 
    if n_Bui <= 20: 
        currentPath    = os.path.dirname("__file__")                
        UserPath       = os.path.join(currentPath, 'userProfiles')
        int_gains_FileName = "gains_profile.csv"
        int_g          = os.path.join(UserPath, int_gains_FileName)
        internal_g             = int_gains_nBui20(int_g) 

    if month == 1: 
        day_s = 0
    if month == 2:
        day_s = 744
    if month == 3: 
        day_s = 1416
    if month == 4:
        day_s = 2160     
    if month == 5:
        day_s = 2880
    if month == 6:
        day_s = 3624
    if month == 7: 
        day_s = 4344
    if month == 8: 
        day_s = 5088      
    if month == 9: 
        day_s = 5832
    if month == 10: 
        day_s = 6552
    if month == 11:
        day_s = 7296
    if month == 12:
        day_s = 8016      
    
    start        = int((day_s + (day-1)*24)/timestep)
    
    samples = int(duration/timestep)
    stop         = int(start + samples)

    int_gains_Bui1        = []
    int_gains_Bui2        = []
    int_gains_Bui3        = []
    int_gains_Bui4        = []
    int_gains_Bui5        = []
    int_gains_Bui6        = []
    int_gains_Bui7        = []
    int_gains_Bui8        = []
    int_gains_Bui9        = []
    int_gains_Bui10       = []
    int_gains_Bui11       = []
    int_gains_Bui12       = []
    int_gains_Bui13       = []
    int_gains_Bui14       = []
    int_gains_Bui15       = []
    int_gains_Bui16       = []
    int_gains_Bui17       = []
    int_gains_Bui18       = []
    int_gains_Bui19       = []
    int_gains_Bui20       = []
    
    for i in range(start,stop): 
        int_gains_Bui1.append(internal_g.gains_Bui1[i])
        int_gains_Bui2.append(internal_g.gains_Bui2[i])
        int_gains_Bui3.append(internal_g.gains_Bui3[i])
        int_gains_Bui4.append(internal_g.gains_Bui4[i])
        int_gains_Bui5.append(internal_g.gains_Bui5[i])
        int_gains_Bui6.append(internal_g.gains_Bui6[i])
        int_gains_Bui7.append(internal_g.gains_Bui7[i])
        int_gains_Bui8.append(internal_g.gains_Bui8[i])
        int_gains_Bui9.append(internal_g.gains_Bui9[i])
        int_gains_Bui10.append(internal_g.gains_Bui10[i])
        int_gains_Bui11.append(internal_g.gains_Bui11[i])
        int_gains_Bui12.append(internal_g.gains_Bui12[i])
        int_gains_Bui13.append(internal_g.gains_Bui13[i])
        int_gains_Bui14.append(internal_g.gains_Bui14[i])
        int_gains_Bui15.append(internal_g.gains_Bui15[i])
        int_gains_Bui16.append(internal_g.gains_Bui16[i])
        int_gains_Bui17.append(internal_g.gains_Bui17[i])
        int_gains_Bui18.append(internal_g.gains_Bui18[i])
        int_gains_Bui19.append(internal_g.gains_Bui19[i])
        int_gains_Bui20.append(internal_g.gains_Bui20[i])

    return(int_gains_Bui1,
           int_gains_Bui2,
           int_gains_Bui3,
           int_gains_Bui4,
           int_gains_Bui5,
           int_gains_Bui6,
           int_gains_Bui7,
           int_gains_Bui8,
           int_gains_Bui9,
           int_gains_Bui10,
           int_gains_Bui11,
           int_gains_Bui12,
           int_gains_Bui13,
           int_gains_Bui14,
           int_gains_Bui15,
           int_gains_Bui16,
           int_gains_Bui17,
           int_gains_Bui18,
           int_gains_Bui19,
           int_gains_Bui20)
    
class occ_profile_nBui20(object): 
    def __init__(self, th_sp_FileName):
         occ_labels = ['hour', 'n_Bui1', 'n_Bui2', 'n_Bui3', 'n_Bui4', 'n_Bui5', 'n_Bui6', 'n_Bui7', 'n_Bui8', 'n_Bui9',
                         'n_Bui10', 'n_Bui11', 'n_Bui12', 'n_Bui13', 'n_Bui14', 'n_Bui15', 'n_Bui16', 'n_Bui17', 'n_Bui18',
                         'n_Bui19', 'n_Bui20']
          
         # Import csv file
         self.occupancy_data = pd.read_csv(th_sp_FileName, skiprows=1, header=None, names=occ_labels, low_memory=(False))
         
    @property
    def occupancy_Bui1(self):        
         return self.occupancy_data['n_Bui1'].tolist()
     
    @property
    def occupancy_Bui2(self):        
          return self.occupancy_data['n_Bui2'].tolist()
    
    @property
    def occupancy_Bui3(self):        
          return self.occupancy_data['n_Bui3'].tolist()
      
    @property
    def occupancy_Bui4(self):        
          return self.occupancy_data['n_Bui4'].tolist()
      
    @property
    def occupancy_Bui5(self):        
         return self.occupancy_data['n_Bui5'].tolist()
     
    @property
    def occupancy_Bui6(self):        
          return self.occupancy_data['n_Bui6'].tolist()
    
    @property
    def occupancy_Bui7(self):        
          return self.occupancy_data['n_Bui7'].tolist()
      
    @property
    def occupancy_Bui8(self):        
          return self.occupancy_data['n_Bui8'].tolist()
      
    @property
    def occupancy_Bui9(self):        
         return self.occupancy_data['n_Bui9'].tolist()
     
    @property
    def occupancy_Bui10(self):        
          return self.occupancy_data['n_Bui10'].tolist()
    
    @property
    def occupancy_Bui11(self):        
          return self.occupancy_data['n_Bui11'].tolist()
      
    @property
    def occupancy_Bui12(self):        
          return self.occupancy_data['n_Bui12'].tolist()
      
    @property
    def occupancy_Bui13(self):        
         return self.occupancy_data['n_Bui13'].tolist()
     
    @property
    def occupancy_Bui14(self):        
          return self.occupancy_data['n_Bui14'].tolist()
    
    @property
    def occupancy_Bui15(self):        
          return self.occupancy_data['n_Bui15'].tolist()
      
    @property
    def occupancy_Bui16(self):        
          return self.occupancy_data['n_Bui16'].tolist()
      
    @property
    def occupancy_Bui17(self):        
          return self.occupancy_data['n_Bui17'].tolist()
      
    @property
    def occupancy_Bui18(self):        
         return self.occupancy_data['n_Bui18'].tolist()
     
    @property
    def occupancy_Bui19(self):        
          return self.occupancy_data['n_Bui19'].tolist()
    
    @property
    def occupancy_Bui20(self):        
          return self.occupancy_data['n_Bui20'].tolist()

def occProfile(n_Bui,day,month,duration,timestep): 
    if n_Bui <= 20: 
        currentPath    = os.path.dirname("__file__")                
        UserPath       = os.path.join(currentPath, 'userProfiles')
        occ_profile_FileName = "occ_profile.csv"
        occ_p          = os.path.join(UserPath, occ_profile_FileName)
        occ_profile             = occ_profile_nBui20(occ_p) 

    if month == 1: 
        day_s = 0
    if month == 2:
        day_s = 744
    if month == 3: 
        day_s = 1416
    if month == 4:
        day_s = 2160     
    if month == 5:
        day_s = 2880
    if month == 6:
        day_s = 3624
    if month == 7: 
        day_s = 4344
    if month == 8: 
        day_s = 5088      
    if month == 9: 
        day_s = 5832
    if month == 10: 
        day_s = 6552
    if month == 11:
        day_s = 7296
    if month == 12:
        day_s = 8016      
    
    start        = int((day_s + (day-1)*24)/timestep)
    
    samples = int(duration/timestep)
    stop         = int(start + samples)

    occ_profile_Bui1        = []
    occ_profile_Bui2        = []
    occ_profile_Bui3        = []
    occ_profile_Bui4        = []
    occ_profile_Bui5        = []
    occ_profile_Bui6        = []
    occ_profile_Bui7        = []
    occ_profile_Bui8        = []
    occ_profile_Bui9        = []
    occ_profile_Bui10       = []
    occ_profile_Bui11       = []
    occ_profile_Bui12       = []
    occ_profile_Bui13       = []
    occ_profile_Bui14       = []
    occ_profile_Bui15       = []
    occ_profile_Bui16       = []
    occ_profile_Bui17       = []
    occ_profile_Bui18       = []
    occ_profile_Bui19       = []
    occ_profile_Bui20       = []

    for i in range(start,stop): 
        occ_profile_Bui1.append(occ_profile.occupancy_Bui1[i])
        occ_profile_Bui2.append(occ_profile.occupancy_Bui2[i])
        occ_profile_Bui3.append(occ_profile.occupancy_Bui3[i])
        occ_profile_Bui4.append(occ_profile.occupancy_Bui4[i])
        occ_profile_Bui5.append(occ_profile.occupancy_Bui5[i])
        occ_profile_Bui6.append(occ_profile.occupancy_Bui6[i])
        occ_profile_Bui7.append(occ_profile.occupancy_Bui7[i])
        occ_profile_Bui8.append(occ_profile.occupancy_Bui8[i])
        occ_profile_Bui9.append(occ_profile.occupancy_Bui9[i])
        occ_profile_Bui10.append(occ_profile.occupancy_Bui10[i])
        occ_profile_Bui11.append(occ_profile.occupancy_Bui11[i])
        occ_profile_Bui12.append(occ_profile.occupancy_Bui12[i])
        occ_profile_Bui13.append(occ_profile.occupancy_Bui13[i])
        occ_profile_Bui14.append(occ_profile.occupancy_Bui14[i])
        occ_profile_Bui15.append(occ_profile.occupancy_Bui15[i])
        occ_profile_Bui16.append(occ_profile.occupancy_Bui16[i])
        occ_profile_Bui17.append(occ_profile.occupancy_Bui17[i])
        occ_profile_Bui18.append(occ_profile.occupancy_Bui18[i])
        occ_profile_Bui19.append(occ_profile.occupancy_Bui19[i])
        occ_profile_Bui20.append(occ_profile.occupancy_Bui20[i])

    return(occ_profile_Bui1,
           occ_profile_Bui2,
           occ_profile_Bui3,
           occ_profile_Bui4,
           occ_profile_Bui5,
           occ_profile_Bui6,
           occ_profile_Bui7,
           occ_profile_Bui8,
           occ_profile_Bui9,
           occ_profile_Bui10,
           occ_profile_Bui11,
           occ_profile_Bui12,
           occ_profile_Bui13,
           occ_profile_Bui14,
           occ_profile_Bui15,
           occ_profile_Bui16,
           occ_profile_Bui17,
           occ_profile_Bui18,
           occ_profile_Bui19,
           occ_profile_Bui20)

class Tsp_nBui20(object): 
    def __init__(self, th_sp_FileName):
         Tsp_labels = ['hour', 'n_Bui1', 'n_Bui2', 'n_Bui3', 'n_Bui4', 'n_Bui5', 'n_Bui6', 'n_Bui7', 'n_Bui8', 'n_Bui9',
                         'n_Bui10', 'n_Bui11', 'n_Bui12', 'n_Bui13', 'n_Bui14', 'n_Bui15', 'n_Bui16', 'n_Bui17', 'n_Bui18',
                         'n_Bui19', 'n_Bui20']
          
         # Import csv file
         self.Tsp_data = pd.read_csv(th_sp_FileName, skiprows=1, header=None, names=Tsp_labels, low_memory=(False))
         
    @property
    def Tsp_Bui1(self):        
         return self.Tsp_data['n_Bui1'].tolist()
     
    @property
    def Tsp_Bui2(self):        
          return self.Tsp_data['n_Bui2'].tolist()
    
    @property
    def Tsp_Bui3(self):        
          return self.Tsp_data['n_Bui3'].tolist()
      
    @property
    def Tsp_Bui4(self):        
          return self.Tsp_data['n_Bui4'].tolist()
      
    @property
    def Tsp_Bui5(self):        
         return self.Tsp_data['n_Bui5'].tolist()
     
    @property
    def Tsp_Bui6(self):        
          return self.Tsp_data['n_Bui6'].tolist()
    
    @property
    def Tsp_Bui7(self):        
          return self.Tsp_data['n_Bui7'].tolist()
      
    @property
    def Tsp_Bui8(self):        
          return self.Tsp_data['n_Bui8'].tolist()
      
    @property
    def Tsp_Bui9(self):        
         return self.Tsp_data['n_Bui9'].tolist()
     
    @property
    def Tsp_Bui10(self):        
          return self.Tsp_data['n_Bui10'].tolist()
    
    @property
    def Tsp_Bui11(self):        
          return self.Tsp_data['n_Bui11'].tolist()
      
    @property
    def Tsp_Bui12(self):        
          return self.Tsp_data['n_Bui12'].tolist()
      
    @property
    def Tsp_Bui13(self):        
         return self.Tsp_data['n_Bui13'].tolist()
     
    @property
    def Tsp_Bui14(self):        
          return self.Tsp_data['n_Bui14'].tolist()
    
    @property
    def Tsp_Bui15(self):        
          return self.Tsp_data['n_Bui15'].tolist()
      
    @property
    def Tsp_Bui16(self):        
          return self.Tsp_data['n_Bui16'].tolist()
      
    @property
    def Tsp_Bui17(self):        
          return self.Tsp_data['n_Bui17'].tolist()
      
    @property
    def Tsp_Bui18(self):        
         return self.Tsp_data['n_Bui18'].tolist()
     
    @property
    def Tsp_Bui19(self):        
          return self.Tsp_data['n_Bui19'].tolist()
    
    @property
    def Tsp_Bui20(self):        
          return self.Tsp_data['n_Bui20'].tolist()

def th_cooling(n_Bui,day,month,duration,timestep): 
    if n_Bui <= 20: 
        currentPath    = os.path.dirname("__file__")                
        UserPath       = os.path.join(currentPath, 'userProfiles')
        Thsp_FileName = "T_sp_cooling.csv"
        Thsp          = os.path.join(UserPath, Thsp_FileName)
        Tsp             = Tsp_nBui20(Thsp) 

    if month == 1: 
        day_s = 0
    if month == 2:
        day_s = 744
    if month == 3: 
        day_s = 1416
    if month == 4:
        day_s = 2160     
    if month == 5:
        day_s = 2880
    if month == 6:
        day_s = 3624
    if month == 7: 
        day_s = 4344
    if month == 8: 
        day_s = 5088      
    if month == 9: 
        day_s = 5832
    if month == 10: 
        day_s = 6552
    if month == 11:
        day_s = 7296
    if month == 12:
        day_s = 8016      
    
    start        = int((day_s + (day-1)*24)/timestep)
    
    samples = int(duration/timestep)
    stop         = int(start + samples)

    Tsp_Bui1        = []
    Tsp_Bui2        = []
    Tsp_Bui3        = []
    Tsp_Bui4        = []
    Tsp_Bui5        = []
    Tsp_Bui6        = []
    Tsp_Bui7        = []
    Tsp_Bui8        = []
    Tsp_Bui9        = []
    Tsp_Bui10       = []
    Tsp_Bui11       = []
    Tsp_Bui12       = []
    Tsp_Bui13       = []
    Tsp_Bui14       = []
    Tsp_Bui15       = []
    Tsp_Bui16       = []
    Tsp_Bui17       = []
    Tsp_Bui18       = []
    Tsp_Bui19       = []
    Tsp_Bui20       = []

    for i in range(start,stop): 
        Tsp_Bui1.append(Tsp.Tsp_Bui1[i])
        Tsp_Bui2.append(Tsp.Tsp_Bui2[i])
        Tsp_Bui3.append(Tsp.Tsp_Bui3[i])
        Tsp_Bui4.append(Tsp.Tsp_Bui4[i])
        Tsp_Bui5.append(Tsp.Tsp_Bui5[i])
        Tsp_Bui6.append(Tsp.Tsp_Bui6[i])
        Tsp_Bui7.append(Tsp.Tsp_Bui7[i])
        Tsp_Bui8.append(Tsp.Tsp_Bui8[i])
        Tsp_Bui9.append(Tsp.Tsp_Bui9[i])
        Tsp_Bui10.append(Tsp.Tsp_Bui10[i])
        Tsp_Bui11.append(Tsp.Tsp_Bui11[i])
        Tsp_Bui12.append(Tsp.Tsp_Bui12[i])
        Tsp_Bui13.append(Tsp.Tsp_Bui13[i])
        Tsp_Bui14.append(Tsp.Tsp_Bui14[i])
        Tsp_Bui15.append(Tsp.Tsp_Bui15[i])
        Tsp_Bui16.append(Tsp.Tsp_Bui16[i])
        Tsp_Bui17.append(Tsp.Tsp_Bui17[i])
        Tsp_Bui18.append(Tsp.Tsp_Bui18[i])
        Tsp_Bui19.append(Tsp.Tsp_Bui19[i])
        Tsp_Bui20.append(Tsp.Tsp_Bui20[i])

    return(Tsp_Bui1,
           Tsp_Bui2,
           Tsp_Bui3,
           Tsp_Bui4,
           Tsp_Bui5,
           Tsp_Bui6,
           Tsp_Bui7,
           Tsp_Bui8,
           Tsp_Bui9,
           Tsp_Bui10,
           Tsp_Bui11,
           Tsp_Bui12,
           Tsp_Bui13,
           Tsp_Bui14,
           Tsp_Bui15,
           Tsp_Bui16,
           Tsp_Bui17,
           Tsp_Bui18,
           Tsp_Bui19,
           Tsp_Bui20)

def th_heating(n_Bui,day,month,duration,timestep): 
    if n_Bui <= 20: 
        currentPath    = os.path.dirname("__file__")                
        UserPath       = os.path.join(currentPath, 'userProfiles')
        Thsp_FileName = "T_sp_heating.csv"
        Thsp          = os.path.join(UserPath, Thsp_FileName)
        Tsp             = Tsp_nBui20(Thsp) 

    if month == 1: 
        day_s = 0
    if month == 2:
        day_s = 744
    if month == 3: 
        day_s = 1416
    if month == 4:
        day_s = 2160     
    if month == 5:
        day_s = 2880
    if month == 6:
        day_s = 3624
    if month == 7: 
        day_s = 4344
    if month == 8: 
        day_s = 5088      
    if month == 9: 
        day_s = 5832
    if month == 10: 
        day_s = 6552
    if month == 11:
        day_s = 7296
    if month == 12:
        day_s = 8016      
    
    start        = int((day_s + (day-1)*24)/timestep)
    
    samples = int(duration/timestep)
    stop         = int(start + samples)

    Tsp_Bui1        = []
    Tsp_Bui2        = []
    Tsp_Bui3        = []
    Tsp_Bui4        = []
    Tsp_Bui5        = []
    Tsp_Bui6        = []
    Tsp_Bui7        = []
    Tsp_Bui8        = []
    Tsp_Bui9        = []
    Tsp_Bui10       = []
    Tsp_Bui11       = []
    Tsp_Bui12       = []
    Tsp_Bui13       = []
    Tsp_Bui14       = []
    Tsp_Bui15       = []
    Tsp_Bui16       = []
    Tsp_Bui17       = []
    Tsp_Bui18       = []
    Tsp_Bui19       = []
    Tsp_Bui20       = []

    for i in range(start,stop): 
        Tsp_Bui1.append(Tsp.Tsp_Bui1[i])
        Tsp_Bui2.append(Tsp.Tsp_Bui2[i])
        Tsp_Bui3.append(Tsp.Tsp_Bui3[i])
        Tsp_Bui4.append(Tsp.Tsp_Bui4[i])
        Tsp_Bui5.append(Tsp.Tsp_Bui5[i])
        Tsp_Bui6.append(Tsp.Tsp_Bui6[i])
        Tsp_Bui7.append(Tsp.Tsp_Bui7[i])
        Tsp_Bui8.append(Tsp.Tsp_Bui8[i])
        Tsp_Bui9.append(Tsp.Tsp_Bui9[i])
        Tsp_Bui10.append(Tsp.Tsp_Bui10[i])
        Tsp_Bui11.append(Tsp.Tsp_Bui11[i])
        Tsp_Bui12.append(Tsp.Tsp_Bui12[i])
        Tsp_Bui13.append(Tsp.Tsp_Bui13[i])
        Tsp_Bui14.append(Tsp.Tsp_Bui14[i])
        Tsp_Bui15.append(Tsp.Tsp_Bui15[i])
        Tsp_Bui16.append(Tsp.Tsp_Bui16[i])
        Tsp_Bui17.append(Tsp.Tsp_Bui17[i])
        Tsp_Bui18.append(Tsp.Tsp_Bui18[i])
        Tsp_Bui19.append(Tsp.Tsp_Bui19[i])
        Tsp_Bui20.append(Tsp.Tsp_Bui20[i])

    return(Tsp_Bui1,
           Tsp_Bui2,
           Tsp_Bui3,
           Tsp_Bui4,
           Tsp_Bui5,
           Tsp_Bui6,
           Tsp_Bui7,
           Tsp_Bui8,
           Tsp_Bui9,
           Tsp_Bui10,
           Tsp_Bui11,
           Tsp_Bui12,
           Tsp_Bui13,
           Tsp_Bui14,
           Tsp_Bui15,
           Tsp_Bui16,
           Tsp_Bui17,
           Tsp_Bui18,
           Tsp_Bui19,
           Tsp_Bui20)

class DHWdraw_nBui20(object): 
    def __init__(self, th_sp_FileName):
         Tsp_labels = ['n_Bui1', 'n_Bui2', 'n_Bui3', 'n_Bui4', 'n_Bui5', 'n_Bui6', 'n_Bui7', 'n_Bui8', 'n_Bui9',
                         'n_Bui10', 'n_Bui11', 'n_Bui12', 'n_Bui13', 'n_Bui14', 'n_Bui15', 'n_Bui16', 'n_Bui17', 'n_Bui18',
                         'n_Bui19', 'n_Bui20']
          
         # Import csv file
         self.Tsp_data = pd.read_csv(th_sp_FileName, skiprows=1, header=None, names=Tsp_labels, low_memory=(False))
         
    @property
    def DHWdraw_Bui1(self):        
         return self.Tsp_data['n_Bui1'].tolist()
     
    @property
    def DHWdraw_Bui2(self):        
          return self.Tsp_data['n_Bui2'].tolist()
    
    @property
    def DHWdraw_Bui3(self):        
          return self.Tsp_data['n_Bui3'].tolist()
      
    @property
    def DHWdraw_Bui4(self):        
          return self.Tsp_data['n_Bui4'].tolist()
      
    @property
    def DHWdraw_Bui5(self):        
         return self.Tsp_data['n_Bui5'].tolist()
     
    @property
    def DHWdraw_Bui6(self):        
          return self.Tsp_data['n_Bui6'].tolist()
    
    @property
    def DHWdraw_Bui7(self):        
          return self.Tsp_data['n_Bui7'].tolist()
      
    @property
    def DHWdraw_Bui8(self):        
          return self.Tsp_data['n_Bui8'].tolist()
      
    @property
    def DHWdraw_Bui9(self):        
         return self.Tsp_data['n_Bui9'].tolist()
     
    @property
    def DHWdraw_Bui10(self):        
          return self.Tsp_data['n_Bui10'].tolist()
    
    @property
    def DHWdraw_Bui11(self):        
          return self.Tsp_data['n_Bui11'].tolist()
      
    @property
    def DHWdraw_Bui12(self):        
          return self.Tsp_data['n_Bui12'].tolist()
      
    @property
    def DHWdraw_Bui13(self):        
         return self.Tsp_data['n_Bui13'].tolist()
     
    @property
    def DHWdraw_Bui14(self):        
          return self.Tsp_data['n_Bui14'].tolist()
    
    @property
    def DHWdraw_Bui15(self):        
          return self.Tsp_data['n_Bui15'].tolist()
      
    @property
    def DHWdraw_Bui16(self):        
          return self.Tsp_data['n_Bui16'].tolist()
      
    @property
    def DHWdraw_Bui17(self):        
          return self.Tsp_data['n_Bui17'].tolist()
      
    @property
    def DHWdraw_Bui18(self):        
         return self.Tsp_data['n_Bui18'].tolist()
     
    @property
    def DHWdraw_Bui19(self):        
          return self.Tsp_data['n_Bui19'].tolist()
    
    @property
    def DHWdraw_Bui20(self):        
          return self.Tsp_data['n_Bui20'].tolist()

def drawProfile(n_Bui,day,month,duration,timestep): 
    mins = int(timestep*60)
    if n_Bui <= 20: 
        currentPath    = os.path.dirname("__file__")                
        UserPath       = os.path.join(currentPath, 'userProfiles')
        draw_FileName = 'DHWdraw_profile' + str(mins) +'mins_step.csv'
        draw          = os.path.join(UserPath, draw_FileName)
        profile             = DHWdraw_nBui20(draw) 

    if month == 1: 
        day_s = 0
    if month == 2:
        day_s = 744
    if month == 3: 
        day_s = 1416
    if month == 4:
        day_s = 2160     
    if month == 5:
        day_s = 2880
    if month == 6:
        day_s = 3624
    if month == 7: 
        day_s = 4344
    if month == 8: 
        day_s = 5088      
    if month == 9: 
        day_s = 5832
    if month == 10: 
        day_s = 6552
    if month == 11:
        day_s = 7296
    if month == 12:
        day_s = 8016      
    
    start        = int((day_s + (day-1)*24)/timestep)
    
    samples = int(duration/timestep)
    stop         = int(start + samples)

    draw_Bui1        = []
    draw_Bui2        = []
    draw_Bui3        = []
    draw_Bui4        = []
    draw_Bui5        = []
    draw_Bui6        = []
    draw_Bui7        = []
    draw_Bui8        = []
    draw_Bui9        = []
    draw_Bui10       = []
    draw_Bui11       = []
    draw_Bui12       = []
    draw_Bui13       = []
    draw_Bui14       = []
    draw_Bui15       = []
    draw_Bui16       = []
    draw_Bui17       = []
    draw_Bui18       = []
    draw_Bui19       = []
    draw_Bui20       = []

    for i in range(start,stop): 
        draw_Bui1.append(profile.DHWdraw_Bui1[i])
        draw_Bui2.append(profile.DHWdraw_Bui2[i])
        draw_Bui3.append(profile.DHWdraw_Bui3[i])
        draw_Bui4.append(profile.DHWdraw_Bui4[i])
        draw_Bui5.append(profile.DHWdraw_Bui5[i])
        draw_Bui6.append(profile.DHWdraw_Bui6[i])
        draw_Bui7.append(profile.DHWdraw_Bui7[i])
        draw_Bui8.append(profile.DHWdraw_Bui8[i])
        draw_Bui9.append(profile.DHWdraw_Bui9[i])
        draw_Bui10.append(profile.DHWdraw_Bui10[i])
        draw_Bui11.append(profile.DHWdraw_Bui11[i])
        draw_Bui12.append(profile.DHWdraw_Bui12[i])
        draw_Bui13.append(profile.DHWdraw_Bui13[i])
        draw_Bui14.append(profile.DHWdraw_Bui14[i])
        draw_Bui15.append(profile.DHWdraw_Bui15[i])
        draw_Bui16.append(profile.DHWdraw_Bui16[i])
        draw_Bui17.append(profile.DHWdraw_Bui17[i])
        draw_Bui18.append(profile.DHWdraw_Bui18[i])
        draw_Bui19.append(profile.DHWdraw_Bui19[i])
        draw_Bui20.append(profile.DHWdraw_Bui20[i])

    return(draw_Bui1,
           draw_Bui2,
           draw_Bui3,
           draw_Bui4,
           draw_Bui5,
           draw_Bui6,
           draw_Bui7,
           draw_Bui8,
           draw_Bui9,
           draw_Bui10,
           draw_Bui11,
           draw_Bui12,
           draw_Bui13,
           draw_Bui14,
           draw_Bui15,
           draw_Bui16,
           draw_Bui17,
           draw_Bui18,
           draw_Bui19,
           draw_Bui20)