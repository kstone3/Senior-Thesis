import sys
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pandas as pd
from matplotlib.animation import FuncAnimation
from IPython.core.display import HTML
import numpy as np
import math
from math import inf
from datetime import datetime,timedelta
import calendar
import warnings
import matplotlib
matplotlib.use("TkAgg")

# logfile = open("powellOutput.txt", "w")

# # 2) Redirect stdout (and stderr if you want all errors to go here, too)
# original_stdout = sys.stdout
# original_stderr = sys.stderr
# sys.stdout = logfile
# sys.stderr = logfile

class glacierSim():
    def __init__(self, ela=1880,ela_1900=1650,valley_length=3668, time=500,save=10,gamma=0.01,quiet=True, tune_factors=[-0.004,-0.002,0.0065,1.6,2.2], initial_ice=None, start_time=0):
        #MODEL VARS:
        self.valley_length = valley_length
        self.start_ela = ela #in m
        self.curr_ela=ela
        self.ela_1900=ela_1900
        self.num_cells = 50 #set number of cells
        self.dx = self.valley_length/(self.num_cells-1) #cell width in m
        if start_time==0: self.ice = np.zeros(self.num_cells) #initialize ice
        else: self.ice = np.array(initial_ice)
        
        #PLOTTING VARS:
        self.ice_line_list=[] #stores info for plotting ice
        self.snow_line_list=[] #stores info for plotting snow
        self.ela_line_list=[] #stores info for plotting ela
        self.title_list=[] #stores info for plotting title
        self.x = np.linspace(0.5 * self.dx, self.valley_length - 0.5 * self.dx,self.num_cells) #used for plotting topography and ice thickness
        self.topo =[] #initialize topography
        self.quiet=quiet #set to true to printout values while model is running
        
        #ICE FLUX VARS:
        self.q = np.zeros(self.num_cells+1) #initialize ice flux, need to have num_cells+1 to offset it from ice thickness for calculations
        self.g = 9.81 #gravity constant in m/yr^2
        self.p = 917 #density of ice
        self.ice_slope = np.zeros(self.num_cells, dtype=np.longdouble) #initialize ice_slope
        
        #TIME VARS:
        self.run_time=start_time*365.25 #DAYS
        self.prev_display=0 #used for displaying model data while running when quiet=false
        self.time = time #simulation time in years
        self.current_date=datetime(1484,1,1)+timedelta(days=start_time*365.25)
        self.save = save*365.25 #timestep interval in days
        self.frames = ((int)((self.time-start_time)/(self.save/365.25)))+1 #number of frames the animation will run for
        self.timestep_list=[] #days
        
        #MASS BALANCE VARS:
        self.b_max = float(-inf) #maximum yearly mass balance value for whole run
        self.b_min = float(inf) #minimum yearly mass balance valeu for whole run
        self.gamma = gamma #for mass balance equation
        self.b=np.zeros(self.num_cells) #initialize mass balance
        self.snow_depth=np.zeros(self.num_cells) #snow depth along glacier in m
        self.weather_dates=[] #dates for weather data
        self.temps=[] #daily temperatures
        self.precip=[] #daily precipitation
        
        #TUNE FACTORS:
        self.ice_melt_factor=tune_factors[0] #factor to change how much the ice melts per degree C
        self.snow_melt_factor=tune_factors[1] #factor to change how much the snow melts per degree C
        # self.accum_factor=tune_factors[2] #factor to change how much snow gets turned into ice
        # self.snow_conv_factor=tune_factors[3] #factor to change how much precip gets turned to snow
        # self.snow_melt_amplitude=tune_factors[4] #creates a curve to increase snow melt factor in the summer, peaks in august
        # self.ice_melt_amplitude=tune_factors[5] #creates a curve to increase ice melt factor in summer, peaks in august
        self.temp_lapse_rate=tune_factors[2] #temp lapse rate in C/m]
        self.accumfactor_lower=tune_factors[3] #lower bound to change how much snow gets turned into ice
        self.accumfactor_upper=tune_factors[4] #upper bound to change how much snow gets turned into ice
        
        #VERIF VARS:
        self.annual_mb=[] #annual mass balance verif data
        self.winter_mb=[] #winter (positive) mass balance verif data
        self.summer_mb=[] #summer (negative) mass balance verif data
        self.thickness_change_verif=np.zeros(4) #thickness change for verif
        self.front_variation_verif=np.zeros(100) #front variation verification data
        self.ela_verif=np.zeros(100) #ela's for verification
        self.runoff_verif=[] #holds volume validation data
        self.thickness_1986_verif=0 #thickness in 1986 for verif
        self.date_index=np.zeros(100) #stores dates used for vol change verif
        
        #CALCULATED VERIF VARS:
        self.volume_by_year=np.zeros(41) #volume change of glacier per year
        self.daily_volume_change=np.zeros(180) #tracks daily volume change of glacier, eventually summed up for monthly volume change
        self.calculated_annual_mb=np.zeros(40, dtype=np.float64) #used to verify annual mass balance
        self.calculated_winter_mb=np.zeros(40, dtype=np.float64) #used to verify winter (positive) mass balance
        self.calculated_summer_mb=np.zeros(40, dtype=np.float64) #used to verify summer (negative) mass balance
        self.front_variation_calc=np.zeros(100) #model calculated front variation
        self.thickness_change=np.zeros(4) #model thickness change
        self.year_mb=np.zeros(len(self.b)) #keeps track of the mb for the current year to calculate ela line
        self.ice_1986=np.zeros(self.num_cells) #stores the ice in 1986
        self.volume_change=[] #tracks the volume change of glacier
        self.ice_volume=0 #volume of glacier in m^3
        self.snow_melt_vol=0 #volume of melted snow in m^3, used to factor in snowmelt to volume calculations for stream flow data comparisons
        self.glacier_extent=0 #length of glacier in m
        self.ela_list=[] #list of ela values over time, used for verif
        self.ice_thickness_over_time=[] #tracks avg ice thickness over time
        
        #PREV VARS:
        self.prev_thickness=np.mean(self.ice) #previous avg ice thickness used to calculate thickness change
        self.prev_front=np.max(self.x[self.ice>1]) if np.any(self.ice>1) else 0 #previous front location to calculate front var change
        self.prev_volume=0 #previous glacier volume for calculating volume change
        self.initial_volume=0 #initial glacier volume for calculating volume change
        
        #WIDTH VARS:
        self.bins=[] #elevation bins for widths
        self.years=[] #years for width data
        self.areas=[] #area of glacier in m^2, used to calculate width
        self.widths=np.zeros(self.num_cells) #keeps track of glacier widths
        self.widths_over_time=[]  #used to plot width change over time
        
    def init(self, ax,ela=1880,ela_1900=1650,valley_length=3668, time=500,save=10,gamma=0.008, quiet=True, tune_factors=[-0.004,-0.002,0.0065,1.6,2.2], initial_ice=None, start_time=0):
        self.__init__(ela, ela_1900,valley_length, time, save, gamma, quiet, tune_factors, initial_ice, start_time)
        self.load_verif_data()
        self.calc_topo()
        self.calc_widths()
        self.load_mb_data()
        try: curr_ela=self.topo[np.where((self.b[:-1] >= 0) & (self.b[1:] < 0))[0]][0]
        except: curr_ela=self.topo[-1]
        self.ela_list.append(curr_ela)

    def haversine(self, lat1, lon1, lat2, lon2):
        a = math.sin(math.radians(lat2 - lat1) / 2) ** 2 + math.cos(math.radians(lat1)) * math.cos( math.radians(lat2)) * math.sin(math.radians(lon2 - lon1) / 2) ** 2
        return 6371000*(2 * math.atan2(math.sqrt(a), math.sqrt(1 - a)))
        
    def calc_topo(self):
        latitudes = []
        longitudes = []
        topo = []
        df = pd.read_csv('C:/Users/bookn/Downloads/Senior-Thesis/SIA/Data/centerlineBed.csv')
        latitudes = df.iloc[:, 2].astype(float).tolist()  # Latitude is the second column (index 2)
        longitudes = df.iloc[:, 1].astype(float).tolist()  # Longitude is the third column (index 1)
        topo = df.iloc[:, 0].astype(float).tolist()  # Elevation is the first column (index 0)
        cumulative_distances=[0.0]
        for i in range(1, len(latitudes)): cumulative_distances.append(cumulative_distances[-1] + self.haversine(latitudes[i - 1], longitudes[i - 1], latitudes[i], longitudes[i]) )
        self.x=np.linspace(cumulative_distances[0], cumulative_distances[-1], self.num_cells) #creates self.x values
        self.topo=np.interp(np.linspace(cumulative_distances[0], cumulative_distances[-1], self.num_cells), cumulative_distances, topo) #interpolates bed topo to get elevations at self.x values
        self.valley_length=np.max(self.x)
        self.dx = self.valley_length/(self.num_cells-1) #re-initializes dx using updates valley length
        self.ice_slope[:-1] = abs((np.diff(self.topo)/ self.dx)) #initial ice slope is just topo slope
        #THIS CREATES A FAKE INITIAL ICE FOR ICE FLUX TESTING
        #self.ice=np.exp(-0.5 * ((self.x - 700) / 300) ** 2)
        #scale=300/np.max(self.ice)
        #self.ice*=scale
        #INITIALIZE VOLUME VARIABLES
        self.ice_volume=np.sum(self.ice*self.widths*self.dx) 
        self.initial_volume=self.ice_volume
        self.prev_volume=self.ice_volume.copy()
        self.volume_change.append(self.initial_volume)
        #Not entirely sure why this is setup the way it is but it works
        #print(self.thickness_1986_verif.shape)
        df= pd.read_csv('C:/Users/bookn/Downloads/Senior-Thesis/SIA/Data/centerlineThickness.csv')
        self.thickness_1986_verif=df.iloc[0:, 3].astype(float).to_numpy()-df.iloc[0:, 0].astype(float).to_numpy()
        # latitudes = df.iloc[:, 2].astype(float).to_numpy()  # Latitude is the second column (index 2)
        # longitudes = df.iloc[:, 1].astype(float).to_numpy()  # Longitude is the third column (index 1)
        # if type(self.thickness_1986_verif) is not int:
        #     mask=~np.isnan(self.thickness_1986_verif)
        #     self.thickness_1986_verif=self.thickness_1986_verif[mask]
        #     latitudes=latitudes[mask]
        #     longitudes=longitudes[mask]
        # cumulative_distances=[0.0]
        # for i in range(1, len(latitudes)): cumulative_distances.append(cumulative_distances[-1] + self.haversine(latitudes[i - 1], longitudes[i - 1], latitudes[i], longitudes[i]) )
        self.thickness_1986_verif=np.interp(np.linspace(cumulative_distances[0], cumulative_distances[-1], self.num_cells), cumulative_distances, self.thickness_1986_verif) if type(self.thickness_1986_verif) is not int else np.zeros(self.num_cells)
        
    def calc_widths(self):
        df = pd.read_csv('C:/Users/bookn/Downloads/Senior-Thesis/SIA/Data/Input_SouthCascade_Area_Altitude_Distribution.csv')
        self.bins = df.columns[1:].astype(float).to_numpy()
        self.years = df.iloc[:, 0].astype(float).tolist()
        self.areas = df.iloc[:, 1:].astype(float).values*1000000
        self.bins=np.append(self.bins, self.bins[-1]+50)
        #I guess this thing throws some warnings, but it works fine
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            self.widths = np.array(self.areas[0,:]/np.abs(np.diff(self.x[np.array([np.argmin(np.abs(self.topo+self.ice - bin_value)) for bin_value in self.bins])])))[np.searchsorted(self.bins, (self.topo+self.ice))]
        self.widths[np.isinf(self.widths)]=0 #Due to the way widths are calculated, if the width is infinite it really should be 0
                
    def update_widths(self):
        prev_width=self.widths
        #Re-calculates widths
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            self.widths = np.array(self.areas[(self.current_date.year-1950),:]/np.abs(np.diff(self.x[np.array([np.argmin(np.abs(self.topo+self.ice - bin_value)) for bin_value in self.bins])])))[np.clip((np.searchsorted(self.bins, (self.topo+self.ice))),0,(len(self.areas[0])-1))]
        self.widths[np.isinf(self.widths)]=0
        if np.abs(self.widths-prev_width).any()>0:
            self.widths_over_time.append(self.widths.copy())
            self.ice_thickness_over_time.append(self.ice.copy())
        
    def load_mb_data(self):
        df = pd.read_csv('C:/Users/bookn/Downloads/Senior-Thesis/SIA/Data/Input_SouthCascade_Daily_Weather.csv')
        self.dates = pd.to_datetime(df.iloc[:, 0], format="%Y/%m/%d").tolist()
        self.temps = df.iloc[:, 1].astype(float).to_numpy()
        self.precip = df.iloc[:, 2].apply(lambda x: float(x) if not np.isnan(float(x)) else 0).to_numpy()
        #MAKE NOTE OF THIS IN THESIS, SOME OF THE TEMPERATURE DATA IS INTERPOLATED BECAUSE ITS MISSING
        nan_indices = np.isnan(self.temps)
        x = np.arange(len(self.temps))
        self.temps[nan_indices] = np.interp(x[nan_indices], x[~nan_indices], self.temps[~nan_indices])

    def load_verif_data(self):
        df = pd.read_csv('C:/Users/bookn/Downloads/Senior-Thesis/SIA/Data/Output_SouthCascade_Glacier_Wide_solutions_calibrated.csv', skiprows=25)
        self.annual_mb = df.iloc[:-1, 3].astype(float).tolist()
        self.summer_mb = df.iloc[:-1, 2].astype(float).tolist()
        self.winter_mb = df.iloc[:-1, 1].astype(float).tolist()
        self.ela_verif=df.iloc[5:, 4].astype(float).to_numpy()
        self.calculated_annual_mb=np.array([0] * len(self.annual_mb), dtype=np.float64)
        self.calculated_winter_mb=np.array([0] * len(self.winter_mb), dtype=np.float64)
        self.calculated_summer_mb=np.array([0] * len(self.summer_mb), dtype=np.float64)
        df = pd.read_csv('C:/Users/bookn/Downloads/Senior-Thesis/SIA/Data/daily_average_runoff_with_dates.csv')
        df['date'] = pd.to_datetime(df['date'])
        self.runoff_verif = df.groupby(df['date'].dt.to_period('M'))['runoff'].sum().to_numpy()
        self.date_index = {date: idx for idx, date in enumerate(df['date'])}
        self.daily_runoff=np.zeros(len(self.date_index))
        self.thickness_change_verif = pd.read_csv('C:/Users/bookn/Downloads/Senior-Thesis/SIA/Data/thickness_change.csv').iloc[0:, 11].astype(float).to_numpy()
        self.front_variation_verif = pd.read_csv('C:/Users/bookn/Downloads/Senior-Thesis/SIA/Data/front_variation_change.csv').iloc[0:, 9].astype(float).to_numpy()
        self.front_variation_calc = np.zeros(len(self.front_variation_verif))
        
    def calc_verif(self,timestep):
        #Round current date since all of the lists are based on day or year
        current_date_key = self.current_date.replace(hour=0, minute=0, second=0, microsecond=0)
        #These are the days before the start dates for the volume verif data so need to set the pre_volume to calculate volume change
        if current_date_key==datetime(2002,9,30) or  current_date_key==datetime(2003,6,8): self.prev_volume=self.ice_volume.copy()
        #If date is in the list of volume verification dates, calculate the volume change
        if current_date_key in self.date_index:
            self.daily_runoff[self.date_index[current_date_key]] += (self.prev_volume-self.ice_volume+self.snow_melt_vol)
            self.prev_volume=self.ice_volume.copy()
        #If date is the date before thickness change verification data starts then set prev_thickness to calculate thickness change
        if self.current_date==datetime(1998,12,31): self.prev_thickness=np.mean(self.ice)
        #Calculate thickness change data
        if 1999 <= self.current_date.year < 2004: 
            avg_thickness=np.mean(self.ice)
            self.thickness_change[0]+=avg_thickness-self.prev_thickness
            self.prev_thickness=avg_thickness
        elif 2004 <= self.current_date.year < 2009: 
            avg_thickness=np.mean(self.ice)
            self.thickness_change[1]+=avg_thickness-self.prev_thickness
            self.prev_thickness=avg_thickness
        elif 2009 <= self.current_date.year < 2014: 
            avg_thickness=np.mean(self.ice)
            self.thickness_change[2]+=avg_thickness-self.prev_thickness
            self.prev_thickness=avg_thickness
        if 2014 <= self.current_date.year < 2019: 
            avg_thickness=np.mean(self.ice)
            self.thickness_change[3]+=avg_thickness-self.prev_thickness
            self.prev_thickness=avg_thickness
        #If date is the date before front variation change verification data starts then set prev_front to calculate front variation change
        if self.current_date==datetime(1983,12,31): self.prev_front=np.max(self.x[self.ice > 1])
        #Calculate front variation data
        if 1984<=self.current_date.year<2009:
            if np.any(self.ice>1):
                self.front_variation_calc[int(self.current_date.year-1984)]+=self.prev_front-np.max(self.x[self.ice > 1])
                self.prev_front=np.max(self.x[self.ice > 1])
            else:
                self.front_variation_calc[int(self.current_date.year-1984)]+=self.prev_front-self.x[0]
                self.prev_front=self.x[0]
        #Calcualte ice thickness data for 1986
        if self.current_date.year==1986: self.ice_1986=self.ice[~np.isnan(self.thickness_1986_verif)].copy()
        #This will end up setting the yearly volume to the volume at the end of the year
        if self.current_date.year>=1984: self.volume_by_year[self.current_date.year-1984]=self.ice_volume
        #Add mass balance to verification arrays
        if 1984<=self.current_date.year<2024:
            date=int(self.current_date.year-1984)
            self.calculated_annual_mb[date]+=np.mean(np.array(self.b*timestep)) if self.b.size>0 else 0
            self.calculated_winter_mb[date]+=np.mean(np.array(self.b[self.b>0]*timestep)) if self.b[self.b>0].size>0 else 0
            self.calculated_summer_mb[date]+=np.mean(np.array(self.b[self.b<0]*timestep)) if self.b[self.b<0].size>0 else 0
        #Make sure mass balance verification arrays are correct
        if np.any(self.calculated_winter_mb<0): print("ERROR IN WINTER MB")
        if np.any(self.calculated_summer_mb>0): print("ERROR IN SUMMER MB")
    
    def snow_model(self, index, temps,timestep):
        #If temp is less than 0 precip falls as snow, so add it to snow depth
        #Note: This might need to be changed because if temp is just below 0, the snow_conv_factor will be wildly different then if the temp is way below 0
        new_snow=(self.precip[index]/1000)*self.snow_conv_factor*timestep
        self.snow_depth[temps<0]+=new_snow
        if np.any(new_snow<0):print("NEGATIVE NEW SNOW")
        #Calculates the amount of snow to melt in m based on temp and snow_melt_factor, value is always negative since temp is always negative
        snow_melt=(self.snow_melt_factor*temps[temps>0] if int(self.current_date.month) in [12,1,2] else (self.snow_melt_amplitude/2*(1-math.cos(2*math.pi/8*((self.current_date.month+(self.current_date.day+1)/calendar.monthrange(self.current_date.year, self.current_date.month)[1])-11))) + self.snow_melt_factor)*temps[temps>0])*timestep
        snow_melt_diff=(snow_melt * -1) - self.snow_depth[temps > 0] #Difference between the amount of snow to melt and the current snow depth
        snow_melt[snow_melt_diff>0] = self.snow_depth[temps>0][snow_melt_diff>0]*-1 #Janky way of setting the maximum snow melt to the snow depth
        if np.any(snow_melt>0): print("POSITIVE SNOW MELT")
        #Calculates the volume of melted snow
        #Note: Need to change width to valley width instead of glacier width
        self.snow_melt_vol=np.sum(snow_melt*self.dx*self.widths[temps>0])
        self.snow_depth[temps>0]+=snow_melt
        

    def update_b(self, timestep):
        if self.current_date>=datetime(1984,1,2):
            #Calculate which index to get weather data from
            if self.current_date<datetime(2024,10,1): index=self.dates.index(pd.Timestamp(self.current_date.replace(hour=0, minute=0, second=0, microsecond=0)))
            else: index=self.dates.index(pd.Timestamp(datetime(2024, 9, 30)))
            #Calculates temperatures for every self.x value and varies temp with elevation
            x_temps=self.temps[index]-self.temp_lapse_rate*(self.ice+self.topo-272) #weather station elevation is 272m
            mb=np.zeros_like(x_temps) #initialize mass balance
            #Melts ice for temps greater than 0
            mb[np.where((x_temps>0)&((self.ice+self.topo)>self.curr_ela))[0]]=self.snow_melt_factor*x_temps[np.where((x_temps>0)&((self.ice+self.topo)>self.curr_ela))[0]]
            mb[np.where((x_temps>0)&((self.ice+self.topo)<self.curr_ela))[0]]=self.ice_melt_factor*x_temps[np.where((x_temps>0)&((self.ice+self.topo)<self.curr_ela))[0]]
            #mb[x_temps>0]= self.ice_melt_factor*x_temps[x_temps>0] if int(self.current_date.month) in [12,1,2] else (self.ice_melt_amplitude/2*(1-math.cos(2*math.pi/8*((self.current_date.month+(self.current_date.day+1)/calendar.monthrange(self.current_date.year, self.current_date.month)[1])-11))) + self.ice_melt_factor)*x_temps[x_temps>0]
            #Calculates snow melt and accumulation
            # self.snow_model(index,x_temps,timestep)
            #Can't melt ice where there is snow
            #mb[(x_temps>0)&(self.snow_depth>0.1)]=0
            #Where temps are less than 0, accumulate ice based on snow depth
            #Note: Should this be changed to not factor in temps? Will snow accumulate to ice if the temp is above 0?
            # mb[x_temps<0]=self.snow_depth[x_temps<0]*self.accum_factor
            #IDEA: Try multiplying the accum_factor by elevation to change accumulation along glacier
            accumfactor = self.accumfactor_lower + ((self.current_date.year - 1984) / (2024 - 1984)) * (self.accumfactor_upper - self.accumfactor_lower)
            mb[x_temps<0]=(self.precip[index]/1000)*accumfactor
            #Subtract snow that was just turned into ice
            # subtract_snow = (self.snow_depth[x_temps<0]*self.accum_factor*timestep)
            # Check if any of the potential new values are negative
            # if np.any(subtract_snow > self.snow_depth[x_temps < 0]): print("THERE WILL BE NEGATIVE SNOW DEPTH")
            # self.snow_depth[x_temps<0]-=subtract_snow
            #Double check for accidental negative snow depth
            if np.any(self.snow_depth[x_temps<0]<0): print("NEGATIVE SNOW DEPTH, NEGATIVE TEMP")
            if np.any(self.snow_depth[x_temps>0]<0): print("NEGATIVE SNOW DEPTH, POSITIVE TEMP")
            #Make sure there is no negative snow, everything before this should make this impossible
            # self.snow_depth=np.maximum(self.snow_depth,0)
            #Alert for weird mass balance values
            if np.any(mb[x_temps<0]<0): print("NEGATIVE MB", mb[x_temps<0])
            if np.any(np.isnan(mb[x_temps < 0])): print("NAN MB", self.precip[index])
            if np.any(np.abs(mb)>100): print("MB TOO BIG",mb)
            return mb
        else: 
            #Used to spin up the glacier. Glacier starts retreating in 1900
            return ((self.topo+self.ice-self.curr_ela)*self.gamma)/365.25 #meters per day
            #return ((self.topo+self.ice-self.curr_ela)*self.gamma)/365.25
    
    def calc_q(self):
        #Calculate ice clops
        self.ice_slope[:-1] = -(np.diff(self.ice+self.topo) / self.dx)
        #Check for weird ice and ice_slope values
        if np.any(np.isnan(self.ice_slope)): #and not self.quiet:
            print('NaN detected in ice_slope:', self.ice_slope)
            print("MASS BALANCE: ", self.b)
            print("TIME: ", self.run_time)
            plt.plot(self.x, self.ice+self.topo)
            plt.plot(self.timestep_list)
            return
        if np.any(np.isnan(self.ice)): #and not self.quiet:
            print('NaN detected in ice:', self.ice)
            print("TIME: ", self.run_time)
            return
        if np.any(np.isinf(self.ice_slope)): # and not self.quiet:
            print('Infinity detected in ice_slope:', self.ice_slope)
            print('Ice: ',self.ice)
            print("Q: ",self.q)
            print("TIME: ", self.run_time)
            plt.plot(self.timestep_list)
            return
        #Calculate ice flux
        #self.q[1:] = (0.2 *(2e-17*2) * (self.p * self.g)**2) * np.sin(np.arctan(self.ice_slope))**2 * (self.ice**5)/5
        #self.q[1:]=2e-17* ((self.p*self.g*np.sin(np.arctan(self.ice_slope)))**3)*((self.ice**5)/5) #YEARS
        self.q[1:]=5.87e-19* ((self.p*self.g*np.sin(np.arctan(self.ice_slope)))**3)*((self.ice**5)/5) #DAYS
        #Check for weird q values
        if np.any(np.isnan(self.q)): #and not self.quiet:
            print('NaN detected in q:', self.q)
            print(self.ice_slope)
            print(self.ice)
            print("TIME: ", self.run_time)
            #plt.plot(self.x, self.ice+self.topo)
            plt.plot(self.timestep_list)
            return
        #Print model status every self.save years if quiet is false
        if (self.prev_display==0.0 or ((self.run_time>=(self.prev_display+(self.save/365.25))*365.25) and self.run_time<(self.time*365.25))) and not self.quiet:
            print("TIME: ", math.floor(self.run_time/365.25))
            print("ELA: ", self.curr_ela)
            print("ICE: ",self.ice)
            print("ICE VOLUME: ",self.ice_volume/1e9, " km^3")
            print("Difference from inital volume: ", self.ice_volume/1e9-self.initial_volume/1e9, " km^3")
            print("SLOPE: ",self.ice_slope)
            print("MASS BALANCE: ", self.b)
            print("Q: ",self.q)
            print("DQDX: ",np.diff(self.q)/self.dx)
            print("SUM Q: ", np.sum(self.q))
            print("SUM DQDX: ",np.sum(np.diff(self.q)/self.dx))
            print()
            self.prev_display=self.run_time/365.25
        #Calculate the change in ice thickness and return it
        return self.ice_slope,(self.b-(np.diff(self.q)/self.dx)) 

    def report_final_values(self,u):
        ice_extent = self.x[self.ice > 2]
        if ice_extent.size > 0:
            self.glacier_extent=np.max(ice_extent)
            print("Final glacier length: " + str(np.max(ice_extent)) + 'm')
        else: print("Final glacier length: 0m (no ice extent)")
        print("Final max ice thickness: " + str(np.max(self.ice)) + 'm')
        #print('Final max velocity: ' +str(np.max(u)) + "m/yr")
        # if(self.default_b):
        #     print('B min: ' + str(self.b_min))
        #     print('B max: ' + str(self.b_max))

    def run_model(self,i):
        #Used to plot initial values
        if i==0:
            curr_ela_plt=[self.start_ela]*self.num_cells
            self.ice_line_list.append((self.x, self.ice + self.topo, "c", "Ice"))
            self.snow_line_list.append((self.x, self.snow_depth + self.topo+self.ice, "c", "Snow"))
            self.ela_line_list.append((self.x, curr_ela_plt, "r", "dashed", "ELA"))
            self.title_list.append('Year = ' + str(self.current_date.year))
            return
        iter_time=0.0 #days
        timestep=0.0 #days
        while(iter_time<self.save):
            #End while loop if the difference between the save and current iter_time is less than 1, this is done to make sure it doesn't run over
            if (self.save-iter_time)<1: break
            ice_slope,dqdx=self.calc_q()
            #Calculate glacier velocity for timestep calculation
            u = (5.87e-19)*((self.p*self.g*np.sin(np.arctan(ice_slope)))**3)*(((self.ice)**4)/5)
            #Calculate the timestep based on glacier velocity, the 0.2 multiplier can be changed, increase it to increase the timestep
            #Timestep is constrained between 0.0001 and 1 and rouned to the 5th decimal place
            timestep = round(np.clip(((self.dx / np.max(u)) * 0.2), 0.0001, 1),5) if np.any(u > 0) else 1
            self.timestep_list.append(timestep)
            self.b=self.update_b(timestep)
            #Set new ice thickness
            self.ice = np.maximum((self.ice + dqdx * timestep), 0.0)
            self.ice_volume=np.sum(self.ice*self.dx*self.widths)
            #Update time variables
            iter_time+=float(timestep)
            self.run_time+=timestep
            self.current_date+=timedelta(days=float(timestep))
            #Calculate verification data
            self.calc_verif(timestep)
            #Update mass balance for current year
            self.year_mb+=self.b*timestep
            #Set new b_min and b_max
            self.b_max = max(np.max(self.b*365.25),self.b_max)
            self.b_min = min(np.min(self.b*365.25),self.b_min)
            #Update width if width data is available for the current date, updates on the 1st of the year
            if self.current_date>=datetime(1950,1,1) and self.current_date.day==1 and self.current_date.month==1: self.update_widths()
            if self.current_date.year==1900: self.curr_ela=self.ela_1900
        #If this if-statement runs somehow the model ran over which is a problem. The first if statement in the while loop should prevent this
        if(iter_time>self.save):
            print("BROKEN ABORT ABORT", iter_time,self.save)
            return
        #Calculate last ice_flux for this run
        ice_slope,dqdx=self.calc_q()
        #This if statement is also to make sure it doesn't run over. Since self.save should never be a fraction of a year then the current run should always end on 12/31
        if(self.current_date.day!=1):
            #Calculate the timestep so that it doesn't go over
            if(self.current_date.day==31): timestep=float((datetime((self.current_date.year+1),1,1)-self.current_date).total_seconds()/ (24 * 3600))
            else: timestep=float(self.save-iter_time)
            self.b=self.update_b(timestep)
            #Set new ice thickness
            self.ice = np.maximum((self.ice + dqdx * timestep), 0.0)
            self.ice_volume=np.sum(self.ice*self.dx*self.widths)
            #Update time variables
            iter_time+=float(timestep)
            self.run_time+=timestep
            self.current_date+=timedelta(days=float(timestep))
            #Calculate verification data
            self.calc_verif(timestep)
            #Update mass balance for current year. This actually ends up being the mb for self.save years
            self.year_mb+=self.b*timestep
            #Set new b_min and b_max
            self.b_max = max(np.max(self.b*365.25),self.b_max)
            self.b_min = min(np.min(self.b*365.25),self.b_min)
        #If the weather data is being used for the mass balance then calculate ela
        if self.current_date.year>=1984:
            #The try-except block is there incase the yearly mass balance is all positive or all negative
            try: curr_ela=float(self.topo[np.where((self.year_mb[:-1] >= 0) & (self.year_mb[1:] < 0))[0][0]])
            except:
                #If the try-except block fails then set ela to the top or bottom of the glacier depending on if the mass balance is all positive or all negative
                if np.all(self.year_mb>0): curr_ela=float(self.topo[0])
                else: curr_ela=float(self.topo[-1])
            #Add to ela list after 1989 because that's when ela verification data starts
            if self.current_date.year>1989: self.ela_list.append(curr_ela)
            #Reset year_mb to 0 for the next year
            self.year_mb.fill(0)
        #Otherwise set ela to the start ela
        else: curr_ela=self.curr_ela
        #Make sure to reset year_mb to 0 on the first of the year
        if self.current_date.month==1 and self.current_date.day==1: self.year_mb.fill(0)
        curr_ela_plt=[curr_ela]*self.num_cells
        #If the current date is the end date then report final values
        if(self.current_date.year==1484+self.time): self.report_final_values(u)
        #Update width if width data is available for the current date, updates on the 1st of the year
        if self.current_date>=datetime(1950,1,1) and self.current_date.day==1 and self.current_date.month==1: self.update_widths()
        #Add data to plotting lists
        self.ice_line_list.append((self.x, self.ice + self.topo, "c", "Ice"))
        self.snow_line_list.append((self.x, self.snow_depth + self.topo+self.ice, "c", "Snow"))
        self.ela_line_list.append((self.x, curr_ela_plt, "r", "dashed", "ELA"))
        self.title_list.append('Year = ' + str(self.current_date.year))

    #Used to plot the model data after the model has run
    def func(self,i):
        ax.clear() # Remove this to make lines stack up
        ax.set_ylim(min(self.topo) - 100, max(self.topo) + 100)
        ax.set_xlim(0, float(self.valley_length))
        ax.set_ylabel("Height (m)")
        ax.set_xlabel("Distance (m)")
        ax.set_aspect('equal', adjustable='datalim')
        ax.plot(self.x, self.topo, color="b", label="Topography")
        ax.set_title(self.title_list[i])
        self.ela_line=ax.plot(self.ela_line_list[i][0],self.ela_line_list[i][1], color=self.ela_line_list[i][2], linestyle=self.ela_line_list[i][3], label=self.ela_line_list[i][4])
        self.ice_line=ax.plot(self.ice_line_list[i][0],self.ice_line_list[i][1], color=self.ice_line_list[i][2], label=self.ice_line_list[i][3])
        self.snow_line=ax.plot(self.snow_line_list[i][0],self.snow_line_list[i][1], color=self.snow_line_list[i][2], label=self.snow_line_list[i][3])
        ax.legend()
        return self.ice_line_list,self.snow_line_list,self.ela_line_list
    

#Initialize plotting stuff outside of class
fig, ax = plt.subplots() #initialize plotting variables
_ = plt.close(fig) #used to prevent an empty plot from displaying

def optimize(parameter, input_params):
    try:
        # print("Optimizing: ", parameter)
        #print("ACCUMFACTOR: ", input_params[0])
        # print("ICE MELTFACTOR: ", input_params[0])
        # print("SNOW MELTFACTOR: ", input_params[1])
        # print("SNOW CONVERSION FACTOR: ", input_params[3])
        # print("SNOW MELT AMPLITUDE: ", input_params[4])
        # print("ICE MELT AMPLITUDE: ", input_params[5])
        # print("LAPSE RATE: ", input_params[1])
        print("ACCUMFACTOR LOWER: ", input_params[0])
        print("ACCUMFACTOR UPPER: ", input_params[1])
        print(input_params)
        #print("GAMMA: ", input_params[0])
        ela = 1880
        ela_1900=1910
        time = 540
        save = 10 #Needs to be 1 to calculate ela list
        gamma=0.016
        # accumfactor=0.1 #bounds approx 0.1-0.5
        # ice_meltfactor= 0.005 #bounds approx 0.005-0.012
        # snow_meltfactor=0.002 #bounds approx 0.002-0.006
        # snow_conv_factor=5 #bounds 5-15
        # snow_melt_amplitude=0.004
        # ice_melt_amplitude=0.007
        # tune_factors=[ice_meltfactor, snow_meltfactor, accumfactor, snow_conv_factor,snow_melt_amplitude,ice_melt_amplitude]
        #input_params=[-6.21530477e-03, -3.16373793e-03,  3.90806503e-03,  6.37238199e+00,-5.37087541e-03,-1.14816075e-03]
        input_params=[-0.003,-0.002,0.0065,input_params[0],input_params[1]]
        quiet = True
        start_time = 500
        ice = [ 53.89550985, 61.2302675, 68.52805603, 72.16752233, 78.19477103,
            86.57434438, 94.52278703, 113.00567764, 124.65045342, 131.1336047,
            132.61805723, 126.05975829, 117.01403765, 110.72201024, 107.36448442,
            113.23158091, 127.19316445, 150.80920841, 147.75800077, 135.19461921,
            132.62404257, 130.58089845, 128.95030221, 129.67226275, 131.1275745,
            137.76555425, 144.54125286, 149.5876901, 152.88274456, 152.52189309,
            142.71479768, 122.66643947, 105.65745331, 89.73132543, 84.64598526,
            84.39967405, 78.27781775, 70.24575207, 57.81117783, 43.58883439,
            34.64893057, 15.29705107, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
        model = glacierSim(ela=ela, ela_1900=ela_1900,time=time, save=save,gamma=gamma,quiet=quiet, tune_factors=input_params, initial_ice=ice, start_time=start_time)
        model.init(ax,ela=ela, ela_1900=ela_1900,time=time, save=save,gamma=gamma,quiet=quiet, tune_factors=input_params, initial_ice=ice, start_time=start_time)
        for i in range(0,model.frames):
            model.run_model(i)
        try:
            if parameter == 'summer':
                print("Function result: ",np.sum(np.abs(model.calculated_summer_mb - model.summer_mb)/np.abs(model.summer_mb))*100)
                return np.sum(np.abs(model.calculated_summer_mb - model.summer_mb)/np.abs(model.summer_mb))*100
            elif parameter == 'winter': 
                print("Function result: ",np.sum(np.abs(model.calculated_winter_mb-model.winter_mb)/model.winter_mb*100))
                return np.sum(np.abs(model.calculated_winter_mb-model.winter_mb)/model.winter_mb*100)
            elif parameter == 'annual': 
                print("Function result: ",np.sum((np.array(model.calculated_annual_mb) - np.array(model.annual_mb)) ** 2))
                return np.sum((np.array(model.calculated_annual_mb) - np.array(model.annual_mb)) ** 2)
            elif parameter == 'summer_winter':
                print("Function result: ",np.sum((np.array(model.calculated_summer_mb) - np.array(model.summer_mb)) ** 2) + np.sum((np.array(model.calculated_winter_mb) - np.array(model.winter_mb)) ** 2))
                return np.sum((np.array(model.calculated_summer_mb) - np.array(model.summer_mb)) ** 2) + np.sum((np.array(model.calculated_winter_mb) - np.array(model.winter_mb)) ** 2)
            elif parameter == 'ela':
                print("Function result: ", np.sum((model.ela_list-model.ela_verif)/np.abs(model.ela_verif)*100))
                return np.sum((model.ela_list-model.ela_verif)/np.abs(model.ela_verif)*100)
            elif parameter=='extent':
                print("Function result: ",abs(model.glacier_extent-3251)/3251)
                return abs(model.glacier_extent-3251)/3251
            elif parameter == 'front_var':
                print("Function result: ",np.sum((model.front_variation_calc-model.front_variation_verif) / np.abs(model.front_variation_verif) * 100))
                return np.sum((model.front_variation_calc-model.front_variation_verif) / np.abs(model.front_variation_verif) * 100)
            elif parameter == 'thick':
                print("Function result: ",np.sum((model.thickness_change-model.thickness_change_verif) /np.abs(model.thickness_change_verif) * 100))
                return np.sum((model.thickness_change-model.thickness_change_verif) /np.abs(model.thickness_change_verif) * 100)
            elif parameter == 'vol_change':
                df = pd.DataFrame({'date': pd.to_datetime(list(model.date_index.keys())),'volume_change': model.daily_volume_change})
                df['date'] = df['date'].dt.to_period('M') 
                monthly_volume_change = df.groupby('date')['volume_change'].sum().to_numpy()
                print("Function result: ",np.sum((monthly_volume_change-model.volume_valid) / np.abs(model.volume_valid)* 100))
                return np.sum((monthly_volume_change-model.volume_valid) / np.abs(model.volume_valid)* 100)
            else: raise ValueError("Invalid parameter. Choose from 'summer', 'winter', 'annual', 'summer_winter', 'vol_change'.")
        except Exception as e:
            print("Error during function calculation:", e) 
            print("Function result: INF")
            return inf
    except Exception as e:
        print("Error during model execution:", e) 
        print("Function result: INF")
        return inf

# Initialize the model
# accumfactor=0.001 #bounds approx 0.1-0.5, might want to change to 0.001,0.004
# ice_meltfactor= -0.005 #bounds approx -0.005--0.012
# snow_meltfactor=-0.002 #bounds approx -0.002--0.006
# snow_conv_factor=5 #bounds 5-15
# snow_melt_amplitude=-0.004
# ice_melt_amplitude=-0.007
#initial_guess=[-6.21530477e-03, -3.16373793e-03,  3.90806503e-03,  6.37238199e+00,-5.37087541e-03,-1.14816075e-03,0.0154]
accumfactor=0.01 #bounds approx 0.001-0.005
ice_meltfactor= -0.01 #bounds approx 0.005-0.012
snow_meltfactor=-0.01 #bounds approx 0.002-0.006
snow_conv_factor=6#bounds 5-15
snow_melt_amplitude=-0.001
ice_melt_amplitude=-0.001
#initial_guess=[ice_meltfactor, snow_meltfactor, accumfactor, snow_conv_factor,snow_melt_amplitude,ice_melt_amplitude]
#bounds = [(-0.012,-0.005), (-0.006,-0.002), (0.001, 0.005), (5, 15),(-0.01,-0.001),(-0.01,-0.001)]
#bounds=[(-0.02,-0.001),(-0.015,-0.001),(0.001,0.01),(5,10),(0,0),(0,0)]
# bounds=[(0.013,0.014)]
# bounds=[(-0.01,-0.00001),(0.0001,0.01)]
# gamma=0.0065
bounds=[(1.3,1.6),(1.5,2)]
# bounds=[(0.005,0.007)]
#result = differential_evolution(lambda x: optimize('summer_winter', x), bounds)
#initial_guess=[-9.57979633e-03, -9.85185254e-03 , 1.05581961e-02,  6.34612522e+00,-1.00242107e-03, -1.05457432e-03]
# initial_guess=[-0.012,-0.006,0.001,5,0,0]
# initial_guess=[-0.01,0.0065]
initial_guess=[1.5,1.9]
opt_method='Nelder-Mead'
with open(f"../Results/{opt_method}-Results.txt", "a") as file:
    file.write("-------------------Optimize Winter MB-----------\n")
    # for opt_var in ['ela','front_var','thick','vol_change']:
    for opt_var in ['winter']:
        result = minimize(lambda x: optimize(opt_var, x),initial_guess,method=opt_method,bounds=bounds,options={'disp': True})
        # result_ice_melt, result_snow_melt, result_accum, result_snow_conv, result_snow_amp, result_ice_amp=result.x
        result_accum=result.x
        # print("Final Gamma: ", result.x)
        # print("Final objective function value:", result.fun)
        # print("OPTIMIZED PARAMS FOR:", opt_var)
        print("Optimized accumulation factor:", result_accum)
        # print("Optimized ice meltfactor:",result_ice_melt)
        # print("Optimized snow meltfactor:",result_snow_melt)
        # print("Optimized snow conversion factor:",result_snow_conv)
        # print("Optimized snow melt amplitude:",result_snow_amp)
        # print("Optimized ice melt amplitude:",result_ice_amp)
        # print("Optimized lapse rate: ", result_lapse_rate)
        print("Final objective function value:", result.fun)
        print(result.x)
        # file.write(f"Optimized parameters for {opt_var}:\n")
        file.write(f"Optimized accumulation factor: {result_accum}\n")
        # file.write(f"Optimized ice meltfactor: {result_ice_melt}\n")
        # file.write(f"Optimized snow meltfactor: {result_snow_melt}\n")
        # file.write(f"Optimized snow conversion factor: {result_snow_conv}\n")
        # file.write(f"Optimized snow melt amplitude: {result_snow_amp}\n")
        # file.write(f"Optimized ice melt amplitude: {result_ice_amp}\n")
        file.write(f"Final objective function value: {result.fun}\n")
        # file.write(f"Optimized lapse rate: {result_lapse_rate}\n")
        file.write(f"{result.x}\n")
        # file.write(f"Final Gamma:  {result.x}\n")
        # file.write(f"Final objective function value: {result.fun}\n")
        file.flush()
# sys.stdout = original_stdout
# sys.stderr = original_stderr

# # 4) Close the file
# logfile.close()