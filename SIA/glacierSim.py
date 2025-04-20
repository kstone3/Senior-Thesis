import pandas as pd
import numpy as np
import math
from math import inf
from datetime import datetime,timedelta
import random
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

class glacierSim():
    def __init__(self, ela=1880,ela_1900=1650,time=500,save=10,gamma=0.01,quiet=True, tune_factors=[-0.004,-0.002,0.0065,1.6,2.2,5,0.006], initial_ice=None, start_time=0, input_files=None):
        #MODEL VARS:
        self.valley_length=0 #in m, defined in calc_topo
        self.start_ela = ela #in m
        self.curr_ela=ela
        self.ela_1900=ela_1900
        self.num_cells = 50 #set number of cells
        self.dx = 0 #cell width in m, will be calculated in calc_topo
        if start_time==0: self.ice = np.zeros(self.num_cells) #initialize ice
        else: self.ice = np.array(initial_ice)
        
        #PLOTTING VARS:
        self.ice_line_list=[] #stores info for plotting ice
        self.snow_line_list=[] #stores info for plotting snow
        self.ela_line_list=[] #stores info for plotting ela
        self.title_list=[] #stores info for plotting title
        self.x = np.zeros(self.num_cells) #used for plotting topography and ice thickness, defined in calc_topo
        self.topo =[] #initialize topography
        self.quiet=quiet #set to true to printout values while model is running
        
        #ICE FLUX VARS:
        self.q = np.zeros(self.num_cells+1) #initialize ice flux, need to have num_cells+1 to offset it from ice thickness for calculations
        self.g = 9.81 #gravity constant in m/yr^2
        self.p = 917 #density of ice
        self.ice_slope = np.zeros(self.num_cells) #initialize ice_slope
        
        #TIME VARS:
        self.time=start_time*365.25 #DAYS
        self.prev_display=0 #used for displaying model data while running when quiet=false
        self.time = time #simulation time in years
        # self.current_date=datetime(1484,1,1)+timedelta(days=start_time*365.25)
        if start_time==0: self.current_date=datetime(1484,1,1)
        else: self.current_date=datetime(1984,1,1)
        self.save = save*365.25 #timestep interval in days
        self.frames = ((int)((self.time-start_time)/(self.save/365.25)))+1 #number of frames the animation will run for
        self.timestep_list=[] #days
        self.month_count=np.zeros((self.time-start_time)*12+1)

        #MASS BALANCE VARS:
        self.b_max = float(-inf) #maximum yearly mass balance value for whole run
        self.b_min = float(inf) #minimum yearly mass balance valeu for whole run
        self.gamma = gamma #for mass balance equation
        self.b=np.zeros(self.num_cells) #initialize mass balance
        self.weather_dates=[] #dates for weather data
        self.temps=[] #daily temperatures
        self.precip=[] #daily precipitation
        
        #SNOW VARS:
        self.snow_depth=np.zeros(self.num_cells) #snow depth along glacier in m
        self.snow_melt_amt=np.zeros(self.num_cells) #snow melt amount in m
        self.monthly_precip_vol=np.zeros((self.time-start_time)*12+1) #holds daily snow volume data for all data
        self.snow_depth_list=[]
        self.avalanche_dates=[self.current_date]*4
        self.monthly_snow_depth=np.zeros((self.time-start_time)*12+1)
        
        #TUNE FACTORS:
        self.ice_melt_factor=tune_factors[0] #factor to change how much the ice melts per degree C
        self.snow_melt_factor=tune_factors[1] #factor to change how much the snow melts per degree C
        self.temp_lapse_rate=tune_factors[2] #temp lapse rate in C/m]
        self.accumfactor_lower=tune_factors[3] #lower bound to change how much snow gets turned into ice
        self.accumfactor_upper=tune_factors[4] #upper bound to change how much snow gets turned into ice
        self.avalanche_percent=tune_factors[5] #factor to change how much precip gets turned to snow
        self.precip_conv_factor=tune_factors[6] #precipitation lapse rate in m/m
        
        #VERIF VARS:
        #MB VARS
        self.annual_mb=[] #annual mass balance verif data
        self.winter_mb=[] #winter (positive) mass balance verif data
        self.summer_mb=[] #summer (negative) mass balance verif data
        #GLAICER VARS
        self.thickness_change_verif=np.zeros(4) #thickness change for verif
        self.front_variation_verif=np.zeros(100) #front variation verification data
        self.ela_verif=np.zeros(100) #ela's for verification
        self.thickness_1958_verif=0 #thickness in 1958 for verif
        self.thickness_1986_verif=0 #thickness in 1986 for verif
        self.thickness_2021_verif=0 #thickness in 2021 for verif
        #RUNOFF VARS
        self.measured_runoff_training=[] #holds daily runoff data for training
        self.measured_runoff_verif=[] #holds daily runoff data for verification
        self.measured_runoff_all=[]
        self.date_index_training=[] #holds date index for training data
        self.date_index_verif=[] #holds date index for verification data
        self.date_index_all=[] #holds date index for all data
        self.month_index_all=[]
        
        #CALCULATED VERIF VARS:
        #MB VARS
        self.calculated_annual_mb=np.zeros(40, dtype=np.float64) #used to verify annual mass balance
        self.calculated_winter_mb=np.zeros(40, dtype=np.float64) #used to verify winter (positive) mass balance
        self.calculated_summer_mb=np.zeros(40, dtype=np.float64) #used to verify summer (negative) mass balance
        self.year_mb=np.zeros(len(self.b)) #keeps track of the mb for the current year to calculate ela line
        #GLACIER VARS
        self.front_variation_calc=np.zeros(100) #model calculated front variation
        self.thickness_change=np.zeros(4) #model thickness change
        self.ice_1958=np.zeros(self.num_cells) #stores the ice in 1958
        self.ice_1986=np.zeros(self.num_cells) #stores the ice in 1986
        self.ice_2021=np.zeros(self.num_cells) #stores the ice in 2021
        self.glacier_extent=0 #length of glacier in m
        self.ice_thickness_over_time=[] #tracks avg ice thickness over time
        self.ela_list=[] #list of ela values over time, used for verif
        #RUNOFF VARS
        self.glacial_melt=0
        self.rain_vol_per_step=0
        self.snow_melt_vol=0 #volume of melted snow in m^3, used to factor in snowmelt to volume calculations for stream flow data comparisons
        self.daily_runoff_training=[] #holds daily runoff data for training
        self.daily_runoff_verif=[] #holds daily runoff data for verification
        self.daily_runoff_all_data=[]
        self.monthly_runoff_all=np.zeros((self.time-start_time)*12+1) #holds daily runoff data for all data
        self.monthly_glacier_melt=np.zeros((self.time-start_time)*12+1) #holds daily runoff data for all data
        self.precip_accum=np.zeros(41)
        self.glacier_melt_runoff_data=[]
        
        #PREV VARS:
        self.prev_thickness=np.mean(self.ice) #previous avg ice thickness used to calculate thickness change
        self.prev_front=np.max(self.x[self.ice>1]) if np.any(self.ice>1) else 0 #previous front location to calculate front var change
        
        #WIDTH VARS: 
        self.bins=[] #elevation bins for widths
        self.years=[] #years for width data
        self.areas=[] #area of glacier in m^2, used to calculate width
        self.year_area=np.zeros(self.num_cells)
        self.glacier_area=0
        
        #INPUT FILE VARS:
        if input_files is not None: self.input_files=input_files #list of input files for glacier data
        else: print("No input files provided")
        
    def init(self, ela=1880,ela_1900=1650,valley_length=3668, time=500,save=10,gamma=0.008, quiet=True, tune_factors=[-0.004,-0.002,0.0065,1.6,2.2,5,0.006], initial_ice=None, start_time=0, input_files=None): 
        self.__init__(ela, ela_1900, time, save, gamma, quiet, tune_factors, initial_ice, start_time, input_files)
        self.load_verif_data()
        self.calc_topo()
        self.load_area_data()
        self.load_weather_data()
        self.ela_list.append(self.topo[-1]) 
    
    #Used to calculate distance between two lat lon points to interpolate bed topo
    def haversine(self, lat1, lon1, lat2, lon2):
        a = math.sin(math.radians(lat2 - lat1) / 2) ** 2 + math.cos(math.radians(lat1)) * math.cos( math.radians(lat2)) * math.sin(math.radians(lon2 - lon1) / 2) ** 2
        return 6371000*(2 * math.atan2(math.sqrt(a), math.sqrt(1 - a)))
        
    def calc_topo(self):
        latitudes = []
        longitudes = []
        topo = []
        #Load in the bed topo data
        df = pd.read_csv(self.input_files['bed'])
        latitudes = df.iloc[:, 2].astype(float).tolist()  # Latitude is the second column (index 2)
        longitudes = df.iloc[:, 1].astype(float).tolist()  # Longitude is the third column (index 1)
        topo = df.iloc[:, 0].astype(float).tolist()  # Elevation is the first column (index 0)
        cumulative_distances=[0.0]
        #Interpolate the bed topo data into my topo array with a set number of cells
        for i in range(1, len(latitudes)): cumulative_distances.append(cumulative_distances[-1] + self.haversine(latitudes[i - 1], longitudes[i - 1], latitudes[i], longitudes[i]) )
        self.x=np.linspace(cumulative_distances[0], cumulative_distances[-1], self.num_cells) #creates self.x values
        self.topo=np.interp(np.linspace(cumulative_distances[0], cumulative_distances[-1], self.num_cells), cumulative_distances, topo) #interpolates bed topo to get elevations at self.x values
        #Reinitialize the valley length, dx and ice slope
        self.valley_length=np.max(self.x)
        self.dx = self.valley_length/(self.num_cells-1) #re-initializes dx using updates valley length
        self.ice_slope[:-1] = abs((np.diff(self.topo)/ self.dx)) #initial ice slope is just topo slope
        #THIS CREATES A FAKE INITIAL ICE FOR ICE FLUX TESTING
        #self.ice=np.exp(-0.5 * ((self.x - 700) / 300) ** 2)
        #scale=300/np.max(self.ice)
        #self.ice*=scale
        #Load in the 2021 and 1986 glacier extent data for verification
        #This is done here instead of in the load_verif function below because 
        #it needs the cumulative_distances to interpolate the data onto the model arrays
        data_1958 = np.loadtxt(self.input_files['glacier_1958'], delimiter=',', skiprows=1)
        self.thickness_1958_verif = data_1958[:, 3] - data_1958[:, 0]
        data_1986 = np.loadtxt(self.input_files['glacier_1986'], delimiter=',', skiprows=1)
        self.thickness_1986_verif = data_1986[:, 3] - data_1986[:, 0]
        data_2021 = np.loadtxt(self.input_files['glacier_2021'], delimiter=',', skiprows=1)
        self.thickness_2021_verif = data_2021[:, 3] - data_2021[:, 0]
        self.thickness_1958_verif=(np.interp(np.linspace(cumulative_distances[0], cumulative_distances[-1], self.num_cells), cumulative_distances, self.thickness_1958_verif) if type(self.thickness_1958_verif) is not int else np.zeros(self.num_cells))
        self.thickness_1958_verif = np.maximum(self.thickness_1958_verif, 0)  # Ensure no negative thickness
        self.thickness_1986_verif=np.interp(np.linspace(cumulative_distances[0], cumulative_distances[-1], self.num_cells), cumulative_distances, self.thickness_1986_verif) if type(self.thickness_1986_verif) is not int else np.zeros(self.num_cells)
        self.thickness_1986_verif = np.maximum(self.thickness_1986_verif, 0)  # Ensure no negative thickness
        self.thickness_2021_verif=np.maximum(np.interp(np.linspace(cumulative_distances[0], cumulative_distances[-1], self.num_cells), cumulative_distances, self.thickness_2021_verif) if type(self.thickness_2021_verif) is not int else np.zeros(self.num_cells),0)
        self.thickness_2021_verif = np.maximum(self.thickness_2021_verif, 0)  # Ensure no negative thickness
        
    #Calculate the initial glacier area
    def load_area_data(self):
        #Load in the area data
        df = pd.read_csv(self.input_files['area'])
        self.bins = df.columns[2:].astype(float).to_numpy()
        self.years = df.iloc[:, 0].astype(float).tolist()
        self.areas = df.iloc[:, 2:].astype(float).to_numpy()*1000000 #convert to m^2
        self.bin_bounds=np.concatenate(([self.bins[0] - 25], (self.bins[:-1] + self.bins[1:]) / 2, [self.bins[-1] + 25]))
        #Call update_areas() to initialize the area arrays
        if self.current_date.year>=1984:self.update_areas()
                
    def update_areas(self):
        #Get the area for the current year, area data starts in 1950
        area=self.areas[(self.current_date.year-1950),:]
        bin_indices = np.digitize((self.topo+self.ice), self.bins)
        for i in range(len(area)):
            mask = (bin_indices == i)
            cell_indices = np.where(mask)[0]
            n_cells = len(cell_indices)
            if n_cells > 0:
                base_area = area[i] / n_cells
                self.year_area[mask] = base_area
                self.year_area[cell_indices[0]] += (area[i] - (base_area * n_cells))
        if np.sum(self.year_area)<np.sum(area): self.year_area+=(np.sum(area)-np.sum(self.year_area))/len(self.year_area)
        if round(np.sum(self.year_area))!=round(np.sum(area)):
            print("Warning: Year area does not equal area")
            print("Year area: ", round(np.sum(self.year_area)))
            print("Area: ", round(np.sum(area)))
            print("Difference: ", round(np.sum(self.year_area)-np.sum(area)))
    
    def load_weather_data(self):
        #Load temperature and precipitation data
        df = pd.read_csv(self.input_files['temp_precip'])
        self.weather_dates = pd.to_datetime(df.iloc[:, 0], format="%Y/%m/%d").tolist()
        self.temps = df.iloc[:, 1].astype(float).to_numpy()
        self.precip = df.iloc[:, 2].apply(lambda x: float(x) if not np.isnan(float(x)) else 0).to_numpy()
        #MAKE NOTE OF THIS IN THESIS, SOME OF THE TEMPERATURE DATA IS INTERPOLATED BECAUSE ITS MISSING
        #If the temperature data is nan(missing data) then interpolate it
        nan_indices = np.isnan(self.temps)
        x = np.arange(len(self.temps))
        self.temps[nan_indices] = np.interp(x[nan_indices], x[~nan_indices], self.temps[~nan_indices])

    def load_verif_data(self):
        #Load the mass balance and ELA data
        #Skips the first 25 rows because they are before 1984
        df = pd.read_csv(self.input_files['mass_balance'], skiprows=25)
        self.annual_mb = df.iloc[:-1, 3].astype(float).to_numpy()
        self.summer_mb = df.iloc[:-1, 2].astype(float).to_numpy()
        self.winter_mb = df.iloc[:-1, 1].astype(float).to_numpy()
        self.ela_verif=df.iloc[5:, 4].astype(float).to_numpy()
        #Initialize the calculated mass balance arrays
        self.calculated_annual_mb=np.zeros(len(self.annual_mb))
        self.calculated_winter_mb=np.zeros(len(self.winter_mb))
        self.calculated_summer_mb=np.zeros(len(self.summer_mb))
        #Load the runoff data. 
        #I only have data from 1992-2007 and there are many gaps in that record so I keep track of what days I have data for 
        #and I keep track of that so I can only add the runoff data for those days to my calculated runoff arrays
        #Also the runoff data is summed by month so that there is less variability in the data
        df = pd.read_csv(self.input_files['runoff'])
        df['date'] = pd.to_datetime(df['date'])
        #Runoff training data variables
        training_start_date = "1997-08-15"
        training_end_date = "2001-07-15"
        training_data = df[(df["date"] >= training_start_date) & (df["date"] <= training_end_date)]
        self.measured_runoff_training = training_data.groupby(training_data["date"].dt.to_period("M"))["runoff"].sum().to_numpy()
        self.date_index_training = {date: idx for idx, date in enumerate(training_data['date'])}
        self.daily_runoff_training = np.zeros(len(self.date_index_training))
        #Runoff verification data variables
        verif_data = df[(df["date"] < training_start_date) | (df["date"] > training_end_date)]
        self.measured_runoff_verif = verif_data.groupby(verif_data["date"].dt.to_period("M"))["runoff"].sum().to_numpy()
        self.date_index_verif = {date: idx for idx, date in enumerate(verif_data['date'])}
        self.daily_runoff_verif = np.zeros(len(self.date_index_verif))
        #Runoff all data variables
        self.measured_runoff_all = df.groupby(df['date'].dt.to_period('M'))['runoff'].sum().to_numpy()
        self.date_index_all = {date: idx for idx, date in enumerate(df['date'])}
        self.month_index_all = sorted(set(date.strftime("%Y-%m") for date in self.date_index_all.keys()))
        self.daily_runoff_all_data = np.zeros(len(self.date_index_all))
        self.glacier_melt_runoff_data=np.zeros(len(self.date_index_all))
        #Load the thickness change and front variation data
        self.thickness_change_verif = pd.read_csv(self.input_files['thickness_change']).iloc[0:, 11].astype(float).to_numpy()
        self.front_variation_verif = pd.read_csv(self.input_files['front_variation']).iloc[0:, 9].astype(float).to_numpy()
        self.front_variation_calc = np.zeros(len(self.front_variation_verif))
        #Load the basin area data (used for the snow melt calculations)
        self.basin_areas, self.mean_snow_bin_elev=np.loadtxt(self.input_files['basin_area'], delimiter=',', usecols=(1, 2), unpack=True)
        self.snow_depth=np.zeros(len(self.basin_areas))
        self.snow_melt_amt=np.zeros(len(self.basin_areas))
        
    def calc_verif(self,timestep):
        #Round current date since all of the lists are based on day or year
        current_date_key = self.current_date.replace(hour=0, minute=0, second=0, microsecond=0)
        #These are the days before the start dates for the volume verif data so need to set the pre_volume to calculate Runoff
        # if current_date_key==datetime(2002,9,30) or  current_date_key==datetime(2003,6,8): self.prev_volume=self.ice_volume.copy()
        #If date is in the list of volume verification dates, calculate the total runoff
        if current_date_key in self.date_index_all: 
            self.daily_runoff_all_data[self.date_index_all[current_date_key]] += (self.glacial_melt+self.snow_melt_vol+self.rain_vol_per_step)
            self.glacier_melt_runoff_data[self.date_index_all[current_date_key]]+=self.glacial_melt
        if current_date_key in self.date_index_training: 
            self.daily_runoff_training[self.date_index_training[current_date_key]] += (self.glacial_melt+self.snow_melt_vol+self.rain_vol_per_step)
        if current_date_key in self.date_index_verif: 
            self.daily_runoff_verif[self.date_index_verif[current_date_key]] += (self.glacial_melt+self.snow_melt_vol+self.rain_vol_per_step)
        #Add the runoff (or glacial melt) data to the correct bin of the monthly runoff arrays.
        #These arrays have a bin for each month for each year from 1984-2024
        self.monthly_runoff_all[(self.current_date.year - 1984) * 12 + (self.current_date.month - 1)] += (self.glacial_melt+self.snow_melt_vol+self.rain_vol_per_step)
        self.monthly_glacier_melt[(self.current_date.year - 1984) * 12 + (self.current_date.month - 1)] += self.glacial_melt
        self.monthly_snow_depth[(self.current_date.year - 1984) * 12 + (self.current_date.month - 1)] += np.mean(self.snow_depth)
        self.monthly_precip_vol[(self.current_date.year - 1984) * 12 + (self.current_date.month - 1)] += self.snow_melt_vol+self.rain_vol_per_step
        self.month_count[(self.current_date.year - 1984) * 12 + (self.current_date.month - 1)]+=1
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
        #Calcualte ice thickness data for 1986 and 2021
        if self.current_date.year==1958: self.ice_1958=self.ice.copy()
        if self.current_date.year==1986: self.ice_1986=self.ice.copy()
        if self.current_date.year==2021: self.ice_2021=self.ice.copy()
        #Add mass balance to verification arrays
        if 1984<=self.current_date.year<2024:
            date=int(self.current_date.year-1984)
            self.calculated_annual_mb[date]+=np.mean(np.array(self.b*timestep)) if self.b.size>0 else 0
            self.calculated_winter_mb[date]+=np.mean(np.array(self.b[self.b>0]*timestep)) if self.b[self.b>0].size>0 else 0
            self.calculated_summer_mb[date]+=np.mean(np.array(self.b[self.b<0]*timestep)) if self.b[self.b<0].size>0 else 0
        #Make sure mass balance verification arrays are correct
        if np.any(self.calculated_winter_mb<0): print("ERROR IN WINTER MB")
        if np.any(self.calculated_summer_mb>0): print("ERROR IN SUMMER MB")
    
    #This function performs all the calculations for the snow model. This models the snow off the glacier
    def snow_model(self, index,timestep):
        #Get the temperature for the day either using the monthly lapse rate list or a fixed lapse rate
        if type(self.temp_lapse_rate) is int: snow_temps=self.temps[index]+self.temp_lapse_rate*(self.mean_snow_bin_elev-272)
        else: snow_temps=self.temps[index]+self.temp_lapse_rate[self.current_date.month-1]*(self.mean_snow_bin_elev-272)
        #Calculate snow accumulation using the precipitation for the day multiplied by the precipitation conversion factor
        self.snow_depth[snow_temps<=0]+=(self.precip[index]/1000*self.precip_conv_factor)*timestep
        #Calculate the amount of snow to melt using the snow melt factor and the temperature for the day
        self.snow_melt_amt.fill(0)
        self.snow_melt_amt[snow_temps>0]=self.snow_melt_factor*snow_temps[snow_temps>0]*timestep
        #Make sure snow melt doesn't exceed snow depth
        self.snow_melt_amt=-np.minimum(np.abs(self.snow_melt_amt),self.snow_depth)
        if np.any(self.snow_melt_amt>0): print("SNOW MELT POSITIVE")
        #Melt the snow
        self.snow_depth+=self.snow_melt_amt
        #Calculate the volume of snow melt using the snow melt amount and the basin areas minus the glacier areas
        self.snow_melt_vol=np.sum(np.abs(self.snow_melt_amt[:-1])*np.maximum((self.basin_areas[:-1]-self.areas[(self.current_date.year-1950),:]),0))+(np.abs(self.snow_melt_amt[-1])*self.basin_areas[-1])
        #Also calculate the amount of rain that fell. This is calculated for the entire basin (including the glacier)
        self.rain_vol_per_step=np.sum((self.precip[index]/1000*self.precip_conv_factor)*timestep*(self.basin_areas[snow_temps>0]))
        #Now avalanche the snow. 
        #This avalances a user defined percentage of the snow off the top 4 bins (the number of bins is set by the length of avalanche_dates) and distributes it evenly over the bottom 6 bins.
        matching_indices = [i for i, d in enumerate(self.avalanche_dates) if d == self.current_date.replace(hour=0, minute=0, second=0, microsecond=0)]
        for elem in matching_indices:
            snow_to_move=self.snow_depth[-elem]*self.avalanche_percent
            self.snow_depth[-elem]-=snow_to_move
            for i in range(7):
                self.snow_depth[i]+=snow_to_move/6
        self.snow_depth_list.append(np.mean(self.snow_depth))
        if self.snow_melt_vol<0: print("NEGATIVE SNOW MELT")

    def update_b(self, timestep):
        if self.current_date>=datetime(1984,1,2):
            #Calculate which index to get weather data from
            if self.current_date<datetime(2024,10,1): index=self.weather_dates.index(pd.Timestamp(self.current_date.replace(hour=0, minute=0, second=0, microsecond=0)))
            else: index=self.weather_dates.index(pd.Timestamp(datetime(2024, 9, 30)))
            #Calculates temperatures for every self.x value and varies temp with elevation
            if type(self.temp_lapse_rate) is int: x_temps=self.temps[index]+self.temp_lapse_rate*(self.ice+self.topo-272) #weather station elevation is 272m
            else: x_temps=self.temps[index]+self.temp_lapse_rate[self.current_date.month-1]*(self.ice+self.topo-272)
            mb=np.zeros_like(x_temps) #initialize mass balance
            #Melts ice for temps greater than 0
            #Above the ELA the snow melt factor is used and below the ELA a linear gradient is defined 
            #starting with the snow melt factor at the ELA and gradually shifts to the ice melt factor at the base of the glacier
            melt_arr=self.snow_melt_factor+((self.curr_ela-(self.ice+self.topo))/(self.curr_ela-np.nanmin(self.topo+self.ice)))*(self.ice_melt_factor-self.snow_melt_factor)
            melt_arr[melt_arr<self.ice_melt_factor]=self.ice_melt_factor
            mb[np.where((x_temps>0)&((self.ice+self.topo)>=self.curr_ela))[0]]=self.snow_melt_factor*x_temps[np.where((x_temps>0)&((self.ice+self.topo)>self.curr_ela))[0]]
            # mb[np.where((x_temps>0)&((self.ice+self.topo)<self.curr_ela))[0]]=melt_arr[np.where((x_temps>0)&((self.ice+self.topo)<self.curr_ela))[0]]*x_temps[np.where((x_temps>0)&((self.ice+self.topo)<self.curr_ela))[0]]
            mb[np.where((x_temps>0)&((self.ice+self.topo)<self.curr_ela))[0]]=melt_arr[np.where((x_temps>0)&((self.ice+self.topo)<self.curr_ela))[0]]*x_temps[np.where((x_temps>0)&((self.ice+self.topo)<self.curr_ela))[0]]
            #Calculates snow melt and accumulation
            self.snow_model(index,timestep)
            #The accumulation factor (precipitation to ice) is also defined by a linear gradient.
            #This gradient changes by year and starts as accumfactor_lower in 1984 and gradually increases to accumfactor_upper in 2024
            accumfactor = self.accumfactor_lower + ((self.current_date.year - 1984) / (2024 - 1984)) * (self.accumfactor_upper-self.accumfactor_lower)
            mb[x_temps<0]=self.precip[index]/1000*self.precip_conv_factor*accumfactor*timestep
            #Make sure the model can't melt more ice than the glacier has
            mb[np.where((mb<0)&(np.abs(mb)>self.ice))] = -self.ice[np.where((mb<0)&(np.abs(mb)>self.ice))]
            #Calculate how much the glacier melted using the mass balance and area data
            self.glacial_melt=np.sum(np.abs(mb[mb<0])*self.year_area[mb<0])*timestep
            self.precip_accum[self.current_date.year-1984]+=np.sum((self.precip[index]/1000)*self.precip_conv_factor*timestep)
            #Alert for weird mass balance values
            if np.any(mb[x_temps<0]<0): print("NEGATIVE MB", mb[x_temps<0])
            if np.any(np.isnan(mb[x_temps < 0])): print("NAN MB", self.precip[index])
            if np.any(np.abs(mb)>100): print("MB TOO BIG",mb)
            return mb
        else: 
            #Used to spin up the glacier. Glacier starts retreating in 1900
            # print(self.gamma)
            # print(self.curr_ela)
            # print(self.ice)
            # print(self.topo)
            return ((self.topo+self.ice-self.curr_ela)*self.gamma)/365.25 #meters per day
    
    def calc_q(self):
        #Calculate ice slope
        self.ice_slope[:-1] = -(np.diff(self.ice+self.topo) / self.dx)
        #Check for weird ice and ice_slope values
        if np.any(np.isnan(self.ice_slope)): #and not self.quiet:
            print('NaN detected in ice_slope:', self.ice_slope)
            print("MASS BALANCE: ", self.b)
            print("TIME: ", self.time)
            plt.plot(self.x, self.ice+self.topo)
            plt.plot(self.timestep_list)
            return
        if np.any(np.isnan(self.ice)): #and not self.quiet:
            print('NaN detected in ice:', self.ice)
            print("TIME: ", self.time)
            return
        if np.any(np.isinf(self.ice_slope)): # and not self.quiet:
            print('Infinity detected in ice_slope:', self.ice_slope)
            print('Ice: ',self.ice)
            print("Q: ",self.q)
            print("TIME: ", self.time)
            plt.plot(self.timestep_list)
            return
        #Calculate ice flux
        #self.q[1:]=2e-17* ((self.p*self.g*np.sin(np.arctan(self.ice_slope)))**3)*((self.ice**5)/5) #per year
        self.q[1:]=5.87e-19* ((self.p*self.g*np.sin(np.arctan(self.ice_slope)))**3)*((self.ice**5)/5) #per day
        #Check for weird q values
        if np.any(np.isnan(self.q)): #and not self.quiet:
            print('NaN detected in q:', self.q)
            print(self.ice_slope)
            print(self.ice)
            print("TIME: ", self.time)
            #plt.plot(self.x, self.ice+self.topo)
            plt.plot(self.timestep_list)
            return
        #Print model status every self.save years if quiet is false
        if (self.prev_display==0.0 or ((self.time>=(self.prev_display+(self.save/365.25))*365.25) and self.time<(self.time*365.25))) and not self.quiet:
            print("TIME: ", math.floor(self.time/365.25))
            print("ELA: ", self.curr_ela)
            print("ICE: ",self.ice)
            print("SLOPE: ",self.ice_slope)
            print("MASS BALANCE: ", self.b)
            print("Q: ",self.q)
            print("DQDX: ",np.diff(self.q)/self.dx)
            print("SUM Q: ", np.sum(self.q))
            print("SUM DQDX: ",np.sum(np.diff(self.q)/self.dx))
            print()
            self.prev_display=self.time/365.25
        #Calculate the change in ice thickness and return it
        return (self.b-(np.diff(self.q)/self.dx)) 

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
            self.ela_line_list.append((self.x, curr_ela_plt, "r", "dashed", "ELA"))
            self.title_list.append('Year = ' + str(self.current_date.year))
            return
        iter_time=0.0 #days
        timestep=0.0 #days
        while(iter_time<self.save):
            #End while loop if the difference between the save and current iter_time is less than 1, this is done to make sure it doesn't run over
            if (self.save-iter_time)<1: break
            #Calculate glacier velocity for timestep calculation
            u = (5.87e-19)*((self.p*self.g*np.sin(np.arctan(self.ice_slope)))**3)*(((self.ice)**4)/5)
            #Calculate the timestep based on glacier velocity, the 0.2 multiplier can be changed, increase it to increase the timestep
            #Timestep is constrained between 0.0001 and 1 and rouned to the 5th decimal place, its pretty much always 1
            timestep = round(np.clip(((self.dx / np.nanmax(u)) * 0.2), 0.0001, 1),4) if np.any(u > 0) else 1
            self.timestep_list.append(timestep)
            self.b=self.update_b(timestep)
            dqdx=self.calc_q()
            #Set new ice thickness
            self.ice = np.maximum((self.ice + dqdx * timestep), 0.0)
            #Update time variables
            iter_time+=float(timestep)
            self.time+=timestep
            self.current_date+=timedelta(days=float(timestep))
            #Calculate verification data
            self.calc_verif(timestep)
            #Update mass balance for current year
            self.year_mb+=self.b*timestep
            #Set new b_min and b_max
            self.b_max = max(np.max(self.b*365.25),self.b_max)
            self.b_min = min(np.min(self.b*365.25),self.b_min)
            #Update area if area data is available for the current date, updates on the 1st of the year
            # if self.current_date>=datetime(1950,1,1) and self.current_date.day==1 and self.current_date.month==1: self.update_areas()
            #This ela_1900 is used for the model spinup to set the glacier into retreat before 1984
            if self.current_date.year==1900: self.curr_ela=self.ela_1900
            #Calculate the avalanche dates for the upcoming year on the 1st of the year
            if self.current_date.month==1 and 1<=self.current_date.day<2:
                for i in range(len(self.avalanche_dates)):
                    self.avalanche_dates[i]=((self.current_date+timedelta(days=2)) + timedelta(days=random.randint(0, (datetime(self.current_date.year,3,31)-(self.current_date+timedelta(days=1))).days))).replace(hour=0, minute=0, second=0, microsecond=0)
        #If this if-statement runs somehow the model ran over which is a problem. The first if statement in the while loop should prevent this
        if(iter_time>self.save):
            print("BROKEN ABORT ABORT", iter_time,self.save)
            return
        #This if statement is also to make sure it doesn't run over. Since self.save should never be a fraction of a year then the current run should always end on 12/31
        if(self.current_date.day!=1):
            #Calculate the timestep so that it doesn't go over
            if(self.current_date.day==31): timestep=float((datetime((self.current_date.year+1),1,1)-self.current_date).total_seconds()/ (24 * 3600))
            else: timestep=float(self.save-iter_time)
            #Calculate last ice_flux for this run
            dqdx=self.calc_q()
            self.timestep_list.append(timestep)
            self.b=self.update_b(timestep)
            #Set new ice thickness
            self.ice = np.maximum((self.ice + dqdx * timestep), 0.0)
            #Update time variables
            iter_time+=float(timestep)
            self.time+=timestep
            self.current_date+=timedelta(days=float(timestep))
            #Calculate verification data
            self.calc_verif(timestep)
            #Update mass balance for current year
            self.year_mb+=self.b*timestep
            #Set new b_min and b_max
            self.b_max = max(np.max(self.b*365.25),self.b_max)
            self.b_min = min(np.min(self.b*365.25),self.b_min)
            #Update area if area data is available for the current date, updates on the 1st of the year
            # if self.current_date>=datetime(1950,1,1) and self.current_date.day==1 and self.current_date.month==1: self.update_areas()
            if self.current_date.year==1900: self.curr_ela=self.ela_1900
        #If the weather data is being used (starts in 1984) for the mass balance then calculate ela
        if self.current_date.year>=1984:
            #If all of the mass balance is positive the the ELA is at the bottom of the glacier
            if np.all(self.year_mb>0): curr_ela=float(self.topo[-1])
            #If all of the mass balance is negative then the ELA is at the top of the glacier
            elif np.all(self.year_mb<0): curr_ela=float(self.topo[0])
            #Otherwise find the ELA using the mass balance value that is on the border of positive and negative
            #If for some bizarre reason there is multiple transitions from positive to negative this will chose the one highest on the glacier
            #But that shouldn't be possible
            else: curr_ela = float(self.topo[np.where(np.sign(self.year_mb[:-1]) != np.sign(self.year_mb[1:]))[0][0]])
            #Add to ela list after 1989 because that's when ela verification data starts
            if self.current_date.year>1989: self.ela_list.append(curr_ela)
            #Reset year_mb to 0 for the next year
            # self.year_mb.fill(0)
            #Calculate the avalanche dates for the upcoming year on the 1st of the year
            if self.current_date.month==1 and 0<=self.current_date.day<=1:
                delta = datetime(self.current_date.year,3,31)-(self.current_date+timedelta(days=1))
                for i in range(len(self.avalanche_dates)):
                    self.avalanche_dates[i]=((self.current_date+timedelta(days=2)) + timedelta(days=random.randint(0, delta.days))).replace(hour=0, minute=0, second=0, microsecond=0)
        #Otherwise set ela to the start ela
        else: curr_ela=self.curr_ela
        #Make sure to reset year_mb to 0 on the first of the year
        if self.current_date.month==1 and 1<=self.current_date.day<=2: self.year_mb.fill(0)
        curr_ela_plt=[curr_ela]*self.num_cells
        #If the current date is the end date then report final values
        if(self.current_date.year==1484+self.time): self.report_final_values(u)
        #Update width if width data is available for the current date, updates on the 1st of the year
        if self.current_date>=datetime(1950,1,1): self.update_areas()
        #Add data to plotting lists
        self.ice_line_list.append((self.x, self.ice + self.topo, "c", "Ice"))
        # self.snow_line_list.append((self.x, self.snow_depth + self.topo+self.ice, "c", "Snow"))
        self.ela_line_list.append((self.x, curr_ela_plt, "r", "dashed", "ELA"))
        self.title_list.append('Year = ' + str(self.current_date.year))

    #Used to plot the model data after the model has run
    def func(self,i,ax):
        ax.clear() # Remove this to make lines stack up
        # ax.set_ylim(min(self.topo) - 100, max(self.topo) + 100)
        ax.set_xlim(0, float(self.valley_length))
        ax.set_ylabel("Height (m)")
        ax.set_xlabel("Distance (m)")
        ax.set_aspect('equal', adjustable='datalim')
        ax.plot(self.x, self.topo, color="b", label="Topography")
        ax.set_title(self.title_list[i])
        self.ela_line=ax.plot(self.ela_line_list[i][0],self.ela_line_list[i][1], color=self.ela_line_list[i][2], linestyle=self.ela_line_list[i][3], label=self.ela_line_list[i][4])
        self.ice_line=ax.plot(self.ice_line_list[i][0],self.ice_line_list[i][1], color=self.ice_line_list[i][2], label=self.ice_line_list[i][3])
        # self.snow_line=ax.plot(self.snow_line_list[i][0],self.snow_line_list[i][1], color=self.snow_line_list[i][2], label=self.snow_line_list[i][3])
        ax.legend()
        return self.ice_line_list,self.ela_line_list         