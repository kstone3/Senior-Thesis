from datetime import datetime, timedelta
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from math import inf
import matplotlib
matplotlib.use("TkAgg")

class glacierSim():
    def __init__(self, ela=2000,valley_length=3668, time=500,save=10,gamma=0.01,quiet=True, ice_meltfactor=0.00006, snow_meltfactor=0.00001, accumfactor=0.00006, initial_ice=None, start_time=0):
        #print("OBJ INIT")
        self.valley_length = valley_length
        self.run_time=start_time*365.25 #DAYS
        self.prev_display=0
        self.start_ela = ela #in m
        self.curr_ela=ela
        self.num_cells = 50 #set number of cells
        self.dx = self.valley_length/(self.num_cells-1) #cell width in m
        self.time = time #simulation time in years
        self.current_date=datetime(984,1,1)+timedelta(days=start_time*365.25)
        self.save = save #timestep interval in years
        self.frames = ((int)((self.time-start_time)/self.save))+1 #number of frames the animation will run for
        if start_time==0: self.ice = np.zeros(self.num_cells) #initialize ice
        else: self.ice = np.array(initial_ice)
        self.q = np.zeros(self.num_cells+1) #initialize ice flux, need to have num_cells+1 to offset it from ice thickness for calculations
        self.x = np.linspace(0.5 * self.dx, self.valley_length - 0.5 * self.dx,self.num_cells) #used for plotting topography and ice thickness
        self.topo =[] #initialize topography
        self.g = 9.81 #gravity constant in m/yr^2
        self.p = 917 #density of ice
        self.b_max = float(-inf)
        self.b_min = float(inf)
        self.gamma = gamma #for mass balance equation
        self.ice_volume=0
        self.volume_change=[]
        self.yearly_volume_change=np.zeros(41)
        self.volume_validation=np.zeros(180)
        self.volume_data=[]
        self.prev_volume=0
        self.initial_volume=0
        self.timestep_list=[] #days
        self.initial_run=False
        self.b=np.zeros(self.num_cells)
        self.glacier_extent=0
        self.snow_depth=np.zeros(self.num_cells)
        self.ice_slope = np.zeros(self.num_cells, dtype=np.longdouble) #initialize ice_slope
        self.quiet=quiet
        self.dates=[]
        self.temps=[]
        self.precip=[]
        self.annual_mb=[]
        self.winter_mb=[]
        self.summer_mb=[]
        self.bins=[]
        self.years=[]
        self.ela_list=[]
        self.areas=[]
        self.calculated_annual_mb=np.zeros(40, dtype=np.float64)
        self.calculated_winter_mb=np.zeros(40, dtype=np.float64)
        self.calculated_summer_mb=np.zeros(40, dtype=np.float64)
        self.annual_mb_arr=np.zeros(len(self.b))
        self.ice_melt_factor=ice_meltfactor
        self.snow_melt_factor=snow_meltfactor
        self.accum_factor=accumfactor
        self.widths=np.zeros(self.num_cells)
        self.widths_over_time=[]
        self.ice_thickness_over_time=[]
        self.date_index={date: idx for idx, date in enumerate(list(pd.date_range("2002-10-01", "2002-12-05").to_pydatetime()) + list(pd.date_range("2003-06-09", "2003-09-30").to_pydatetime()))}
    
    def init(self, ax,ela=6700,valley_length=3668, time=500,save=10,gamma=0.008, quiet=True, ice_meltfactor=0.00006, snow_meltfactor=0.00001,accumfactor=0.00006, initial_ice=None, start_time=0):
        #print("INIT")
        self.__init__(ela, valley_length, time, save, gamma, quiet, ice_meltfactor, snow_meltfactor, accumfactor, initial_ice, start_time)
        self.calc_topo()
        self.calc_widths()
        self.load_mb_data()
        self.load_verif_data()
        try: curr_ela=self.topo[np.where((self.b[:-1] >= 0) & (self.b[1:] < 0))[0]][0]
        except: curr_ela=self.topo[-1]
        self.ela_list.append(curr_ela)
        ax.clear()
        ax.set_ylim(min(self.topo) - 100, max(self.topo) + 100)
        ax.set_xlim(0, self.valley_length)
        ax.set_ylabel("Height (m)")
        ax.set_xlabel("Distance (m)")
        ax.set_aspect('equal', adjustable='datalim')
        ax.plot(self.x, self.topo, color="b", label="Topography")
        self.line = ax.plot(self.x, self.ice + self.topo, color="c", label="Ice Thickness")
        self.snow_line = ax.plot(self.x, self.snow_depth + self.topo+self.ice, color="g", label="Snow")
        self.ela_line = ax.axhline(y=float(curr_ela), color="r", linestyle="dashed", label="ELA")
        ax.legend()

    def haversine(self, lat1, lon1, lat2, lon2):
        a = math.sin(math.radians(lat2 - lat1) / 2) ** 2 + math.cos(math.radians(lat1)) * math.cos( math.radians(lat2)) * math.sin(math.radians(lon2 - lon1) / 2) ** 2
        return 6371000*(2 * math.atan2(math.sqrt(a), math.sqrt(1 - a)))
        
    def calc_topo(self):
        #print("CALC TOPO")
        latitudes = []
        longitudes = []
        topo = []
        df = pd.read_csv('Data/centerlineBed.csv',skiprows=1)
        latitudes = df.iloc[:, 2].astype(float).tolist()  # Latitude is the second column (index 2)
        longitudes = df.iloc[:, 1].astype(float).tolist()  # Longitude is the third column (index 1)
        topo = df.iloc[:, 0].astype(float).tolist()  # Elevation is the first column (index 0)
        cumulative_distances=[0.0]
        for i in range(1, len(latitudes)): cumulative_distances.append(cumulative_distances[-1] + self.haversine(latitudes[i - 1], longitudes[i - 1], latitudes[i], longitudes[i]) )
        self.x=np.linspace(cumulative_distances[0], cumulative_distances[-1], self.num_cells)
        self.topo=np.interp(np.linspace(cumulative_distances[0], cumulative_distances[-1], self.num_cells), cumulative_distances, topo)
        self.valley_length=max(np.max(self.x),100)
        self.dx = self.valley_length/(self.num_cells-1)
        self.default_b = True
        self.ice_slope[:-1] = abs((np.diff(self.topo)/ self.dx))
        #self.ice=np.exp(-0.5 * ((self.x - 700) / 300) ** 2)
        #scale=300/np.max(self.ice)
        #self.ice*=scale
        self.ice_volume=np.sum(self.ice*self.widths*self.dx)
        self.initial_volume=self.ice_volume
        self.prev_volume=self.ice_volume
        self.volume_change.append(self.initial_volume)
        
    def calc_widths(self):
        df = pd.read_csv('Data/Input_SouthCascade_Area_Altitude_Distribution.csv')
        self.bins = df.columns[1:].astype(float).to_numpy()
        self.years = df.iloc[:, 0].astype(float).tolist()
        self.areas = df.iloc[:, 1:].astype(float).values
        self.widths=np.array([self.areas[0][np.array([np.argmin(np.abs(self.bins[:, None]-np.round(self.topo)), axis=0)])]][0][0])/0.05*1000
        
    def update_widths(self):
        prev_width=self.widths
        self.widths=np.array([self.areas[self.current_date.year-1950][np.array([np.argmin(np.abs(self.bins[:, None]-np.round(self.topo)), axis=0)])]][0][0])/0.05*1000 # type: ignore
        if np.abs(self.widths-prev_width).any()>0:
            self.widths_over_time.append(self.widths.copy())
            self.ice_thickness_over_time.append(self.ice.copy())
        
    def load_mb_data(self):
        df = pd.read_csv('Data/Input_SouthCascade_Daily_Weather.csv')
        self.dates = pd.to_datetime(df.iloc[:, 0], format="%Y/%m/%d").tolist()
        self.temps = df.iloc[:, 1].astype(float).tolist()
        self.precip = df.iloc[:, 2].apply(lambda x: float(x) if not np.isnan(float(x)) else 0).to_numpy()

    def load_verif_data(self):
        df = pd.read_csv('Data/Output_SouthCascade_Glacier_Wide_solutions_calibrated.csv', skiprows=25)
        self.annual_mb = df.iloc[:, 3].astype(float).tolist()
        self.summer_mb = df.iloc[:, 2].astype(float).tolist()
        self.winter_mb = df.iloc[:, 1].astype(float).tolist()
        self.calculated_annual_mb=np.array([0] * len(self.annual_mb), dtype=np.float64)
        self.calculated_winter_mb=np.array([0] * len(self.winter_mb), dtype=np.float64)
        self.calculated_summer_mb=np.array([0] * len(self.summer_mb), dtype=np.float64)
        self.volume_data = pd.read_csv('Data/daily_average_runoff_with_dates.csv')['Volume'].to_numpy()
        self.volume_validation=np.zeros(len(self.volume_data))
    
    def snow_model(self, index, temps):
        x_temps=self.temps[index]+-0.004*(self.ice+self.topo-272)
        self.snow_depth[x_temps<0]+=((self.precip[index]*15)/1000) - (self.snow_depth[x_temps<0]*self.accum_factor)
        self.snow_depth[x_temps>0]+=self.snow_melt_factor*temps[x_temps>0]
    
    def update_b(self):
        if self.current_date>=datetime(1984,1,1):
            if self.current_date<datetime(2024,10,1): index=self.dates.index(pd.Timestamp(self.current_date.replace(hour=0, minute=0, second=0, microsecond=0)))
            else: index=self.dates.index(pd.Timestamp(datetime(2024, 9, 30)))
            x_temps=self.temps[index]+-0.004*(self.ice+self.topo-272)
            mb=np.zeros_like(x_temps)
            mb[x_temps>0]=self.ice_melt_factor*x_temps[x_temps>0]
            mb[(x_temps>0)&(self.snow_depth>0)]=0
            #mb[x_temps<0]=self.accum_factor*self.precip[index]/1000*self.dx
            self.snow_model(index,x_temps)
            mb[x_temps<0]=self.snow_depth[x_temps<0]*self.accum_factor
            # if np.any(mb[x_temps<0]<0): print("NEGATIVE MB", mb[x_temps<0])
            # if np.any(np.isnan(mb[x_temps < 0])): print("NAN MB", self.precip[index])
            # if np.any(mb[x_temps<0]>100): print("MB TOO BIG", mb[x_temps<0])
            #print(math.floor(self.run_time/365.25)-1000)
            #date=math.floor((self.run_time-1000*365.25)/365.25)
            date=int(self.current_date.year-1984)
            self.calculated_annual_mb[date]+=np.sum(np.array(mb))
            self.calculated_winter_mb[date]+=np.sum(np.array(mb[mb>0]))
            self.calculated_summer_mb[date]+=np.sum(np.array(mb[mb<0]))
            if np.any(self.calculated_winter_mb<0): print("ERROR IN WINTER MB")
            if np.any(self.calculated_summer_mb>0): print("ERROR IN WINTER MB")
            #print(self.calculated_annual_mb)
            return mb
        else: return ((self.topo+self.ice-self.curr_ela)*self.gamma)/365.25 #meters per day
        
    # def update_b(self):
    #     if self.run_time>(1000*365.25):
    #         target_date = datetime(1984, 1, 2) + timedelta(days=math.floor(self.run_time-(1000*365.25)))
    #         index=self.dates.index(target_date)
    #         x_temps=self.temps[index]+-0.004*(self.ice+self.topo-272)
    #         mb=np.zeros_like(x_temps)
    #         mb[x_temps>0]=self.melt_factor*x_temps[x_temps>0]
    #         mb[x_temps<0]=self.accum_factor*self.precip[index]/1000*self.dx
    #         # if np.any(mb[x_temps<0]<0): print("NEGATIVE MB", mb[x_temps<0])
    #         # if np.any(np.isnan(mb[x_temps < 0])): print("NAN MB", self.precip[index])
    #         # if np.any(mb[x_temps<0]>100): print("MB TOO BIG", mb[x_temps<0])
    #         #print(math.floor(self.run_time/365.25)-1000)
    #         self.calculated_annual_mb[math.floor((self.run_time-1000*365.25)/365.25)]+=np.sum(np.array(mb))
    #         self.calculated_winter_mb[math.floor((self.run_time-1000*365.25)/365.25)]+=np.sum(np.array(mb[mb>0]))
    #         self.calculated_summer_mb[math.floor((self.run_time-1000*365.25)/365.25)]+=np.sum(np.array(mb[mb<0]))
    #         if np.any(self.calculated_winter_mb<0): print("ERROR IN WINTER MB")
    #         if np.any(self.calculated_summer_mb>0): print("ERROR IN WINTER MB")
    #         #print(self.calculated_annual_mb)
    #         return mb
    #     else: return ((self.topo+self.ice-self.curr_ela)*self.gamma)/365.25 #meters per day
    
    def calc_q(self):
        self.ice_slope[:-1] = -(np.diff(self.ice+self.topo) / self.dx) #calculate ice slope
        if np.any(np.isnan(self.ice_slope)): #and not self.quiet:
            print('NaN detected in ice_slope:', self.ice_slope)
            print("MASS BALANCE: ", self.b)
            print("TIME: ", self.run_time)
            #plt.plot(self.x, self.ice+self.topo)
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
        #self.q[1:] = (0.2 *(2e-17*2) * (self.p * self.g)**2) * np.sin(np.arctan(self.ice_slope))**2 * (self.ice**5)/5
        #self.q[1:]=2e-17* ((self.p*self.g*np.sin(np.arctan(self.ice_slope)))**3)*((self.ice**5)/5) #YEARS
        self.q[1:]=5.87e-19* ((self.p*self.g*np.sin(np.arctan(self.ice_slope)))**3)*((self.ice**5)/5) #DAYS
        if np.any(np.isnan(self.q)): #and not self.quiet:
            print('NaN detected in q:', self.q)
            print(self.ice_slope)
            print(self.ice)
            print("TIME: ", self.run_time)
            #plt.plot(self.x, self.ice+self.topo)
            plt.plot(self.timestep_list)
            return
        if (self.prev_display==0.0 or ((self.run_time>=(self.prev_display+self.save)*365.25) and self.run_time<(self.time*365.25))) and not self.quiet:
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

        return self.ice_slope,(self.b-(np.diff(self.q)/self.dx)) #calculate the change in ice thickness and return it

    def report_final_values(self,u):
        ice_extent = self.x[self.ice > 1]
        if ice_extent.size > 0:
            self.glacier_extent=np.max(ice_extent)
            print("Final glacier length: " + str(np.max(ice_extent)) + 'm')
        else: print("Final glacier length: 0m (no ice extent)")
        print("Final max ice thickness: " + str(np.max(self.ice)) + 'm')
        print('Final max velocity: ' +str(np.max(u)) + "m/yr")
        if(self.default_b):
            print('B min: ' + str(self.b_min))
            print('B max: ' + str(self.b_max))

    def run_model(self,i):
        if i==0 and not self.initial_run:
            self.initial_run=True
            try: curr_ela=float(self.topo[np.where((self.annual_mb_arr[:-1] >= 0) & (self.annual_mb_arr[1:] < 0))[0][0]])
            except: curr_ela=float(self.topo[-1])
            self.ela_list.append(curr_ela)
            ax.clear()
            ax.set_ylim(min(self.topo) - 100, max(self.topo) + 100)
            ax.set_xlim(0, float(self.valley_length))
            ax.set_ylabel("Height (m)")
            ax.set_xlabel("Distance (m)")
            ax.set_aspect('equal', adjustable='datalim')
            ax.plot(self.x, self.topo, color="b", label="Topography")
            #ax.set_title('Time = ' + str(round((self.run_time/365.25),(str(self.save)[::-1].find('.')+1))) + ' years')
            ax.set_title('Time = ' + str(self.current_date.year) + ' years')
            self.ela_line = ax.axhline(y=curr_ela, color="r", linestyle="dashed", label="ELA")
            self.line = ax.plot(self.x, self.ice + self.topo, color="c", label="Ice")
            self.snow_line = ax.plot(self.x, self.snow_depth + self.topo+self.ice, color="g", label="Snow")
            ax.legend()
        elif i>0.0:
            iter_time=0.0 #years
            timestep=0.0 #days
            while(iter_time<self.save):
                if (self.save-iter_time)<(1/365.25): break
                ice_slope,dqdx=self.calc_q()
                u = (5.87e-19)*((self.p*self.g*np.sin(np.arctan(ice_slope)))**3)*(((self.ice)**4)/5)
                timestep = round(np.clip(((self.dx / np.max(u)) * 0.2), 0.0001, 1),5) if np.any(u > 0) else 1
                self.timestep_list.append(timestep)
                self.ice = np.maximum((self.ice + dqdx * timestep), 0.0)
                self.ice_volume=np.sum(self.ice*self.dx*self.widths)
                #self.volume_change.append(abs(self.ice_volume/1e9-self.initial_volume/1e9))
                current_date_key = self.current_date.replace(hour=0, minute=0, second=0, microsecond=0)
                if current_date_key==datetime(2002,9,30) or  current_date_key==datetime(2003,6,8): self.prev_volume=self.ice_volume
                if current_date_key in self.date_index:
                    self.volume_validation[self.date_index[current_date_key]] += (self.ice_volume-self.prev_volume)
                    self.prev_volume=self.ice_volume
                iter_time+=timestep/365.25
                self.run_time+=timestep
                self.current_date+=timedelta(days=float(timestep))
                self.b=self.update_b()
                if self.current_date.year>=1984: self.yearly_volume_change[self.current_date.year-1984]=self.ice_volume/1e9
                if math.floor(self.run_time%365.25)==0: self.annual_mb_arr.fill(0)
                self.annual_mb_arr+=self.b
                self.b_max = max(np.max(self.b*365.25),self.b_max)
                self.b_min = min(np.min(self.b*365.25),self.b_min)
                if self.current_date>=datetime(1950,1,1): self.update_widths()
            if(iter_time>self.save):
                print("BROKEN ABORT ABORT")
                return
            ice_slope,dqdx=self.calc_q()
            timestep=(self.save-iter_time)*365.25 
            self.timestep_list.append(timestep)
            self.ice = np.maximum((self.ice + dqdx * timestep), 0.0)
            self.ice_volume=np.sum(self.ice*self.dx*self.widths)
            #self.volume_change.append(abs(self.ice_volume/1e9-self.initial_volume/1e9))
            self.volume_change.append(self.ice_volume/1e9)
            self.run_time+=timestep
            self.current_date+=timedelta(days=float(timestep))
            current_date_key = self.current_date.replace(hour=0, minute=0, second=0, microsecond=0)
            if current_date_key==datetime(2002,9,30) or  current_date_key==datetime(2003,6,8): self.prev_volume=self.ice_volume
            if current_date_key in self.date_index:
                self.volume_validation[self.date_index[current_date_key]] += (self.ice_volume-self.prev_volume)
                self.prev_volume=self.ice_volume
            self.b=self.update_b()
            if math.floor(self.run_time%365.25)==0: self.annual_mb_arr.fill(0)
            self.annual_mb_arr+=self.b
            try: curr_ela=float(self.topo[np.where((self.annual_mb_arr[:-1] >= 0) & (self.annual_mb_arr[1:] < 0))[0][0]])
            except: curr_ela=float(self.topo[-1])
            self.ela_list.append(curr_ela)
            ax.clear()
            ax.set_ylim(min(self.topo) - 100, max(self.topo) + 100)
            ax.set_xlim(0, float(self.valley_length))
            ax.set_ylabel("Height (m)")
            ax.set_xlabel("Distance (m)")
            ax.set_aspect('equal', adjustable='datalim')
            ax.plot(self.x, self.topo, color="b", label="Topography")
            #ax.set_title('Time = ' + str(round((self.run_time/365.25),(str(self.save)[::-1].find('.')+1))) + ' years')
            ax.set_title('Time = ' + str(self.current_date.year) + ' years')
            self.ela_line = ax.axhline(y=curr_ela, color="r", linestyle="dashed", label="ELA")
            self.line = ax.plot(self.x, self.ice + self.topo, color="c", label="Ice")
            self.snow_line = ax.plot(self.x, self.snow_depth + self.topo+self.ice, color="g", label="Snow")
            ax.legend()
            if(round(self.run_time,2)==self.time*365.25): self.report_final_values(u)
        return self.line, self.ela_line, self.snow_line

fig, ax = plt.subplots() #initialize plotting variables
_ = plt.close(fig) #used to prevent an empty plot from displaying