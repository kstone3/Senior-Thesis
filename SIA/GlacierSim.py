
import os
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython.core.display import HTML
import numpy as np
import math
from math import inf
import csv
import matplotlib
matplotlib.use("TkAgg")

class glacierSim():
    def __init__(self, ela=2000,valley_length=3668.0, time=500,save=10,gamma=0.01,quiet=True):
        #print("OBJ INIT")
        self.valley_length = valley_length
        self.run_time=0.0
        self.prev_display=0
        self.start_ela = ela #in m
        self.curr_ela=ela
        self.num_cells = 50 #set number of cells
        self.dx = self.valley_length/(self.num_cells-1) #cell width in m
        self.time = time #simulation time in years
        self.save = save #timestep interval in years
        self.frames = ((int)(self.time/self.save))+1 #number of frames the animation will run for
        self.ice = np.zeros(self.num_cells) #initialize ice
        self.q = np.zeros(self.num_cells+1) #initialize ice flux, need to have num_cells+1 to offset it from ice thickness for calculations
        self.x = np.linspace(0.5 * self.dx, self.valley_length - 0.5 * self.dx,self.num_cells) #used for plotting topography and ice thickness
        self.topo = [] #initialize topography
        self.g = 9.81 #gravity constant in m/yr^2
        self.p = 917 #density of ice
        self.b_max = float(-inf)
        self.b_min = float(inf)
        self.gamma = gamma #for mass balance equation
        self.ice_volume=0
        self.volume_change=[]
        self.initial_volume=0
        self.timestep_list=[]
        self.initial_run=False
        self.b=np.zeros(self.num_cells)
        self.glacier_extent=0
        self.ice_slope = np.zeros(self.num_cells, dtype=np.longdouble) #initialize ice_slope
        self.quiet=quiet
        self.glacier_extent=0
    
    def init(self, ax,ela=6700,valley_length=3668, time=500,save=10,gamma=0.008):
        #print("INIT")
        self.__init__(ela, valley_length, time, save, gamma)
        self.calc_topo()
        ax.clear()
        ax.set_ylim(min(self.topo) - 100, max(self.topo) + 100)
        ax.set_xlim(0, self.valley_length)
        ax.set_ylabel("Height (m)")
        ax.set_xlabel("Distance (m)")
        ax.set_aspect('equal', adjustable='datalim')
        ax.plot(self.x, self.topo, color="b", label="Topography")
        self.line = ax.plot(self.x, self.ice + self.topo, color="c", label="Ice Thickness")
        self.ela_line = ax.axhline(y=self.start_ela, color="r", linestyle="dashed", label="ELA")
        ax.legend()


    def haversine(self, lat1, lon1, lat2, lon2):
        a = math.sin(math.radians(lat2 - lat1) / 2) ** 2 + math.cos(math.radians(lat1)) * math.cos( math.radians(lat2)) * math.sin(math.radians(lon2 - lon1) / 2) ** 2
        return 6371000*(2 * math.atan2(math.sqrt(a), math.sqrt(1 - a)))
        
    def calc_topo(self):
        #print("CALC TOPO")
        latitudes = []
        longitudes = []
        topo = []
        with open('centerlineBed.csv', 'r', newline='') as file:
            reader = csv.reader(file)
            next(reader)  # Skip the header line
            for row in reader:
                latitudes.append(float(row[2]))    # Latitude is the second column
                longitudes.append(float(row[1]))   # Longitude is the third column
                topo.append(float(row[0]))
        cumulative_distances=[0.0]
        for i in range(1, len(latitudes)): cumulative_distances.append(cumulative_distances[-1] + self.haversine(latitudes[i - 1], longitudes[i - 1], latitudes[i], longitudes[i]) )
        self.x=np.linspace(cumulative_distances[0], cumulative_distances[-1], self.num_cells)
        self.topo=np.interp(np.linspace(cumulative_distances[0], cumulative_distances[-1], self.num_cells), cumulative_distances, topo)
        self.valley_length=max(np.max(self.x),100.0)
        self.dx = self.valley_length/(self.num_cells-1)
        self.default_b = True
        self.ice_slope[:-1] = abs((np.diff(self.topo)/ self.dx))
        #self.ice=np.exp(-0.5 * ((self.x - 700) / 300) ** 2)
        #scale=300/np.max(self.ice)
        #self.ice*=scale
        self.ice_volume=self.dx*np.sum(self.ice)
        self.initial_volume=self.ice_volume
        
    def update_b(self): return (self.topo+self.ice-self.curr_ela)*self.gamma
    
    def calc_q(self):
        self.ice_slope[:-1] = -(np.diff(self.ice+self.topo) / self.dx) #calculate ice slope
        if np.any(np.isnan(self.ice_slope)) and not self.quiet:
            print('NaN detected in ice_slope:', self.ice_slope)
            print("TIME: ", self.run_time)
            #plt.plot(self.x, self.ice+self.topo)
            plt.plot(self.timestep_list)
            return
        if np.any(np.isnan(self.ice)) and not self.quiet:
            print('NaN detected in ice:', self.ice)
            print("TIME: ", self.run_time)
            return
        if np.any(np.isinf(self.ice_slope)):
            print('Infinity detected in ice_slope:', self.ice_slope)
            print('Ice: ',self.ice)
            print("Q: ",self.q)
            print("TIME: ", self.run_time)
            plt.plot(self.timestep_list)
            return
        #self.q[1:] = (0.2 *(2e-17*2) * (self.p * self.g)**2) * np.sin(np.arctan(self.ice_slope))**2 * (self.ice**5)/5
        self.q[1:]=2e-17* ((self.p*self.g*np.sin(np.arctan(self.ice_slope)))**3)*((self.ice**5)/5)
        if np.any(np.isnan(self.q)) and not self.quiet:
            print('NaN detected in q:', self.q)
            print(self.ice_slope)
            print(self.ice)
            print("TIME: ", self.run_time)
            #plt.plot(self.x, self.ice+self.topo)
            plt.plot(self.timestep_list)
            return
        if (self.prev_display==0.0 or self.run_time>=self.prev_display+self.save) and self.run_time<self.time and not self.quiet:
            print("TIME: ", self.run_time)
            print("ELA: ", self.curr_ela)
            print("ICE: ",self.ice)
            print("ICE VOLUME: ",self.ice_volume/1e9, " km^3")
            print("Difference from inital volume: ", self.ice_volume/1e9-self.initial_volume/1e9, " km^3")
            self.volume_change.append(abs(self.ice_volume/1e9-self.initial_volume/1e9))
            print("SLOPE: ",self.ice_slope)
            print("MASS BALANCE: ", self.b)
            print("Q: ",self.q)
            print("DQDX: ",np.diff(self.q)/self.dx)
            print("SUM Q: ", np.sum(self.q))
            print("SUM DQDX: ",np.sum(np.diff(self.q)/self.dx))
            print()
            self.prev_display=self.run_time

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
            ax.clear()
            ax.set_ylim(min(self.topo) - 100, max(self.topo) + 100)
            ax.set_xlim(0.0, float(self.valley_length))
            ax.set_ylabel("Height (m)")
            ax.set_xlabel("Distance (m)")
            ax.set_aspect('equal', adjustable='datalim')
            ax.plot(self.x, self.topo, color="b", label="Topography")
            ax.set_title('Time = ' + str(round(self.run_time,(str(self.save)[::-1].find('.')+1))) + ' years')
            self.ela_line = ax.axhline(y=self.start_ela, color="r", linestyle="dashed", label="ELA")
            self.line = ax.plot(self.x, self.ice + self.topo, color="c", label="Ice")
            ax.legend()
        elif i>0.0:
            iter_time=0.0
            timestep=0.0
            while(iter_time<(self.save-timestep)):
                ice_slope,dqdx=self.calc_q()
                u = (2e-17)*((self.p*self.g*np.sin(np.arctan(ice_slope)))**3)*(((self.ice)**4)/5)
                timestep = np.clip(((self.dx / np.max(u)) * 0.2), 0.0001, 0.01) if np.any(u > 0) else 0.01
                self.timestep_list.append(timestep)
                self.ice = np.maximum((self.ice + dqdx * timestep), 0.0)
                self.ice_volume=self.dx*np.sum(self.ice)
                iter_time+=timestep
                self.run_time+=timestep
                self.b=self.update_b()
                self.b_max = max(np.max(self.b),self.b_max)
                self.b_min = min(np.min(self.b),self.b_min)
                
            ice_slope,dqdx=self.calc_q()
            timestep=abs(iter_time-self.save)
            self.timestep_list.append(timestep)
            self.ice = np.maximum((self.ice + dqdx * timestep), 0.0)
            self.ice_volume=self.dx*np.sum(self.ice)
            self.run_time+=timestep
            self.b=self.update_b()
            
            ax.clear()
            ax.set_ylim(min(self.topo) - 100, max(self.topo) + 100)
            ax.set_xlim(0, float(self.valley_length))
            ax.set_ylabel("Height (m)")
            ax.set_xlabel("Distance (m)")
            ax.set_aspect('equal', adjustable='datalim')
            ax.plot(self.x, self.topo, color="b", label="Topography")
            ax.set_title('Time = ' + str(round(self.run_time,(str(self.save)[::-1].find('.')+1))) + ' years')
            self.ela_line = ax.axhline(y=self.start_ela, color="r", linestyle="dashed", label="ELA")
            self.line = ax.plot(self.x, self.ice + self.topo, color="c", label="Ice")
            ax.legend()
            if(round(self.run_time,2)==self.time): self.report_final_values(u)
        return (self.line)
    
fig, ax = plt.subplots() #initialize plotting variables
_ = plt.close(fig) #used to prevent an empty plot from displaying
ela=1880
time=1000
save=100
gamma=0.00395
model = glacierSim(ela=ela, time=time, save=save,gamma=gamma)
# plt.plot(glac.x, glac.topo)
# plt.plot(glac.x,glac.ice+glac.topo, color='b')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.show()
#print(model.dx)
anim = FuncAnimation(fig, model.run_model, model.frames, init_func=model.init(ax,ela=ela, time=time, save=save,gamma=gamma), blit=False, repeat=False)
#vid = HTML(anim.to_jshtml())
try:
    os.remove('glacier_simulation.gif')
except:
    print()
anim.save('glacier_simulation.gif', writer='pillow') 
print("ICE: ",model.ice)
print("SLOPE: ",model.ice_slope)
print("MASS BALANCE: ",model.b)
print("DONE")
#plt.plot(model.volume_change)
plt.plot(model.timestep_list)