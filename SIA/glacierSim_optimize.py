import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pandas as pd
import numpy as np
from math import inf
import matplotlib
from glacierSim import glacierSim
import math
matplotlib.use("TkAgg")

def optimize(parameter, input_params):
    try:
        # print("Optimizing: ", parameter)
        #print("ACCUMFACTOR: ", input_params[0])
        print("ICE MELTFACTOR: ", input_params[0])
        print("SNOW MELTFACTOR: ", input_params[1])
        # print("AVALANCHE PERCENTAGE: ", input_params[0])
        # print("SNOW MELT AMPLITUDE: ", input_params[4])
        # print("ICE MELT AMPLITUDE: ", input_params[5])
        # print("LAPSE RATE: ", input_params[1])
        # print("ACCUMFACTOR LOWER: ", input_params[0])
        # print("ACCUMFACTOR UPPER: ", input_params[1])
        # print("PRECIP LAPSE RATE: ", input_params[0])
        # print("GAMMA: ", input_params[0])
        # print("ELA: ", input_params[1])
        # print("ELA 1900: ", input_params[2])
        print(input_params)
        #print("GAMMA: ", input_params[0])
        # gamma=input_params[0] #0.016 #0.012 #0.016
        # ela=input_params[1] #1903 #1903
        # ela_1900=input_params[2] #1930 #1930
        gamma=0.0309
        ela=1903
        ela_1900=1930
        time = 540
        save = 1 #Needs to be 1 to calculate ela list
        # gamma=input_params[0] #0.016 #0.012 #0.016
        # accumfactor_lower=0.66
        # accumfactor_upper=1.7
        accumfactor_lower=0.94 #0.66 #0.66
        accumfactor_upper=1.44 #1.7 #1.7
        # accumfactor_lower=input_params[0] #0.94 #0.66 #0.66
        # accumfactor_upper=input_params[1] #1.44 #1.7 #1.7
        # ice_meltfactor=-0.00314334
        # snow_meltfactor=-0.00315546 #bounds approx 0.002-0.006
        ice_meltfactor=input_params[0] #0.00314334 #0.00315546 #bounds approx 0.002-0.006
        snow_meltfactor=input_params[1] #0.00315546 #0.00315546 #bounds approx 0.002-0.006
        avalanche_percent=0.32
        precip_conv_factor=1.58
        lapse_rate=[-0.00334774, -0.00544884, -0.00577458, -0.00679377, -0.00661499, -0.00627995, -0.00529508, -0.00534911, -0.00495446, -0.00494315, -0.00472614, -0.00452499]
        tune_factors=[ice_meltfactor,snow_meltfactor,lapse_rate,accumfactor_lower,accumfactor_upper,avalanche_percent, precip_conv_factor]
        quiet = True
        start_time = 500
        # ice = [ 53.89550985, 61.2302675, 68.52805603, 72.16752233, 78.19477103,
        #     86.57434438, 94.52278703, 113.00567764, 124.65045342, 131.1336047,
        #     132.61805723, 126.05975829, 117.01403765, 110.72201024, 107.36448442,
        #     113.23158091, 127.19316445, 150.80920841, 147.75800077, 135.19461921,
        #     132.62404257, 130.58089845, 128.95030221, 129.67226275, 131.1275745,
        #     137.76555425, 144.54125286, 149.5876901, 152.88274456, 152.52189309,
        #     142.71479768, 122.66643947, 105.65745331, 89.73132543, 84.64598526,
        #     84.39967405, 78.27781775, 70.24575207, 57.81117783, 43.58883439,
        #     34.64893057, 15.29705107, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
        # ice=[64.50665454, 73.34825753, 82.14677938, 87.1479015, 94.7306883, 104.5882329, 113.70294678, 133.07943541, 144.98064389, 151.46297586, 152.82095652, 146.08286833, 136.91565242, 130.65094348, 127.48133943, 133.66396017, 147.77059641, 171.2134459, 167.71477891, 154.67101123, 151.65240415, 149.12263247, 146.956214, 147.08248738, 147.85770221, 153.72402157, 159.62910372, 163.73136372, 166.02527862, 164.61263013, 153.6902369, 132.39018704, 113.84097541, 95.95253908, 88.15995433, 84.34092379, 73.81852984, 59.65209497, 37.52930102, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0,0]
        ice=np.loadtxt("Data/ice19840101.csv", delimiter=',')
        model = glacierSim(ela=ela, ela_1900=ela_1900,end_time=time, save=save,gamma=gamma,quiet=quiet, tune_factors=tune_factors, initial_ice=ice, start_time=start_time, input_files=input_files)
        model.init(ela=ela, ela_1900=ela_1900,end_time=time, save=save,gamma=gamma,quiet=quiet, tune_factors=tune_factors, initial_ice=ice, start_time=start_time, input_files=input_files)
        for i in range(0,model.frames):
            model.run_model(i)
        # fig, ax = plt.subplots() #initialize plotting variables
        # _ = plt.close(fig) #used to prevent an empty plot from displaying
        # plt.rcParams['animation.embed_limit'] = 40
        # def frame_generator(): 
        #     for i in range(0,model.frames): yield i
        # anim = FuncAnimation(fig, model.func, frames=frame_generator, fargs=(ax,), blit=False, repeat=False, save_count=model.frames)
        # vid = HTML(anim.to_jshtml())
        # display(vid)
        # anim.save("animation.mp4", writer="ffmpeg")
        # plt.plot(model.timestep_list)
        # plt.show()
        try:
            if parameter == 'summer':
                summer_min=min(np.nanmin(model.calculated_summer_mb),np.nanmin(model.summer_mb))
                summer_max=max(np.nanmax(model.calculated_summer_mb),np.nanmax(model.summer_mb))
                calc_summer_mb_norm=(model.calculated_summer_mb-summer_min)/(summer_max-summer_min)
                meas_summer_mb_norm= (np.array(model.summer_mb)-summer_min)/(summer_max-summer_min)
                mse=np.mean((calc_summer_mb_norm-meas_summer_mb_norm)**2)
                mse_non_norm=np.mean((model.summer_mb-model.calculated_summer_mb)**2)
                print("Function result: ", mse_non_norm)
                return mse_non_norm
            elif parameter == 'winter':
                # winter_min=min(np.nanmin(model.calculated_winter_mb),np.nanmin(model.winter_mb))
                # winter_max=max(np.nanmax(model.calculated_winter_mb),np.nanmax(model.winter_mb))
                winter_min=np.nanmin(np.append(model.calculated_winter_mb,model.winter_mb))
                winter_max=np.nanmax(np.append(model.calculated_winter_mb,model.winter_mb))
                calc_winter_mb_norm=(model.calculated_winter_mb-winter_min)/(winter_max-winter_min)
                meas_winter_mb_norm= (np.array(model.winter_mb)-winter_min)/(winter_max-winter_min)
                mse=np.mean((calc_winter_mb_norm-meas_winter_mb_norm)**2)
                mse_non_norm=np.mean((model.winter_mb-model.calculated_winter_mb)**2)
                rmse= math.sqrt(np.sum((model.winter_mb-model.calculated_winter_mb)**2)/len(model.winter_mb))/np.mean(model.winter_mb)
                std=np.std((model.winter_mb-model.calculated_winter_mb))
                print("Function result: ", mse_non_norm)
                return mse_non_norm
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
                df = pd.DataFrame({'date': pd.to_datetime(list(model.date_index_training.keys())),'volume_change': model.daily_runoff_training})
                df['date'] = df['date'].dt.to_period('M')
                train_min=min(np.nanmin(df.groupby('date')['volume_change'].sum().to_numpy()),np.nanmin(model.measured_runoff_training))
                train_max=max(np.nanmax(df.groupby('date')['volume_change'].sum().to_numpy()),np.nanmax(model.measured_runoff_training))
                monthly_volume_change_train_normalized = (df.groupby('date')['volume_change'].sum().to_numpy()-train_min)/(train_max-train_min)
                measured_runoff_train_normalized= (model.measured_runoff_training-train_min)/(train_max-train_min)
                mse=np.mean((monthly_volume_change_train_normalized-measured_runoff_train_normalized)**2)
                print("Function result: ",mse)
                return mse
            elif parameter == 'avalanche_percent':
                print("Function result: ",np.mean(model.snow_depth_list))
                return np.mean(model.snow_depth_list)
            elif parameter == 'spinup':
                meas_area_1958=np.sum(model.thickness_1958_verif[~np.isnan(model.thickness_1958_verif)]*model.dx)
                meas_area_1986=np.sum(model.thickness_1986_verif[~np.isnan(model.thickness_1986_verif)]*model.dx)
                calc_area_1958=np.sum(model.ice_1958[~np.isnan(model.thickness_1958_verif)]*model.dx)
                calc_area_1986=np.sum(model.ice_1986[~np.isnan(model.thickness_1986_verif)]*model.dx)
                glacier_min=min(np.nanmin(model.thickness_1958_verif), np.nanmin(model.thickness_1986_verif), np.nanmin(model.ice_1958), np.nanmin(model.ice_1986))
                glacier_max=max(np.nanmax(model.thickness_1958_verif), np.nanmax(model.thickness_1986_verif), np.nanmax(model.ice_1958), np.nanmax(model.ice_1986))
                ice_norm_1958_verif=(model.thickness_1958_verif[~np.isnan(model.thickness_1958_verif)]-glacier_min)/(glacier_max-glacier_min)
                ice_norm_1986_verif=(model.thickness_1986_verif[~np.isnan(model.thickness_1986_verif)]-glacier_min)/(glacier_max-glacier_min)
                ice_norm_1958_calc=(model.ice_1958[~np.isnan(model.thickness_1958_verif)]-glacier_min)/(glacier_max-glacier_min)
                ice_norm_1986_calc=(model.ice_1986[~np.isnan(model.thickness_1986_verif)]-glacier_min)/(glacier_max-glacier_min)
                mse_1958_norm=np.mean((ice_norm_1958_verif-ice_norm_1958_calc)**2)
                mse_1986_norm=np.mean((ice_norm_1986_verif-ice_norm_1986_calc)**2)
                mse_1958=np.mean((model.thickness_1958_verif[~np.isnan(model.thickness_1958_verif)]-model.ice_1958[~np.isnan(model.thickness_1958_verif)])**2)
                mse_1986=np.mean((model.thickness_1986_verif[~np.isnan(model.thickness_1986_verif)]-model.ice_1986[~np.isnan(model.thickness_1986_verif)])**2)
                mse_1958_area=np.mean((meas_area_1958-calc_area_1958)**2)
                mse_1986_area=np.mean((meas_area_1986-calc_area_1986)**2)
                avg_mse_norm=(mse_1958_norm+mse_1986_norm)/2
                avg_mse=(mse_1958+mse_1986)/2
                avg_mse_area=(mse_1958_area+mse_1986_area)/2
                avg_rmse_area=((math.sqrt((meas_area_1958-calc_area_1958)**2)/meas_area_1958*100)+(math.sqrt((meas_area_1986-calc_area_1986)**2)/meas_area_1986*100))/2
                meas_area_1958=model.thickness_1958_verif[~np.isnan(model.thickness_1958_verif)]*model.dx
                meas_area_1986=model.thickness_1986_verif[~np.isnan(model.thickness_1986_verif)]*model.dx
                calc_area_1958=model.ice_1958[~np.isnan(model.thickness_1958_verif)]*model.dx
                calc_area_1986=model.ice_1986[~np.isnan(model.thickness_1986_verif)]*model.dx
                avg_rmse_area2=((math.sqrt(np.sum((meas_area_1958-calc_area_1958)**2)/len(meas_area_1958))/np.mean(meas_area_1958)*100)+(math.sqrt(np.sum((meas_area_1986-calc_area_1986)**2)/len(meas_area_1986))/np.mean(meas_area_1986)*100))/2
                print("Function result: ",avg_rmse_area2)
                return avg_rmse_area2
            else: raise ValueError("Invalid parameter. Choose from 'summer', 'winter', 'annual', 'summer_winter', 'vol_change'.")
        except Exception as e:
            print("Error during function calculation:", e) 
            print("Function result: INF")
            return inf
    except Exception as e:
        print("Error during model execution:", e) 
        print("Function result: INF")
        return inf

#Initialize the input files variables
input_files={}
input_files['bed']='Data/centerlineBed.csv' #FORMAT: csv with columns, elevation, longitude, latitude
# input_files['bed']='Data/rgi7_centerline_bed.csv'
input_files['area']='Data/Input_SouthCascade_Area_Altitude_Distribution.csv' #FORMAT: csv with first row mean elevation of bin, then columns of year, area per bin
input_files['temp_precip']='Data/Input_SouthCascade_Daily_Weather.csv' #FORMAT: csv with columns date, temperature, precipitation
input_files['mass_balance']='Data/Output_SouthCascade_Glacier_Wide_solutions_calibrated.csv' #FORMAT: csv with columns year, winter mass balance, summer mass balance, annual mass balance, ela
input_files['runoff']='Data/runoff_m3_1992-2007.csv' #FORMAT: csv with columns date, runoff
input_files['thickness_change']='Data/thickness_change.csv' #FORMAT: csv with columns date, thickness change
input_files['front_variation']='Data/front_variation_change.csv' #FORMAT: csv with columns date, front variation change
#Fix basin_area file to represent the new polygon in arcgis
input_files['basin_area']='Data/basin_wide_area_elev_bands.csv' #FORMAT: csv with columns area, mean elevation of bin
input_files['glacier_1958']='Data/centerlineThickness_1958.csv' #FORMAT: csv with columns bed elevation, longitude, latitude, elevation in 1986
input_files['glacier_1986']='Data/centerlineThickness_1986.csv' #FORMAT: csv with columns bed elevation, longitude, latitude, elevation in 1986
input_files['glacier_2021']='Data/centerlineThickness_2021.csv' #FORMAT: csv with columns bed elevation, longitude, latitude, elevation in 2021
# input_files['glacier_1958']='Data/rgi7_centerlineThickness_1958.csv' #FORMAT: csv with columns bed elevation, longitude, latitude, elevation in 1986
# input_files['glacier_1986']='Data/rgi7_centerlineThickness_1986.csv' #FORMAT: csv with columns bed elevation, longitude, latitude, elevation in 1986
# input_files['glacier_2021']='Data/rgi7_centerlineThickness_2021.csv' #FORMAT: csv with columns bed elevation, longitude, latitude, elevation in 2021

ela=1903
ela_1900=1930
run_time=502
start_time=0
save=1 #Needs to be 1 for the ela_list to work properly
gamma=0.0309
quiet=True
accumfactor_lower=0.94
accumfactor_upper=1.44
ice_meltfactor=-0.00314334
snow_meltfactor=-0.00315546
avalanche_percent=0.32
precip_conv_factor=1.58
lapse_rate=[-0.00334774, -0.00544884, -0.00577458, -0.00679377, -0.00661499, -0.00627995, -0.00529508, -0.00534911, -0.00495446, -0.00494315, -0.00472614, -0.00452499]
tune_factors=[ice_meltfactor,snow_meltfactor,lapse_rate,accumfactor_lower,accumfactor_upper,avalanche_percent, precip_conv_factor]
bounds=[(-1,0),(-1,0)] #bounds for ela, ela_1900, gamma
initial_guess=[ice_meltfactor,snow_meltfactor]
opt_method='Nelder-Mead'
with open(f"{opt_method}-Results.txt", "a") as file:
    file.write("-------------------Optimize summer-----------\n")
    # for opt_var in ['ela','front_var','thick','vol_change']:
    for opt_var in ['summer']:
        result = minimize(lambda x: optimize(opt_var, x),initial_guess,method=opt_method,bounds=bounds,options={'disp': True})
        result_snow_conv=result.x
        # print("Final Gamma: ", result.x)
        # print("Final objective function value:", result.fun)
        # print("OPTIMIZED PARAMS FOR:", opt_var)
        # print("Optimized accumulation factor:", result_accum)
        # print("Optimized accumulation factor lower:", result_accum_lower)
        # print("Optimized accumulation factor upper:", result_accumt_upper)
        # print("Optimized ice meltfactor:",result_ice_melt)
        # print("Optimized snow meltfactor:",result_snow_melt)
        print("Optimized snow conversion factor:",result_snow_conv)
        # print("Optimized snow melt amplitude:",result_snow_amp)
        # print("Optimized ice melt amplitude:",result_ice_amp)
        # print("Optimized lapse rate: ", result_lapse_rate)
        print("Final objective function value:", result.fun)
        print(result.x)
        # file.write(f"Optimized parameters for {opt_var}:\n")
        # file.write(f"Optimized accumulation factor: {result_accum}\n")
        # file.write(f"Optimized accumulation factor lower: {result_accum_lower}\n")
        # file.write(f"Optimized accumulation factor upper: {result_accumt_upper}\n")
        # file.write(f"Optimized ice meltfactor: {result_ice_melt}\n")
        # file.write(f"Optimized snow meltfactor: {result_snow_melt}\n")
        file.write(f"Optimized snow conversion factor: {result_snow_conv}\n")
        # file.write(f"Optimized snow melt amplitude: {result_snow_amp}\n")
        # file.write(f"Optimized ice melt amplitude: {result_ice_amp}\n")
        file.write(f"Final objective function value: {result.fun}\n")
        # file.write(f"Optimized lapse rate: {result_lapse_rate}\n")
        file.write(f"{result.x}\n")
        # file.write(f"Final Gamma:  {result.x}\n")
        # file.write(f"Final objective function value: {result.fun}\n")
        file.flush()
