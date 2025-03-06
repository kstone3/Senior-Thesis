import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pandas as pd
import numpy as np
from math import inf
import matplotlib
from glacierSim import glacierSim
matplotlib.use("TkAgg")

#Initialize plotting stuff outside of class
fig, ax = plt.subplots() #initialize plotting variables
_ = plt.close(fig) #used to prevent an empty plot from displaying

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
        print(input_params)
        #print("GAMMA: ", input_params[0])
        ela = 1880
        ela_1900=1910
        time = 540
        save = 1 #Needs to be 1 to calculate ela list
        gamma=0.016
        lapse_rate=[-0.00334774, -0.00544884, -0.00577458, -0.00679377, -0.00661499, -0.00627995, -0.00529508, -0.00534911, -0.00495446, -0.00494315, -0.00472614, -0.00452499]
        tune_factors=[input_params[0],input_params[1],lapse_rate,0.66,1.7,0.32,1.58]
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
        model = glacierSim(ela=ela, ela_1900=ela_1900,time=time, save=save,gamma=gamma,quiet=quiet, tune_factors=tune_factors, initial_ice=ice, start_time=start_time)
        model.init(ela=ela, ela_1900=ela_1900,time=time, save=save,gamma=gamma,quiet=quiet, tune_factors=tune_factors, initial_ice=ice, start_time=start_time)
        for i in range(0,model.frames):
            model.run_model(i)
        try:
            if parameter == 'summer':
                summer_min=min(np.nanmin(model.calculated_summer_mb),np.nanmin(model.summer_mb))
                summer_max=max(np.nanmax(model.calculated_summer_mb),np.nanmax(model.summer_mb))
                calc_summer_mb_norm=(model.calculated_summer_mb-summer_min)/(summer_max-summer_min)
                meas_summer_mb_norm= (np.array(model.summer_mb)-summer_min)/(summer_max-summer_min)
                mse=np.mean((calc_summer_mb_norm-meas_summer_mb_norm)**2)
                print("Function result: ", mse)
                return mse
            elif parameter == 'winter':
                winter_min=min(np.nanmin(model.calculated_winter_mb),np.nanmin(model.winter_mb))
                winter_max=max(np.nanmax(model.calculated_winter_mb),np.nanmax(model.winter_mb))
                calc_winter_mb_norm=(model.calculated_winter_mb-winter_min)/(winter_max-winter_min)
                meas_winter_mb_norm= (np.array(model.winter_mb)-winter_min)/(winter_max-winter_min)
                mse=np.mean((calc_winter_mb_norm-meas_winter_mb_norm)**2)
                print("Function result: ", mse)
                return mse
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
            else: raise ValueError("Invalid parameter. Choose from 'summer', 'winter', 'annual', 'summer_winter', 'vol_change'.")
        except Exception as e:
            print("Error during function calculation:", e) 
            print("Function result: INF")
            return inf
    except Exception as e:
        print("Error during model execution:", e) 
        print("Function result: INF")
        return inf

accumfactor=0.01 #bounds approx 0.001-0.005
ice_meltfactor= -0.01 #bounds approx 0.005-0.012
snow_meltfactor=-0.01 #bounds approx 0.002-0.006
snow_conv_factor=6#bounds 5-15
snow_melt_amplitude=-0.001
ice_melt_amplitude=-0.001
bounds=[(-1,0),(-1,0)]
initial_guess=[-0.0039,-0.0025]
opt_method='Nelder-Mead'
with open(f"../Results/{opt_method}-Results.txt", "a") as file:
    file.write("-------------------Optimize summer mb-----------\n")
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
