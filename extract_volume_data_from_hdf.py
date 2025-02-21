import xarray as xr
import pandas as pd

# Path to your NetCDF file containing melt_on_glacier_monthly
file_path = '/mnt/c/Users/bookn/Downloads/run_hydro_w5e5_gcm_merged_TaiESM1_ssp585_bc_2000_2019_Batch_18000_19000.nc'

# 1. Open the NetCDF dataset using xarray
ds = xr.open_dataset(file_path, engine="netcdf4")
# print(ds)

# # 2. Print all coordinate names and their values:
# print("Dataset coordinates:", list(ds.coords))

# # 3. Print the dimensions (size of each coordinate axis):
# print("Dataset dimensions:", ds.dims)

# # 4. Print the data variables in the dataset:
# print("Data variables:", list(ds.data_vars))

# # 5. Print detail for just the variable "melt_on_glacier_monthly":
# print(ds['melt_on_glacier_monthly'])
# 2. Select the variable of interest: melt_on_glacier_monthly
#    and filter by the specified RGI ID (adjust if the coordinate name is different)
melt_data = ds['melt_on_glacier_monthly'].sel(rgi_id='RGI60-02.18778')

# 3. Convert the xarray DataArray to a pandas DataFrame
melt_df = melt_data.to_dataframe().reset_index()

# 4. Save the DataFrame to a CSV file
output_csv_path = 'melt_data_RGI60-02.18778_2.csv'
melt_df.to_csv(output_csv_path, index=False)

print(f"CSV file saved: {output_csv_path}")
