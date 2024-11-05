import rasterio
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt

from shapely.geometry import Point, Polygon

# Load the DEM file
dem_path = '../../Data/SouthCascadeData/SouthCascade_1958.08.13_DEM/SouthCascade_1958.08.13_DEM.tif'
with rasterio.open(dem_path) as dem:
    dem_data = dem.read(1)  # Read the first band
    dem_transform = dem.transform  # Get the transformation info
    dem_crs = dem.crs  # Get the CRS of the DEM

# Load the polygon from the .exp file
exp_file_path = 'south_cascade_glacier.exp'

# Read the coordinates from the .exp file
coords = pd.read_csv(exp_file_path, sep=' ', header=None, names=['x', 'y'])

# Create a list of Points from the coordinates
points = [Point(xy) for xy in zip(coords['x'], coords['y'])]

# Create a Polygon from the Points
if len(points) > 2:  # Need at least 3 points to form a polygon
    polygon = Polygon(points)
else:
    raise ValueError("Not enough points to create a polygon.")

# Set the CRS for the polygon (replace with actual CRS if known)
polygon_crs = "EPSG:32610"  # WGS84 for lat/lon coordinates
polygon_gdf = gpd.GeoDataFrame(geometry=[polygon], crs=polygon_crs)

# Reproject polygon to DEM CRS if necessary
if dem_crs != polygon_gdf.crs:
    print("Reprojecting polygon to DEM CRS.")
    polygon_gdf = polygon_gdf.to_crs(dem_crs)

# Define resolution in degrees (approx. 100 meters in latitude)
resolution = 0.000899  # Approximate degree change for 100 meters

# Get bounds of the polygon and set x and y ranges to create the grid
minx, miny, maxx, maxy = polygon_gdf.total_bounds
print("Polygon bounds:", polygon_gdf.total_bounds)

# Create x and y points for a grid
x_points = np.arange(minx, maxx + resolution, resolution)
y_points = np.arange(miny, maxy + resolution, resolution)
print(x_points)
# Create a meshgrid of points
x_mesh, y_mesh = np.meshgrid(x_points, y_points)
flat_points = np.c_[x_mesh.ravel(), y_mesh.ravel()]

# Create a GeoDataFrame of all grid points and filter by polygon
grid_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(flat_points[:, 0], flat_points[:, 1]), crs=polygon_gdf.crs)
grid_gdf = grid_gdf[grid_gdf.geometry.within(polygon_gdf.unary_union)]
print(grid_gdf)
# Initialize a square 2D array for elevations within the polygon bounds
elevation_grid = np.full((len(y_points), len(x_points)), np.nan)

# Assign elevation values to points within the polygon from the DEM data
for point in grid_gdf.geometry:
    row, col = ~dem_transform * (point.x, point.y)
    row, col = int(row), int(col)

    # Check if the point falls within DEM bounds
    if 0 <= row < dem_data.shape[0] and 0 <= col < dem_data.shape[1]:
        # Find indices of y and x coordinates
        y_idx = np.where(y_points == point.y)[0][0]
        x_idx = np.where(x_points == point.x)[0][0]
        elevation_grid[y_idx, x_idx] = dem_data[row, col]

# Remove rows and columns that are entirely NaN from elevation_grid
non_nan_rows = ~np.isnan(elevation_grid).all(axis=1)
non_nan_cols = ~np.isnan(elevation_grid).all(axis=0)
elevation_grid = elevation_grid[non_nan_rows][:, non_nan_cols]
print(dem_data)
# Update x_points and y_points to match the cleaned elevation_grid
x_points = x_points[non_nan_cols]
y_points = y_points[non_nan_rows]

# Save the elevation grid as a .mat file
elevation_grid_filled = np.nan_to_num(elevation_grid, nan=0.0)

# Convert the elevation grid to a DataFrame
elevation_df = pd.DataFrame(elevation_grid_filled)

# Save the DataFrame to a CSV file
elevation_df.to_csv('elevation_grid.csv', index=False, header=False)
# Save x and y coordinates as separate CSV files
pd.DataFrame(x_points, columns=['x_coord']).to_csv('x_coordinates.csv', index=False, header=False)
pd.DataFrame(y_points, columns=['y_coord']).to_csv('y_coordinates.csv', index=False, header=False)

print("Elevation grid saved as 'elevation_grid.csv', x_coords saved as 'x_coordinates.csv', y_coords saved as 'y_coordinates.csv'")


#from osgeo import gdal, osr
# import numpy as np
# from scipy.io import savemat
# import geopandas as gpd
# import pandas as pd

# # Load the DEM file using GDAL
# dem_path = '../../Data/SouthCascadeData/SouthCascade_1958.08.13_DEM/SouthCascade_1958.08.13_DEM.tif'
# dem = gdal.Open(dem_path)
# dem_data = dem.ReadAsArray()  # Read the DEM as an array
# dem_transform = dem.GetGeoTransform()  # Get the transformation info

# # Load the polygon shapefile
# polygon_path = '../../Data/SouthCascadeData/large_polygon.shp'
# polygon = gpd.read_file(polygon_path)

# # Get the CRS from the polygon and set up the DEM's projection if necessary
# polygon_crs = polygon.crs.to_string()  # Get CRS as a string
# dem_crs = osr.SpatialReference(wkt=dem.GetProjection())

# # Reproject DEM to match polygon CRS if needed
# if dem_crs.ExportToProj4() != polygon_crs:
#     # Create a destination projection from polygon CRS
#     target_srs = osr.SpatialReference()
#     target_srs.ImportFromProj4(polygon_crs)

#     # Set up the reprojection
#     reprojected_dem_path = 'reprojected_dem.tif'
#     gdal.Warp(reprojected_dem_path, dem, dstSRS=polygon_crs, resampleAlg=gdal.GRA_Bilinear)

#     # Reopen the reprojected DEM
#     dem = gdal.Open(reprojected_dem_path)
#     dem_data = dem.ReadAsArray()
#     dem_transform = dem.GetGeoTransform()

# # Define resolution
# resolution = 100
# # Get bounds of the polygon and set x and y ranges to create the grid
# minx, miny, maxx, maxy = polygon.total_bounds

# # Create x and y points for a square grid
# x_points = np.arange(minx, maxx + resolution, resolution)
# y_points = np.arange(miny, maxy + resolution, resolution)

# # Determine the number of points along x and y
# num_points = max(len(x_points), len(y_points))

# # Create square grids by extending the shorter dimension
# x_points = np.linspace(minx, maxx, num_points)
# y_points = np.linspace(miny, maxy, num_points)

# # Create a meshgrid of points
# x_mesh, y_mesh = np.meshgrid(x_points, y_points)
# flat_points = np.c_[x_mesh.ravel(), y_mesh.ravel()]

# # Create a GeoDataFrame of all grid points and filter by polygon
# grid_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(flat_points[:, 0], flat_points[:, 1]), crs=polygon.crs)
# grid_gdf = grid_gdf[grid_gdf.geometry.within(polygon.unary_union)]

# # Initialize a square 2D array for elevations within the polygon bounds
# elevation_grid = np.full((num_points, num_points), np.nan)

# # Assign elevation values to points within the polygon from the DEM data
# for point in grid_gdf.geometry:
#     # Get the pixel coordinates of the point in the reprojected DEM
#     x, y = point.x, point.y
#     col = int((x - dem_transform[0]) / dem_transform[1])
#     row = int((y - dem_transform[3]) / dem_transform[5])

#     # Check if the point falls within DEM bounds
#     if 0 <= row < dem_data.shape[0] and 0 <= col < dem_data.shape[1]:
#         y_idx = np.where(y_points == point.y)[0][0]
#         x_idx = np.where(x_points == point.x)[0][0]
#         elevation_grid[y_idx, x_idx] = dem_data[row, col]

# # Remove rows and columns that are entirely NaN from elevation_grid
# non_nan_rows = ~np.isnan(elevation_grid).all(axis=1)
# non_nan_cols = ~np.isnan(elevation_grid).all(axis=0)
# elevation_grid = elevation_grid[non_nan_rows][:, non_nan_cols]

# # Update x_points and y_points to match the cleaned elevation_grid
# x_points = x_points[non_nan_cols]
# y_points = y_points[non_nan_rows]

# # Save the elevation grid as a .mat file
# savemat('elevation_grid.mat', {'elevation_grid': elevation_grid})
# elevation_grid_filled = np.nan_to_num(elevation_grid, nan=0.0)

# # Convert the elevation grid to a DataFrame
# elevation_df = pd.DataFrame(elevation_grid_filled)

# # Save the DataFrame to a CSV file
# elevation_df.to_csv('elevation_grid.csv', index=False, header=False)
# # Save x and y coordinates as separate CSV files
# pd.DataFrame(x_points, columns=['x_coord']).to_csv('x_coordinates.csv', index=False, header=False)
# pd.DataFrame(y_points, columns=['y_coord']).to_csv('y_coordinates.csv', index=False, header=False)

# print("Elevation grid saved as 'elevation_grid.mat', x_coords saved as 'x_coordinates.csv', y_coords saved as 'y_coordinates.csv'")
