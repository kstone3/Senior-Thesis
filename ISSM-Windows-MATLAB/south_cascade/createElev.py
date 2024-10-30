import rasterio
import numpy as np
from scipy.io import savemat
import geopandas as gpd
import pandas as pd

# Load the DEM file
dem_path = '../../Data/SouthCascadeData/SouthCascade_1958.08.13_DEM/SouthCascade_1958.08.13_DEM.tif'
with rasterio.open(dem_path) as dem:
    dem_data = dem.read(1)  # Read the first band
    dem_transform = dem.transform  # Get the transformation info

# Load the polygon shapefile
polygon_path = '../../Data/SouthCascadeData/large_polygon.shp'
polygon = gpd.read_file(polygon_path)

# Reproject polygon to DEM CRS if necessary
if dem.crs != polygon.crs:
    polygon = polygon.to_crs(dem.crs)

# Define resolution
resolution = 100
# Get bounds of the polygon and set x and y ranges to create the grid
minx, miny, maxx, maxy = polygon.total_bounds

# Create x and y points for a square grid
x_points = np.arange(minx, maxx + resolution, resolution)
y_points = np.arange(miny, maxy + resolution, resolution)

# Determine the number of points along x and y
num_points = max(len(x_points), len(y_points))

# Create square grids by extending the shorter dimension
x_points = np.linspace(minx, maxx, num_points)
y_points = np.linspace(miny, maxy, num_points)

# Create a meshgrid of points
x_mesh, y_mesh = np.meshgrid(x_points, y_points)
flat_points = np.c_[x_mesh.ravel(), y_mesh.ravel()]

# Create a GeoDataFrame of all grid points and filter by polygon
grid_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(flat_points[:, 0], flat_points[:, 1]), crs=polygon.crs)
grid_gdf = grid_gdf[grid_gdf.geometry.within(polygon.unary_union)]

# Initialize a square 2D array for elevations within the polygon bounds
elevation_grid = np.full((num_points, num_points), np.nan)

# Assign elevation values to points within the polygon from the DEM data
for point in grid_gdf.geometry:
    row, col = ~dem_transform * (point.x, point.y)
    row, col = int(row), int(col)

    # Check if the point falls within DEM bounds
    if 0 <= row < dem_data.shape[0] and 0 <= col < dem_data.shape[1]:
        y_idx = np.where(y_points == point.y)[0][0]
        x_idx = np.where(x_points == point.x)[0][0]
        elevation_grid[y_idx, x_idx] = dem_data[row, col]

# Save elevation grid as a .mat file
savemat('elevation_grid.mat', {'elevation_grid': elevation_grid})

# Save x and y coordinates as separate CSV files
pd.DataFrame(x_points, columns=['x_coord']).to_csv('x_coordinates.csv', index=False, header=False)
pd.DataFrame(y_points, columns=['y_coord']).to_csv('y_coordinates.csv', index=False, header=False)
print(elevation_grid.shape)
print("Elevation grid saved as 'elevation_grid.mat', x_coords saved as 'x_coordinates.csv', y_coords saved as 'y_coordinates.csv'")
