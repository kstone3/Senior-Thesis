import rasterio
import pandas as pd
import numpy as np
from scipy.io import savemat
import geopandas as gpd
from shapely.geometry import Polygon, box
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def parse_exp_file(exp_file):
    coords = []
    
    with open(exp_file, 'r') as f:
        for line in f:
            point = tuple(map(float, line.strip().split()))
            coords.append(point)

    if len(coords) >= 4:
        if coords[0] != coords[-1]:
            coords.append(coords[0])
        return gpd.GeoSeries([Polygon(coords)])
    else:
        print(f"Not enough coordinates to form a valid polygon.")
        return gpd.GeoSeries([])

def extract_elevations(dem_file, exp_file, resolution=100):
    with rasterio.open(dem_file) as src:
        dem_data = src.read(1)
        dem_transform = src.transform

    polygons = parse_exp_file(exp_file)

    # Ensure CRS is set
    if polygons.crs is None:
        polygons.set_crs("EPSG:32610", inplace=True)

    # Check polygon validity
    if not polygons.is_valid.all():
        print("Polygon is invalid.")

    print(f"Polygon: {polygons[0]}")
    print(f"Polygon Bounds: {polygons.total_bounds}")

    minx, miny, maxx, maxy = polygons.total_bounds
    x_coords = np.arange(minx, maxx, resolution)
    y_coords = np.arange(miny, maxy, resolution)

    elevation_grid = np.full((len(y_coords), len(x_coords)), 1000.0)
    x_output = x_coords + resolution / 2
    y_output = y_coords + resolution / 2

    # Visualization
    fig, ax = plt.subplots()
    polygons.plot(ax=ax, facecolor='none', edgecolor='red')

    for i,y in enumerate(y_coords):
        for j,x in enumerate(x_coords):
            square = box(x, y, x + resolution, y + resolution)

            # Buffering polygon to handle precision issues
            buffered_polygon = polygons.buffer(0)

            if buffered_polygon.intersects(square).any():
                row, col = ~dem_transform * (x + resolution / 2, y + resolution / 2)
                row, col = int(row), int(col)

                if 0 <= row < dem_data.shape[0] and 0 <= col < dem_data.shape[1]:
                    elevation_value = dem_data[row, col]
                    elevation_grid[i, j] = elevation_value
                    #if ~np.isin((x+resolution/2),x_output) and ~np.isin((y+resolution/2),y_output): 
                    # x_output.append(x + resolution / 2)
                    # y_output.append(y + resolution / 2)

                    # Plot square
                    plt.gca().add_patch(plt.Rectangle((x, y), resolution, resolution, fill=None, edgecolor='blue'))

                else:
                    print(f"Row/Col out of bounds: Row: {row}, Col: {col}")
            else:
                print(f"Square at ({x}, {y}) does not intersect with polygon.")

    print("elevation_grid shape: ",elevation_grid.shape)
    print("x_output shape: ", x_output.shape)
    print("y_output shape: ", y_output.shape)
    plt.xlim(minx, maxx)
    plt.ylim(miny, maxy)
    plt.title('Polygon and Squares Visualization')
    #plt.show()
    # Create DataFrames for output
    #elevation_grid=elevation_grid.astype(int)
    np.savetxt('elevation_grid.csv', elevation_grid, delimiter=',', comments='')
    x_coords_df = pd.DataFrame({'X_Coordinate': x_output})
    y_coords_df = pd.DataFrame({'Y_Coordinate': y_output})

    # Save to CSV
    x_coords_df.to_csv('x_coordinates.csv', index=False, header=False)
    y_coords_df.to_csv('y_coordinates.csv', index=False, header=False)

    # Save to MAT file
    savemat('elevation_grid.mat', {'elevations': elevation_grid})

# Example usage
dem_file = '../../Data/SouthCascadeData/SouthCascade_1958.08.13_DEM/SouthCascade_1958.08.13_DEM.tif'
polygon_file = 'south_cascade_glacier_py.exp'
extract_elevations(dem_file, polygon_file)
