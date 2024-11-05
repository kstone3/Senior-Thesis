import geopandas as gpd
import numpy as np

# Load shapefile
glacier_gdf = gpd.read_file('../../Data/SouthCascadeData/glims_download_58465/glims_polygons.shp')

# Ensure it's in a single polygon format (optional)
glacier_outline = glacier_gdf.unary_union

# Create a new GeoDataFrame to hold the transformed geometry
glacier_gdf_transformed = gpd.GeoDataFrame(geometry=[glacier_outline], crs=glacier_gdf.crs)

# Transform to target CRS (EPSG:32610)
glacier_gdf_transformed = glacier_gdf_transformed.to_crs("EPSG:32610")
# Extract boundary coordinates
if glacier_gdf_transformed.geometry[0].geom_type == 'Polygon':
    coords = np.array(glacier_gdf_transformed.geometry[0].exterior.coords)
else:
    raise ValueError("Expected a single polygon for glacier outline.")

# Write to .exp file
with open('south_cascade_glacier.exp', 'w') as f:
    #f.write(f"{len(coords)}\n")  # Write number of points
    for i, coord in enumerate(coords, start=1):
        x, y, z = coord  # Extract x, y, and z
        f.write(f"{x} {y}\n")