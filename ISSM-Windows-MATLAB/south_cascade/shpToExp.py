import geopandas as gpd
import numpy as np

# Load shapefile
glacier_outline = gpd.read_file('../../Data/SouthCascadeData/glims_download_58465/glims_polygons.shp')

# Ensure it's in a single polygon format (optional)
glacier_outline = glacier_outline.unary_union

# Extract boundary coordinates
if glacier_outline.geom_type == 'Polygon':
    coords = np.array(glacier_outline.exterior.coords)
else:
    raise ValueError("Expected a single polygon for glacier outline.")

# Write to .exp file
with open('south_cascade_glacier.exp', 'w') as f:
    #f.write(f"{len(coords)}\n")  # Write number of points
    for i, coord in enumerate(coords, start=1):
        x, y, z = coord  # Extract x, y, and z
        f.write(f"{x} {y}\n")