import geopandas as gpd
import numpy as np

# Load shapefile
#glacier_gdf = gpd.read_file('../../Data/SouthCascadeData/glims_download_58465/glims_polygons.shp')
glacier_gdf=gpd.read_file('Data/glims_polygons_crs3413.shp')
# Ensure it's in a single polygon format (optional)
glacier_outline = glacier_gdf.unary_union

# Create a new GeoDataFrame to hold the transformed geometry

# Transform to target CRS (EPSG:32610)
#glacier_gdf_transformed = glacier_gdf_transformed.to_crs("EPSG:32610")
# Extract boundary coordinates
if glacier_outline.geom_type == 'Polygon':
    coords = np.array(glacier_outline.exterior.coords)
else:
    raise ValueError("Expected a single polygon for glacier outline.")
# Write to .exp file
# coords[:, 1] = -coords[:, 1]
# coords[1,:] = -coords[1,:]
with open('Data/south_cascade_glacier2.exp', 'w') as f:
    f.write("## Name:DomainOutline\n")
    f.write("## Icon:0\n")
    f.write("# Points Count Value\n")
    f.write(f"{len(coords)} 1.000000\n")
    f.write("# X pos Y pos\n")#f.write(f"{len(coords)}\n")  # Write number of points
    for i, coord in enumerate(coords, start=1):
        x, y,z = coord  # Extract x, y, and z
        f.write(f"{x} {y}\n")