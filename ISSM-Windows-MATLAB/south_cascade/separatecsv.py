import pandas as pd

# Load the CSV file
input_file = 'elevations_within_polygon.csv'  # Replace with your input file name
df = pd.read_csv(input_file)

# Initialize lists to hold the extracted data
x_coords = []
y_coords = []
elevations = []

# Process each row to extract coordinates and elevation
for index, row in df.iterrows():
    # Extract geometry and elevation
    geometry = row['geometry']
    elevation = row['elevation']
    
    # Split the geometry string to extract x and y coordinates
    coords = geometry.replace('POINT (', '').replace(')', '').split()
    x = float(coords[0])
    y = float(coords[1])
    
    # Append the values to the lists
    x_coords.append(x)
    y_coords.append(y)
    elevations.append(elevation)

# Create DataFrames for each coordinate type
x_df = pd.DataFrame(x_coords, columns=['x_coord'])
y_df = pd.DataFrame(y_coords, columns=['y_coord'])
elev_df = pd.DataFrame(elevations, columns=['elevation'])

# Save to separate CSV files
x_df.to_csv('x_coordinates.csv', index=False, header=False)
y_df.to_csv('y_coordinates.csv', index=False, header=False)
elev_df.to_csv('elevations.csv', index=False, header=False)

print("CSV files have been created: x_coordinates.csv, y_coordinates.csv, elevations.csv")
