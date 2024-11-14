import os
from pyproj import Transformer

def load_local_coordinates(filename):
    """Load local coordinates from a space-separated text file, including Z values."""
    local_coords = []
    with open(filename, 'r') as file:
        for line in file:
            if line.strip():  # Skip empty lines
                parts = line.strip().split()

                if len(parts) >= 3:
                    x_local = float(parts[0])
                    y_local = float(parts[1])
                    z_value = float(parts[2])  # Preserve Z value
                    local_coords.append((x_local, y_local, z_value))
                else:
                    print(f"Warning: Skipping incomplete line: {line.strip()}")
    return local_coords

def transform_coordinates(local_coords):
    """Transform local coordinates to EPSG:3413 coordinates, preserving Z."""
    transformed_coords = []
    # Define the transformer from UTM Zone 10N to EPSG:3413
    transformer_to_crs3413 = Transformer.from_crs("epsg:32610", "epsg:3413", always_xy=True)

    # Transformation parameters
    scale_factor = 0.99985
    easting_offset = 641900  # meters
    northing_offset = 5355300  # meters

    for x_local, y_local, z_value in local_coords:
        # Compute UTM Easting and Northing
        easting = x_local * scale_factor + easting_offset
        northing = y_local * scale_factor + northing_offset

        # Convert UTM coordinates to EPSG:3413 coordinates
        x_crs3413, y_crs3413 = transformer_to_crs3413.transform(easting, northing)

        transformed_coords.append((x_crs3413, y_crs3413, z_value))
    return transformed_coords

def write_converted_coordinates(filename, transformed_coords):
    """Write the transformed coordinates to a text file, including Z values."""
    with open(filename, 'w') as file:
        file.write("X_CRS3413\tY_CRS3413\tZ\n")
        for x_crs3413, y_crs3413, z_value in transformed_coords:
            file.write(f"{x_crs3413}\t{y_crs3413}\t{z_value}\n")

def main():
    input_filename = '../../Data/SouthCascadeData/SCG_USGS_100m_grid_xyz_bed.txt'
    output_filename = 'converted_coordinates.txt'

    if not os.path.exists(input_filename):
        print(f"Input file '{input_filename}' does not exist.")
        return

    # Load local coordinates
    local_coords = load_local_coordinates(input_filename)
    if not local_coords:
        print("No valid local coordinates found.")
        return

    # Transform coordinates
    transformed_coords = transform_coordinates(local_coords)

    # Write transformed coordinates to output file
    write_converted_coordinates(output_filename, transformed_coords)
    print(f"Converted coordinates have been saved to '{output_filename}'.")

if __name__ == "__main__":
    main()
