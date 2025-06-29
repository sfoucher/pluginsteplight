#!/usr/bin/env python3
"""
Basic usage example for PythonStep package
"""

import sys
import os
import numpy as np

# Add the parent directory to the path to import PythonStep
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

try:
    import PythonStep as ps
except ImportError as e:
    print(f"Error importing PythonStep: {e}")
    print("Make sure the package is properly installed and all dependencies are available.")
    sys.exit(1)


def main():
    """Main example function"""
    
    # Path to the twoCylinders.txt file
    data_file = "../twoCylinders.txt"
    
    if not os.path.exists(data_file):
        print(f"Data file not found: {data_file}")
        print("Please make sure twoCylinders.txt is in the parent directory.")
        return
    
    print("=== PythonStep Basic Usage Example ===\n")
    
    # Example 1: Create grid from text file
    print("1. Creating grid from text file...")
    try:
        grid = ps.create_grid_from_txt(data_file, resolution=0.05)
        print(f"   ✓ Grid created successfully!")
        print(f"   Grid dimensions: {grid.dimensions}")
        print(f"   Grid resolution: {grid.resolution}")
        print(f"   Grid bounds: {grid.bounds}")
    except Exception as e:
        print(f"   ✗ Error creating grid: {e}")
        return
    
    # Example 2: Access grid data
    print("\n2. Accessing grid data...")
    try:
        data = grid.get_data()
        print(f"   ✓ Grid data shape: {data.shape}")
        print(f"   Data type: {data.dtype}")
        print(f"   Min value: {data.min()}")
        print(f"   Max value: {data.max()}")
        print(f"   Mean value: {data.mean():.2f}")
    except Exception as e:
        print(f"   ✗ Error accessing grid data: {e}")
    
    # Example 3: Get values at specific coordinates
    print("\n3. Getting values at specific coordinates...")
    try:
        # Get dimensions
        nx, ny, nz = grid.dimensions
        
        # Get value at center of grid
        center_x, center_y, center_z = nx // 2, ny // 2, nz // 2
        value = grid.get_value(center_x, center_y, center_z)
        print(f"   ✓ Value at center ({center_x}, {center_y}, {center_z}): {value}")
        
        # Get value at origin
        value = grid.get_value(0, 0, 0)
        print(f"   ✓ Value at origin (0, 0, 0): {value}")
        
    except Exception as e:
        print(f"   ✗ Error getting values: {e}")
    
    # Example 4: Coordinate conversion
    print("\n4. Coordinate conversion...")
    try:
        # Convert pixel to cartesian
        pixel_coords = (10, 20, 30)
        cartesian_coords = grid.pixel_to_cartesian(*pixel_coords)
        print(f"   ✓ Pixel {pixel_coords} → Cartesian {cartesian_coords}")
        
        # Convert cartesian to pixel
        cartesian_coords = (1.5, 2.3, 3.1)
        pixel_coords = grid.cartesian_to_pixel(*cartesian_coords)
        print(f"   ✓ Cartesian {cartesian_coords} → Pixel {pixel_coords}")
        
    except Exception as e:
        print(f"   ✗ Error in coordinate conversion: {e}")
    
    # Example 5: Read point cloud data
    print("\n5. Reading point cloud data...")
    try:
        points, normals = ps.read_point_cloud_file(data_file)
        print(f"   ✓ Read {len(points)} points from file")
        print(f"   Points shape: {points.shape}")
        print(f"   Normals shape: {normals.shape}")
        
        # Validate data
        ps.validate_point_cloud(points, normals)
        print(f"   ✓ Point cloud data is valid")
        
        # Calculate bounding box
        min_coords, max_coords = ps.calculate_bounding_box(points, padding=0.1)
        print(f"   Bounding box: min={min_coords}, max={max_coords}")
        
    except Exception as e:
        print(f"   ✗ Error reading point cloud: {e}")
    
    # Example 6: Create grid from points
    print("\n6. Creating grid from points...")
    try:
        grid2 = ps.create_grid_from_points(
            points=points,
            normals=normals,
            resolution=0.1,
            padding=0.5,
            na_value=-9999
        )
        print(f"   ✓ Grid created from points successfully!")
        print(f"   Grid dimensions: {grid2.dimensions}")
        
    except Exception as e:
        print(f"   ✗ Error creating grid from points: {e}")
    
    print("\n=== Example completed successfully! ===")


if __name__ == "__main__":
    main() 