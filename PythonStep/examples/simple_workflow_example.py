#!/usr/bin/env python3
"""
Simple Workflow Example using twoCylinders.txt (No Visualization)

This example demonstrates the core workflow without visualization dependencies:
1. Load point cloud from twoCylinders.txt
2. Create a 3D grid
3. Apply various filtering and analysis workflows
4. Print results

Requirements:
- numpy
- scipy
"""

import os
import sys
import numpy as np

# Add the parent directory to the path to import PythonStep modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from PythonStep.utils import read_point_cloud_file, calculate_bounding_box, validate_point_cloud
from PythonStep.workflows import GridWorkflows


def main():
    print("=== PythonStep Simple Workflow Example ===\n")
    
    # File path
    data_file = "/home/sfoucher/DEV/pluginsteplight/twoCylinders.txt"
    
    if not os.path.exists(data_file):
        print(f"Error: File not found: {data_file}")
        return
    
    # Step 1: Load point cloud data
    print("1. Loading point cloud data...")
    try:
        points, normals = read_point_cloud_file(data_file)
        print(f"   Loaded {len(points)} points")
        print(f"   Point cloud bounds: {points.min(axis=0)} to {points.max(axis=0)}")
        
        # Validate the data
        validate_point_cloud(points, normals)
        print("   ✓ Data validation passed")
        
    except Exception as e:
        print(f"   Error loading data: {e}")
        return
    
    # Step 2: Create a synthetic grid
    print("\n2. Creating synthetic grid...")
    
    # Calculate grid dimensions based on point cloud bounds
    bounds = calculate_bounding_box(points, padding=0.1)
    min_coords, max_coords = bounds
    resolution = 0.1
    
    # Calculate grid dimensions
    dimensions = np.ceil((max_coords - min_coords) / resolution).astype(int)
    print(f"   Grid dimensions: {dimensions}")
    print(f"   Resolution: {resolution}")
    
    # Create a synthetic grid with peaks at cylinder locations
    grid_data = np.zeros(dimensions, dtype=np.float32)
    
    # Add synthetic peaks based on the point cloud structure
    from scipy.cluster.vq import kmeans
    
    # Use k-means to find cluster centers (simulating cylinder centers)
    centroids, _ = kmeans(points, 2)  # 2 cylinders
    
    # Add peaks at cluster centers
    for i, centroid in enumerate(centroids):
        # Convert to grid coordinates
        grid_x = int((centroid[0] - min_coords[0]) / resolution)
        grid_y = int((centroid[1] - min_coords[1]) / resolution)
        grid_z = int((centroid[2] - min_coords[2]) / resolution)
        
        # Ensure coordinates are within bounds
        grid_x = max(0, min(grid_x, dimensions[0] - 1))
        grid_y = max(0, min(grid_y, dimensions[1] - 1))
        grid_z = max(0, min(grid_z, dimensions[2] - 1))
        
        # Add a peak with some spread
        peak_value = 100 + i * 50  # Different heights for each cylinder
        grid_data[grid_x, grid_y, grid_z] = peak_value
        
        # Add some spread around the peak
        for dx in range(-2, 3):
            for dy in range(-2, 3):
                for dz in range(-2, 3):
                    nx, ny, nz = grid_x + dx, grid_y + dy, grid_z + dz
                    if (0 <= nx < dimensions[0] and 
                        0 <= ny < dimensions[1] and 
                        0 <= nz < dimensions[2]):
                        distance = np.sqrt(dx**2 + dy**2 + dz**2)
                        if distance <= 2:
                            grid_data[nx, ny, nz] = max(grid_data[nx, ny, nz], 
                                                       peak_value * np.exp(-distance))
    
    # Add some noise
    np.random.seed(42)
    noise = np.random.normal(0, 5, grid_data.shape)
    grid_data += noise
    grid_data = np.maximum(grid_data, 0)  # Ensure non-negative values
    
    print(f"   ✓ Grid created with shape {grid_data.shape}")
    print(f"   Grid value range: {grid_data.min():.2f} to {grid_data.max():.2f}")
    
    # Step 3: Apply workflows and analyze results
    print("\n3. Applying Python workflows...")
    
    # 3.1 Grid statistics
    print("   3.1 Computing grid statistics...")
    stats = GridWorkflows.grid_statistics(grid_data)
    print(f"   - Total cells: {stats['total_cells']}")
    print(f"   - Non-zero cells: {stats['non_zero_cells']}")
    print(f"   - Sparsity: {stats['sparsity']:.2%}")
    print(f"   - Value range: {stats['min_value']:.2f} to {stats['max_value']:.2f}")
    
    if 'non_zero_percentiles' in stats:
        print(f"   - 75th percentile: {stats['non_zero_percentiles']['75']:.2f}")
        print(f"   - 90th percentile: {stats['non_zero_percentiles']['90']:.2f}")
    
    # 3.2 Threshold filtering
    print("   3.2 Applying threshold filtering...")
    threshold = np.percentile(grid_data[grid_data > 0], 75)  # 75th percentile
    filtered_threshold = GridWorkflows.threshold_based_filtering(grid_data, threshold)
    print(f"   - Threshold: {threshold:.2f}")
    print(f"   - Remaining cells: {np.count_nonzero(filtered_threshold)}")
    print(f"   - Reduction: {100 * (1 - np.count_nonzero(filtered_threshold) / np.count_nonzero(grid_data)):.1f}%")
    
    # 3.3 Neighbor-based filtering (local maxima)
    print("   3.3 Applying neighbor-based filtering...")
    filtered_neighbors = GridWorkflows.neighbor_based_filtering(grid_data, neighbors=2)
    print(f"   - Remaining cells: {np.count_nonzero(filtered_neighbors)}")
    print(f"   - Reduction: {100 * (1 - np.count_nonzero(filtered_neighbors) / np.count_nonzero(grid_data)):.1f}%")
    
    # 3.4 Local maxima detection
    print("   3.4 Detecting local maxima...")
    local_maxima = GridWorkflows.local_maxima_detection(
        grid_data, neighborhood_size=2, min_value=threshold, sort_descending=True
    )
    print(f"   - Found {len(local_maxima)} local maxima")
    
    if local_maxima:
        print("   - Top 10 maxima:")
        for i, (x, y, z, value) in enumerate(local_maxima[:10]):
            print(f"     {i+1:2d}. Position ({x:2d}, {y:2d}, {z:2d}), Value: {value:6.2f}")
    
    # 3.5 Gaussian smoothing
    print("   3.5 Applying Gaussian smoothing...")
    smoothed = GridWorkflows.gaussian_filtering(grid_data, sigma=1.0)
    print(f"   - Smoothed grid range: {smoothed.min():.2f} to {smoothed.max():.2f}")
    
    # 3.6 Connected components analysis
    print("   3.6 Connected components analysis...")
    labeled_grid, num_components = GridWorkflows.connected_components_analysis(
        grid_data, min_size=10
    )
    print(f"   - Found {num_components} connected components")
    
    # 3.7 Distance transform
    print("   3.7 Computing distance transform...")
    distance_transform = GridWorkflows.distance_transform(grid_data, metric='euclidean')
    print(f"   - Distance transform range: {distance_transform.min():.2f} to {distance_transform.max():.2f}")
    
    # Step 4: Summary and results
    print("\n4. Workflow Summary:")
    print(f"   - Input: {len(points)} points from {os.path.basename(data_file)}")
    print(f"   - Grid: {grid_data.shape} with resolution {resolution}")
    print(f"   - Processing: Applied {len(local_maxima)} local maxima detection")
    print(f"   - Filtering: Reduced from {np.count_nonzero(grid_data)} to {np.count_nonzero(filtered_threshold)} cells")
    print(f"   - Components: Found {num_components} connected components")
    
    # Step 5: Save results to file
    print("\n5. Saving results...")
    try:
        # Save local maxima to file
        maxima_file = "local_maxima_results.txt"
        with open(maxima_file, 'w') as f:
            f.write("# Local Maxima Results from PythonStep Workflow\n")
            f.write("# Format: x y z value\n")
            for x, y, z, value in local_maxima:
                f.write(f"{x} {y} {z} {value:.2f}\n")
        print(f"   ✓ Local maxima saved to {maxima_file}")
        
        # Save grid statistics
        stats_file = "grid_statistics.txt"
        with open(stats_file, 'w') as f:
            f.write("# Grid Statistics from PythonStep Workflow\n")
            for key, value in stats.items():
                if isinstance(value, dict):
                    f.write(f"{key}:\n")
                    for subkey, subvalue in value.items():
                        f.write(f"  {subkey}: {subvalue}\n")
                else:
                    f.write(f"{key}: {value}\n")
        print(f"   ✓ Grid statistics saved to {stats_file}")
        
    except Exception as e:
        print(f"   ⚠ Error saving results: {e}")
    
    print("\n=== Workflow completed successfully! ===")
    print("Files created:")
    print("  - local_maxima_results.txt: Coordinates and values of detected peaks")
    print("  - grid_statistics.txt: Comprehensive grid statistics")


if __name__ == "__main__":
    main() 