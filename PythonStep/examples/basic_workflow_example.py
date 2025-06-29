#!/usr/bin/env python3
"""
Basic Workflow Example using twoCylinders.txt

This example demonstrates a complete workflow for processing point cloud data:
1. Load point cloud from twoCylinders.txt
2. Create a 3D grid
3. Apply various filtering and analysis workflows
4. Visualize results

Requirements:
- numpy
- scipy
- matplotlib (for visualization)
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Add the parent directory to the path to import PythonStep modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from PythonStep.utils import read_point_cloud_file, calculate_bounding_box, validate_point_cloud
from PythonStep.workflows import GridWorkflows, GridVisualization


def main():
    print("=== PythonStep Basic Workflow Example ===\n")
    
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
    
    # Step 2: Create a synthetic grid (since we don't have the full ray tracing implementation)
    print("\n2. Creating synthetic grid for demonstration...")
    
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
    # Find clusters in the point cloud to simulate cylinder locations
    from scipy.cluster.vq import kmeans, vq
    
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
    
    # Step 3: Apply various workflows
    print("\n3. Applying Python workflows...")
    
    # 3.1 Grid statistics
    print("   3.1 Computing grid statistics...")
    stats = GridWorkflows.grid_statistics(grid_data)
    print(f"   - Total cells: {stats['total_cells']}")
    print(f"   - Non-zero cells: {stats['non_zero_cells']}")
    print(f"   - Sparsity: {stats['sparsity']:.2%}")
    print(f"   - Value range: {stats['min_value']:.2f} to {stats['max_value']:.2f}")
    
    # 3.2 Threshold filtering
    print("   3.2 Applying threshold filtering...")
    threshold = np.percentile(grid_data[grid_data > 0], 75)  # 75th percentile
    filtered_threshold = GridWorkflows.threshold_based_filtering(grid_data, threshold)
    print(f"   - Threshold: {threshold:.2f}")
    print(f"   - Remaining cells: {np.count_nonzero(filtered_threshold)}")
    
    # 3.3 Neighbor-based filtering (local maxima)
    print("   3.3 Applying neighbor-based filtering...")
    filtered_neighbors = GridWorkflows.neighbor_based_filtering(grid_data, neighbors=2)
    print(f"   - Remaining cells: {np.count_nonzero(filtered_neighbors)}")
    
    # 3.4 Local maxima detection
    print("   3.4 Detecting local maxima...")
    local_maxima = GridWorkflows.local_maxima_detection(
        grid_data, neighborhood_size=2, min_value=threshold, sort_descending=True
    )
    print(f"   - Found {len(local_maxima)} local maxima")
    
    if local_maxima:
        print("   - Top 5 maxima:")
        for i, (x, y, z, value) in enumerate(local_maxima[:5]):
            print(f"     {i+1}. Position ({x}, {y}, {z}), Value: {value:.2f}")
    
    # 3.5 Gaussian smoothing
    print("   3.5 Applying Gaussian smoothing...")
    smoothed = GridWorkflows.gaussian_filtering(grid_data, sigma=1.0)
    
    # 3.6 Connected components analysis
    print("   3.6 Connected components analysis...")
    labeled_grid, num_components = GridWorkflows.connected_components_analysis(
        grid_data, min_size=10
    )
    print(f"   - Found {num_components} connected components")
    
    # Step 4: Visualization
    print("\n4. Creating visualizations...")
    
    try:
        # Create a figure with multiple subplots
        fig = plt.figure(figsize=(15, 10))
        
        # Original grid (middle slice)
        ax1 = fig.add_subplot(2, 3, 1)
        slice_data = GridVisualization.create_2d_slice(grid_data, axis=1, slice_index=dimensions[1]//2)
        im1 = ax1.imshow(slice_data.T, cmap='viridis', origin='lower')
        ax1.set_title('Original Grid (Y-slice)')
        plt.colorbar(im1, ax=ax1)
        
        # Threshold filtered
        ax2 = fig.add_subplot(2, 3, 2)
        slice_filtered = GridVisualization.create_2d_slice(filtered_threshold, axis=1, slice_index=dimensions[1]//2)
        im2 = ax2.imshow(slice_filtered.T, cmap='viridis', origin='lower')
        ax2.set_title('Threshold Filtered')
        plt.colorbar(im2, ax=ax2)
        
        # Neighbor filtered
        ax3 = fig.add_subplot(2, 3, 3)
        slice_neighbors = GridVisualization.create_2d_slice(filtered_neighbors, axis=1, slice_index=dimensions[1]//2)
        im3 = ax3.imshow(slice_neighbors.T, cmap='viridis', origin='lower')
        ax3.set_title('Neighbor Filtered')
        plt.colorbar(im3, ax=ax3)
        
        # Smoothed
        ax4 = fig.add_subplot(2, 3, 4)
        slice_smoothed = GridVisualization.create_2d_slice(smoothed, axis=1, slice_index=dimensions[1]//2)
        im4 = ax4.imshow(slice_smoothed.T, cmap='viridis', origin='lower')
        ax4.set_title('Gaussian Smoothed')
        plt.colorbar(im4, ax=ax4)
        
        # Connected components
        ax5 = fig.add_subplot(2, 3, 5)
        slice_components = GridVisualization.create_2d_slice(labeled_grid, axis=1, slice_index=dimensions[1]//2)
        im5 = ax5.imshow(slice_components.T, cmap='tab10', origin='lower')
        ax5.set_title('Connected Components')
        plt.colorbar(im5, ax=ax5)
        
        # 3D scatter plot of local maxima
        ax6 = fig.add_subplot(2, 3, 6, projection='3d')
        if local_maxima:
            maxima_coords = np.array(local_maxima)[:, :3]
            maxima_values = np.array(local_maxima)[:, 3]
            scatter = ax6.scatter(maxima_coords[:, 0], maxima_coords[:, 1], maxima_coords[:, 2], 
                                c=maxima_values, cmap='viridis', s=50)
            ax6.set_title('Local Maxima (3D)')
            plt.colorbar(scatter, ax=ax6)
        
        plt.tight_layout()
        
        # Save the figure
        output_file = "workflow_results.png"
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"   ✓ Visualization saved as {output_file}")
        
        # Show the plot
        plt.show()
        
    except ImportError:
        print("   ⚠ matplotlib not available, skipping visualization")
    except Exception as e:
        print(f"   ⚠ Visualization error: {e}")
    
    # Step 5: Summary
    print("\n5. Workflow Summary:")
    print(f"   - Input: {len(points)} points from {data_file}")
    print(f"   - Grid: {grid_data.shape} with resolution {resolution}")
    print(f"   - Processing: Applied {len(local_maxima)} local maxima detection")
    print(f"   - Output: Filtered grids and analysis results")
    
    print("\n=== Workflow completed successfully! ===")


if __name__ == "__main__":
    main() 