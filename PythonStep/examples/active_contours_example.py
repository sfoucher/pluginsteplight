#!/usr/bin/env python3
"""
Example script demonstrating the Active Contours workflow for extracting characteristic curves
from 3D grid data.

This example shows how to:
1. Create a synthetic 3D grid with curve-like structures
2. Apply the active contours workflow to extract curves
3. Visualize the results
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os

# Add the parent directory to the path to import PythonStep
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from PythonStep.workflows import ActiveContoursWorkflow, OpenActiveContours
from PythonStep.utils import read_point_cloud_file, create_grid_from_points


def create_synthetic_curve_data():
    """
    Create synthetic 3D grid data with curve-like structures
    
    Returns:
        tuple: (grid_data, resolution, true_curves)
    """
    print("Creating synthetic curve data...")
    
    # Grid parameters
    grid_size = (30, 30, 30)
    resolution = 1.0
    
    # Initialize grid
    grid_data = np.zeros(grid_size)
    
    # Create multiple curves
    curves = []
    
    # Curve 1: Straight line
    curve1_points = []
    for t in np.linspace(5, 25, 20):
        x, y, z = int(t), int(10), int(15)
        curve1_points.append([x, y, z])
        # Add high values around the curve
        for dx in range(-2, 3):
            for dy in range(-2, 3):
                for dz in range(-2, 3):
                    nx, ny, nz = x + dx, y + dy, z + dz
                    if (0 <= nx < grid_size[0] and 0 <= ny < grid_size[1] and 0 <= nz < grid_size[2]):
                        dist = np.sqrt(dx**2 + dy**2 + dz**2)
                        if dist <= 2:
                            grid_data[nx, ny, nz] = max(grid_data[nx, ny, nz], 
                                                       50 * np.exp(-dist**2/2))
    curves.append(np.array(curve1_points))
    
    # Curve 2: Helix
    curve2_points = []
    for t in np.linspace(0, 4*np.pi, 30):
        x = int(15 + 8*np.cos(t))
        y = int(15 + 8*np.sin(t))
        z = int(5 + t*2)
        if (0 <= x < grid_size[0] and 0 <= y < grid_size[1] and 0 <= z < grid_size[2]):
            curve2_points.append([x, y, z])
            # Add high values around the curve
            for dx in range(-2, 3):
                for dy in range(-2, 3):
                    for dz in range(-2, 3):
                        nx, ny, nz = x + dx, y + dy, z + dz
                        if (0 <= nx < grid_size[0] and 0 <= ny < grid_size[1] and 0 <= nz < grid_size[2]):
                            dist = np.sqrt(dx**2 + dy**2 + dz**2)
                            if dist <= 2:
                                grid_data[nx, ny, nz] = max(grid_data[nx, ny, nz], 
                                                           40 * np.exp(-dist**2/2))
    curves.append(np.array(curve2_points))
    
    # Curve 3: S-curve
    curve3_points = []
    for t in np.linspace(0, 2*np.pi, 25):
        x = int(5 + 10*np.sin(t))
        y = int(20 + 5*np.cos(t))
        z = int(10 + 5*np.sin(2*t))
        if (0 <= x < grid_size[0] and 0 <= y < grid_size[1] and 0 <= z < grid_size[2]):
            curve3_points.append([x, y, z])
            # Add high values around the curve
            for dx in range(-2, 3):
                for dy in range(-2, 3):
                    for dz in range(-2, 3):
                        nx, ny, nz = x + dx, y + dy, z + dz
                        if (0 <= nx < grid_size[0] and 0 <= ny < grid_size[1] and 0 <= nz < grid_size[2]):
                            dist = np.sqrt(dx**2 + dy**2 + dz**2)
                            if dist <= 2:
                                grid_data[nx, ny, nz] = max(grid_data[nx, ny, nz], 
                                                           35 * np.exp(-dist**2/2))
    curves.append(np.array(curve3_points))
    
    print(f"Created grid with shape {grid_size}")
    print(f"Grid value range: {grid_data.min():.2f} to {grid_data.max():.2f}")
    print(f"Non-zero voxels: {np.sum(grid_data > 0)}")
    
    return grid_data, resolution, curves


def visualize_grid_and_curves(grid_data, extracted_curves, true_curves=None, save_path=None):
    """
    Visualize the grid data and extracted curves
    
    Args:
        grid_data: 3D numpy array
        extracted_curves: List of OpenActiveContours objects
        true_curves: List of true curve points (optional)
        save_path: Path to save the visualization (optional)
    """
    print("Creating visualization...")
    
    # Create figure
    fig = plt.figure(figsize=(15, 10))
    
    # 3D scatter plot
    ax1 = fig.add_subplot(121, projection='3d')
    
    # Plot grid data as scatter points
    x_coords, y_coords, z_coords = np.where(grid_data > 10)
    values = grid_data[x_coords, y_coords, z_coords]
    
    scatter = ax1.scatter(x_coords, y_coords, z_coords, c=values, 
                         cmap='viridis', alpha=0.6, s=1)
    
    # Plot extracted curves
    colors = ['red', 'blue', 'green', 'orange', 'purple']
    for i, curve in enumerate(extracted_curves):
        if len(curve.points) > 0:
            points = curve.get_points_as_array()
            ax1.plot(points[:, 0], points[:, 1], points[:, 2], 
                    color=colors[i % len(colors)], linewidth=3, 
                    label=f'Extracted Curve {i+1}')
    
    # Plot true curves if provided
    if true_curves is not None:
        for i, curve in enumerate(true_curves):
            ax1.plot(curve[:, 0], curve[:, 1], curve[:, 2], 
                    color='black', linestyle='--', linewidth=2, 
                    label=f'True Curve {i+1}')
    
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.set_title('3D Grid Data and Extracted Curves')
    ax1.legend()
    
    # 2D projection plot
    ax2 = fig.add_subplot(122)
    
    # Create 2D projection
    projection = np.max(grid_data, axis=2)
    im = ax2.imshow(projection.T, cmap='viridis', origin='lower')
    
    # Plot curve projections
    for i, curve in enumerate(extracted_curves):
        if len(curve.points) > 0:
            points = curve.get_points_as_array()
            ax2.plot(points[:, 0], points[:, 1], 
                    color=colors[i % len(colors)], linewidth=2, 
                    label=f'Curve {i+1}')
    
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_title('2D Projection (Max Z)')
    ax2.legend()
    
    plt.colorbar(im, ax=ax2)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Visualization saved to {save_path}")
    
    plt.show()


def analyze_curves(extracted_curves):
    """
    Analyze the extracted curves
    
    Args:
        extracted_curves: List of OpenActiveContours objects
    """
    print("\n=== Curve Analysis ===")
    print(f"Number of extracted curves: {len(extracted_curves)}")
    
    for i, curve in enumerate(extracted_curves):
        print(f"\nCurve {i+1}:")
        print(f"  Number of points: {curve.get_n_points()}")
        print(f"  3D length: {curve.length_3d():.2f}")
        
        if len(curve.points) > 0:
            points = curve.get_points_as_array()
            print(f"  Bounding box: X[{points[:, 0].min():.1f}, {points[:, 0].max():.1f}], "
                  f"Y[{points[:, 1].min():.1f}, {points[:, 1].max():.1f}], "
                  f"Z[{points[:, 2].min():.1f}, {points[:, 2].max():.1f}]")


def save_curves_to_file(extracted_curves, filename):
    """
    Save extracted curves to a text file
    
    Args:
        extracted_curves: List of OpenActiveContours objects
        filename: Output filename
    """
    print(f"\nSaving curves to {filename}...")
    
    with open(filename, 'w') as f:
        f.write("# Extracted Active Contours Curves\n")
        f.write("# Format: curve_id point_id x y z\n")
        
        for i, curve in enumerate(extracted_curves):
            if len(curve.points) > 0:
                points = curve.get_points_as_array()
                for j, point in enumerate(points):
                    f.write(f"{i} {j} {point[0]:.6f} {point[1]:.6f} {point[2]:.6f}\n")
    
    print(f"Saved {len(extracted_curves)} curves with total {sum(len(c.points) for c in extracted_curves)} points")


def main():
    """Main function to demonstrate the active contours workflow"""
    print("=== Active Contours Workflow Example ===\n")
    
    # Create synthetic data
    grid_data, resolution, true_curves = create_synthetic_curve_data()
    
    # Create workflow
    print("\nInitializing Active Contours workflow...")
    workflow = ActiveContoursWorkflow(grid_data, resolution)
    
    # Configure workflow parameters
    workflow.n_snakes_max = 5  # Allow up to 5 curves
    workflow.longueur_min = 2.0  # Minimum curve length
    workflow.min_value_for_maxima = 15  # Minimum value for maxima detection
    workflow.maxima_neighborhood_size = 2  # Neighborhood size for maxima detection
    
    print(f"Workflow parameters:")
    print(f"  Max curves: {workflow.n_snakes_max}")
    print(f"  Min length: {workflow.longueur_min}")
    print(f"  Min value for maxima: {workflow.min_value_for_maxima}")
    print(f"  Maxima neighborhood size: {workflow.maxima_neighborhood_size}")
    
    # Extract curves
    print("\nExtracting curves...")
    extracted_curves = workflow.extract_curves()
    
    # Analyze results
    analyze_curves(extracted_curves)
    
    # Visualize results
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    
    visualize_grid_and_curves(
        grid_data, extracted_curves, true_curves,
        save_path=os.path.join(output_dir, "active_contours_visualization.png")
    )
    
    # Save curves to file
    save_curves_to_file(extracted_curves, os.path.join(output_dir, "extracted_curves.txt"))
    
    # Save grid statistics
    stats_file = os.path.join(output_dir, "grid_statistics.txt")
    with open(stats_file, 'w') as f:
        f.write("Grid Statistics:\n")
        f.write(f"Shape: {grid_data.shape}\n")
        f.write(f"Resolution: {resolution}\n")
        f.write(f"Value range: {grid_data.min():.2f} to {grid_data.max():.2f}\n")
        f.write(f"Non-zero voxels: {np.sum(grid_data > 0)}\n")
        f.write(f"Total voxels: {grid_data.size}\n")
        f.write(f"Sparsity: {1.0 - np.sum(grid_data > 0) / grid_data.size:.4f}\n")
    
    print(f"\nGrid statistics saved to {stats_file}")
    print("\n=== Example completed successfully! ===")


if __name__ == "__main__":
    main() 