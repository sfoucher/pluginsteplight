#!/usr/bin/env python3
"""
Example script demonstrating the Active Contours workflow for extracting characteristic curves
from 3D grid data using the twoCylinders.txt point cloud file.

This example shows how to:
1. Load point cloud data from twoCylinders.txt
2. Create a 3D grid from the point cloud
3. Apply the active contours workflow to extract curves
4. Visualize the results
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


def load_two_cylinders_data():
    """
    Load point cloud data from twoCylinders.txt and create a 3D grid
    
    Returns:
        tuple: (grid_data, resolution, point_cloud_info)
    """
    print("Loading twoCylinders.txt point cloud data...")
    
    # Path to the twoCylinders.txt file (assuming it's in the parent directory)
    data_file = os.path.join(os.path.dirname(__file__), '..', '..', 'twoCylinders.txt')
    
    if not os.path.exists(data_file):
        raise FileNotFoundError(f"Could not find twoCylinders.txt at {data_file}")
    
    # Read the point cloud data
    points, normals = read_point_cloud_file(data_file)
    
    print(f"Loaded {len(points)} points from {data_file}")
    print(f"Point cloud bounds:")
    print(f"  X: [{points[:, 0].min():.3f}, {points[:, 0].max():.3f}]")
    print(f"  Y: [{points[:, 1].min():.3f}, {points[:, 1].max():.3f}]")
    print(f"  Z: [{points[:, 2].min():.3f}, {points[:, 2].max():.3f}]")
    
    # Create 3D grid from point cloud
    print("\nCreating 3D grid from point cloud...")
    resolution = 0.5  # Adjust resolution as needed
    grid_data = create_grid_from_points(points, normals, resolution)
    
    # Get grid data as numpy array for analysis
    grid_array = grid_data.get_data()
    print(f"Created grid with shape {grid_array.shape}")
    print(f"Grid value range: {grid_array.min():.2f} to {grid_array.max():.2f}")
    print(f"Non-zero voxels: {np.sum(grid_array > 0)}")
    print(f"Grid sparsity: {1.0 - np.sum(grid_array > 0) / grid_array.size:.4f}")
    
    point_cloud_info = {
        'points': points,
        'normals': normals,
        'bounds': {
            'x_min': points[:, 0].min(),
            'x_max': points[:, 0].max(),
            'y_min': points[:, 1].min(),
            'y_max': points[:, 1].max(),
            'z_min': points[:, 2].min(),
            'z_max': points[:, 2].max()
        }
    }
    
    return grid_data, resolution, point_cloud_info


def visualize_grid_and_curves(grid_data, extracted_curves, point_cloud_info=None, save_path=None):
    """
    Visualize the grid data and extracted curves
    
    Args:
        grid_data: 3D numpy array
        extracted_curves: List of OpenActiveContours objects
        point_cloud_info: Dictionary with point cloud information (optional)
        save_path: Path to save the visualization (optional)
    """
    print("Creating visualization...")
    
    # Use grid_data directly as it's already a numpy array
    grid_array = grid_data
    
    # Create figure
    fig = plt.figure(figsize=(20, 12))
    
    # 3D scatter plot
    ax1 = fig.add_subplot(221, projection='3d')
    
    # Plot grid data as scatter points (only high-value voxels for clarity)
    threshold = np.percentile(grid_array[grid_array > 0], 50) if np.sum(grid_array > 0) > 0 else 1
    x_coords, y_coords, z_coords = np.where(grid_array > threshold)
    values = grid_array[x_coords, y_coords, z_coords]
    
    scatter = ax1.scatter(x_coords, y_coords, z_coords, c=values, 
                         cmap='viridis', alpha=0.6, s=1)
    
    # Plot extracted curves
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'cyan', 'magenta']
    for i, curve in enumerate(extracted_curves):
        if len(curve.points) > 0:
            points = curve.get_points_as_array()
            ax1.plot(points[:, 0], points[:, 1], points[:, 2], 
                    color=colors[i % len(colors)], linewidth=4, 
                    label=f'Extracted Curve {i+1}')
    
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.set_title('3D Grid Data and Extracted Curves')
    ax1.legend()
    
    # 2D projections
    # XY projection
    ax2 = fig.add_subplot(222)
    projection_xy = np.max(grid_array, axis=2)
    im2 = ax2.imshow(projection_xy.T, cmap='viridis', origin='lower')
    
    for i, curve in enumerate(extracted_curves):
        if len(curve.points) > 0:
            points = curve.get_points_as_array()
            ax2.plot(points[:, 0], points[:, 1], 
                    color=colors[i % len(colors)], linewidth=3, 
                    label=f'Curve {i+1}')
    
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_title('XY Projection (Max Z)')
    ax2.legend()
    plt.colorbar(im2, ax=ax2)
    
    # XZ projection
    ax3 = fig.add_subplot(223)
    projection_xz = np.max(grid_array, axis=1)
    im3 = ax3.imshow(projection_xz.T, cmap='viridis', origin='lower')
    
    for i, curve in enumerate(extracted_curves):
        if len(curve.points) > 0:
            points = curve.get_points_as_array()
            ax3.plot(points[:, 0], points[:, 2], 
                    color=colors[i % len(colors)], linewidth=3, 
                    label=f'Curve {i+1}')
    
    ax3.set_xlabel('X')
    ax3.set_ylabel('Z')
    ax3.set_title('XZ Projection (Max Y)')
    ax3.legend()
    plt.colorbar(im3, ax=ax3)
    
    # YZ projection
    ax4 = fig.add_subplot(224)
    projection_yz = np.max(grid_array, axis=0)
    im4 = ax4.imshow(projection_yz.T, cmap='viridis', origin='lower')
    
    for i, curve in enumerate(extracted_curves):
        if len(curve.points) > 0:
            points = curve.get_points_as_array()
            ax4.plot(points[:, 1], points[:, 2], 
                    color=colors[i % len(colors)], linewidth=3, 
                    label=f'Curve {i+1}')
    
    ax4.set_xlabel('Y')
    ax4.set_ylabel('Z')
    ax4.set_title('YZ Projection (Max X)')
    ax4.legend()
    plt.colorbar(im4, ax=ax4)
    
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
        f.write("# Extracted Active Contours Curves from twoCylinders.txt\n")
        f.write("# Format: curve_id point_id x y z\n")
        
        for i, curve in enumerate(extracted_curves):
            if len(curve.points) > 0:
                points = curve.get_points_as_array()
                for j, point in enumerate(points):
                    f.write(f"{i} {j} {point[0]:.6f} {point[1]:.6f} {point[2]:.6f}\n")
    
    print(f"Saved {len(extracted_curves)} curves with total {sum(len(c.points) for c in extracted_curves)} points")


def save_analysis_report(grid_data, extracted_curves, point_cloud_info, filename):
    """
    Save a comprehensive analysis report
    
    Args:
        grid_data: 3D numpy array
        extracted_curves: List of OpenActiveContours objects
        point_cloud_info: Dictionary with point cloud information
        filename: Output filename
    """
    print(f"\nSaving analysis report to {filename}...")
    
    # Use grid_data directly as it's already a numpy array
    grid_array = grid_data
    
    with open(filename, 'w') as f:
        f.write("=== Active Contours Analysis Report ===\n")
        f.write("Data Source: twoCylinders.txt\n\n")
        
        # Point cloud information
        f.write("Point Cloud Information:\n")
        f.write(f"  Number of points: {len(point_cloud_info['points'])}\n")
        f.write(f"  X bounds: [{point_cloud_info['bounds']['x_min']:.3f}, {point_cloud_info['bounds']['x_max']:.3f}]\n")
        f.write(f"  Y bounds: [{point_cloud_info['bounds']['y_min']:.3f}, {point_cloud_info['bounds']['y_max']:.3f}]\n")
        f.write(f"  Z bounds: [{point_cloud_info['bounds']['z_min']:.3f}, {point_cloud_info['bounds']['z_max']:.3f}]\n\n")
        
        # Grid information
        f.write("Grid Information:\n")
        f.write(f"  Grid shape: {grid_array.shape}\n")
        f.write(f"  Resolution: {0.5}\n")  # Hardcoded resolution
        f.write(f"  Value range: {grid_array.min():.2f} to {grid_array.max():.2f}\n")
        f.write(f"  Non-zero voxels: {np.sum(grid_array > 0)}\n")
        f.write(f"  Total voxels: {grid_array.size}\n")
        f.write(f"  Sparsity: {1.0 - np.sum(grid_array > 0) / grid_array.size:.4f}\n\n")
        
        # Curve information
        f.write("Extracted Curves:\n")
        f.write(f"  Number of curves: {len(extracted_curves)}\n")
        
        total_points = 0
        total_length = 0.0
        
        for i, curve in enumerate(extracted_curves):
            f.write(f"\n  Curve {i+1}:\n")
            f.write(f"    Number of points: {curve.get_n_points()}\n")
            f.write(f"    3D length: {curve.length_3d():.2f}\n")
            
            if len(curve.points) > 0:
                points = curve.get_points_as_array()
                f.write(f"    Bounding box: X[{points[:, 0].min():.1f}, {points[:, 0].max():.1f}], "
                       f"Y[{points[:, 1].min():.1f}, {points[:, 1].max():.1f}], "
                       f"Z[{points[:, 2].min():.1f}, {points[:, 2].max():.1f}]\n")
                
                total_points += len(curve.points)
                total_length += curve.length_3d()
        
        f.write(f"\nSummary:\n")
        f.write(f"  Total points across all curves: {total_points}\n")
        f.write(f"  Total length across all curves: {total_length:.2f}\n")
    
    print(f"Analysis report saved to {filename}")


def main():
    """Main function to demonstrate the active contours workflow with twoCylinders.txt"""
    print("=== Active Contours Workflow Example with Two Cylinders ===\n")
    
    # Load data from twoCylinders.txt
    grid_data, resolution, point_cloud_info = load_two_cylinders_data()
    
    # Convert grid to numpy array for the workflow
    grid_array = grid_data.get_data()
    
    # Create workflow
    print("\nInitializing Active Contours workflow...")
    workflow = ActiveContoursWorkflow(grid_array, resolution)
    
    # Configure workflow parameters (adjusted for real data)
    workflow.n_snakes_max = 10  # Allow more curves for real data
    workflow.longueur_min = 1.0  # Minimum curve length
    workflow.min_value_for_maxima = 5  # Lower threshold for real data
    workflow.maxima_neighborhood_size = 3  # Larger neighborhood for real data
    
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
    
    # Create output directory
    output_dir = "output_two_cylinders"
    os.makedirs(output_dir, exist_ok=True)
    
    # Visualize results
    visualize_grid_and_curves(
        grid_array, extracted_curves, point_cloud_info,
        save_path=os.path.join(output_dir, "two_cylinders_active_contours_visualization.png")
    )
    
    # Save curves to file
    save_curves_to_file(extracted_curves, os.path.join(output_dir, "two_cylinders_extracted_curves.txt"))
    
    # Save comprehensive analysis report
    save_analysis_report(
        grid_array, extracted_curves, point_cloud_info,
        os.path.join(output_dir, "two_cylinders_analysis_report.txt")
    )
    
    print(f"\nAll results saved to {output_dir}/")
    print("\n=== Example completed successfully! ===")


if __name__ == "__main__":
    main() 