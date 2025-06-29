#!/usr/bin/env python3
"""
3D Dataset Visualization Example using twoCylinders.txt

This example demonstrates various ways to visualize the 3D point cloud data:
1. 3D scatter plot of points
2. 3D scatter plot with normal vectors
3. Density visualization using 2D projections
4. Interactive 3D visualization (if available)
5. Statistical analysis with visualizations

Requirements:
- numpy
- matplotlib
- scipy
- mplot3d (for 3D plots)
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.colors as mcolors

# Add the parent directory to the path to import PythonStep modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from PythonStep.utils import read_point_cloud_file, validate_point_cloud, normalize_normals


def create_3d_scatter_plot(points, normals=None, title="3D Point Cloud Visualization"):
    """
    Create a 3D scatter plot of the point cloud
    """
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Create scatter plot
    scatter = ax.scatter(points[:, 0], points[:, 1], points[:, 2], 
                        c=points[:, 2], cmap='viridis', s=10, alpha=0.7)
    
    # Add normal vectors if provided
    if normals is not None:
        # Sample every 50th normal to avoid clutter
        sample_indices = np.arange(0, len(normals), 50)
        sample_points = points[sample_indices]
        sample_normals = normals[sample_indices]
        
        # Scale normals for visualization
        scale = 0.1
        normal_ends = sample_points + sample_normals * scale
        
        # Create line segments for normal vectors
        lines = []
        for start, end in zip(sample_points, normal_ends):
            lines.append([start, end])
        
        line_collection = Line3DCollection(lines, colors='red', alpha=0.5, linewidth=1)
        ax.add_collection3d(line_collection)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    
    # Add colorbar
    plt.colorbar(scatter, ax=ax, shrink=0.5, aspect=20, label='Z coordinate')
    
    return fig, ax


def create_density_projections(points, title="Density Projections"):
    """
    Create 2D density projections of the 3D data
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(title, fontsize=16)
    
    # XY projection
    axes[0, 0].hist2d(points[:, 0], points[:, 1], bins=50, cmap='viridis')
    axes[0, 0].set_xlabel('X')
    axes[0, 0].set_ylabel('Y')
    axes[0, 0].set_title('XY Projection (Density)')
    axes[0, 0].set_aspect('equal')
    
    # XZ projection
    axes[0, 1].hist2d(points[:, 0], points[:, 2], bins=50, cmap='viridis')
    axes[0, 1].set_xlabel('X')
    axes[0, 1].set_ylabel('Z')
    axes[0, 1].set_title('XZ Projection (Density)')
    axes[0, 1].set_aspect('equal')
    
    # YZ projection
    axes[1, 0].hist2d(points[:, 1], points[:, 2], bins=50, cmap='viridis')
    axes[1, 0].set_xlabel('Y')
    axes[1, 0].set_ylabel('Z')
    axes[1, 0].set_title('YZ Projection (Density)')
    axes[1, 0].set_aspect('equal')
    
    # 3D scatter plot (top view)
    ax_3d = fig.add_subplot(2, 2, 4, projection='3d')
    scatter = ax_3d.scatter(points[:, 0], points[:, 1], points[:, 2], 
                           c=points[:, 2], cmap='viridis', s=5, alpha=0.6)
    ax_3d.set_xlabel('X')
    ax_3d.set_ylabel('Y')
    ax_3d.set_zlabel('Z')
    ax_3d.set_title('3D Scatter (Top View)')
    
    plt.tight_layout()
    return fig


def create_statistical_visualizations(points, normals):
    """
    Create statistical visualizations of the dataset
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("Statistical Analysis of 3D Dataset", fontsize=16)
    
    # Coordinate distributions
    axes[0, 0].hist(points[:, 0], bins=50, alpha=0.7, label='X', color='red')
    axes[0, 0].hist(points[:, 1], bins=50, alpha=0.7, label='Y', color='green')
    axes[0, 0].hist(points[:, 2], bins=50, alpha=0.7, label='Z', color='blue')
    axes[0, 0].set_xlabel('Coordinate Value')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Coordinate Distributions')
    axes[0, 0].legend()
    
    # Normal vector distributions
    if normals is not None:
        axes[0, 1].hist(normals[:, 0], bins=50, alpha=0.7, label='Nx', color='red')
        axes[0, 1].hist(normals[:, 1], bins=50, alpha=0.7, label='Ny', color='green')
        axes[0, 1].hist(normals[:, 2], bins=50, alpha=0.7, label='Nz', color='blue')
        axes[0, 1].set_xlabel('Normal Component')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('Normal Vector Distributions')
        axes[0, 1].legend()
    
    # Distance from origin
    distances = np.linalg.norm(points, axis=1)
    axes[0, 2].hist(distances, bins=50, alpha=0.7, color='purple')
    axes[0, 2].set_xlabel('Distance from Origin')
    axes[0, 2].set_ylabel('Frequency')
    axes[0, 2].set_title('Distance Distribution')
    
    # Scatter plots of coordinate pairs
    axes[1, 0].scatter(points[:, 0], points[:, 1], s=1, alpha=0.5)
    axes[1, 0].set_xlabel('X')
    axes[1, 0].set_ylabel('Y')
    axes[1, 0].set_title('X vs Y')
    axes[1, 0].set_aspect('equal')
    
    axes[1, 1].scatter(points[:, 0], points[:, 2], s=1, alpha=0.5)
    axes[1, 1].set_xlabel('X')
    axes[1, 1].set_ylabel('Z')
    axes[1, 1].set_title('X vs Z')
    axes[1, 1].set_aspect('equal')
    
    axes[1, 2].scatter(points[:, 1], points[:, 2], s=1, alpha=0.5)
    axes[1, 2].set_xlabel('Y')
    axes[1, 2].set_ylabel('Z')
    axes[1, 2].set_title('Y vs Z')
    axes[1, 2].set_aspect('equal')
    
    plt.tight_layout()
    return fig


def create_cluster_visualization(points, n_clusters=2):
    """
    Create visualization showing clustering of the points
    """
    from scipy.cluster.vq import kmeans, vq
    
    # Perform k-means clustering
    centroids, labels = kmeans(points, n_clusters)
    
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Create scatter plot with colors based on cluster labels
    colors = plt.cm.Set1(np.linspace(0, 1, n_clusters))
    
    for i in range(n_clusters):
        cluster_points = points[labels == i]
        ax.scatter(cluster_points[:, 0], cluster_points[:, 1], cluster_points[:, 2],
                  c=[colors[i]], s=10, alpha=0.7, label=f'Cluster {i+1}')
    
    # Plot centroids
    ax.scatter(centroids[:, 0], centroids[:, 1], centroids[:, 2],
              c='black', s=200, marker='*', label='Centroids')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'K-means Clustering (k={n_clusters})')
    ax.legend()
    
    return fig, centroids, labels


def create_interactive_3d_plot(points, normals=None):
    """
    Create an interactive 3D plot (if matplotlib backend supports it)
    """
    try:
        # Try to use an interactive backend
        plt.ion()
        
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Create scatter plot
        scatter = ax.scatter(points[:, 0], points[:, 1], points[:, 2], 
                           c=points[:, 2], cmap='viridis', s=10, alpha=0.7)
        
        # Add normal vectors if provided
        if normals is not None:
            # Sample every 100th normal to avoid clutter
            sample_indices = np.arange(0, len(normals), 100)
            sample_points = points[sample_indices]
            sample_normals = normals[sample_indices]
            
            # Scale normals for visualization
            scale = 0.1
            normal_ends = sample_points + sample_normals * scale
            
            # Create line segments for normal vectors
            for start, end in zip(sample_points, normal_ends):
                ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], 
                       'r-', alpha=0.5, linewidth=1)
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Interactive 3D Point Cloud (Rotate with mouse)')
        
        # Add colorbar
        plt.colorbar(scatter, ax=ax, shrink=0.5, aspect=20, label='Z coordinate')
        
        print("Interactive plot created. You can rotate it with your mouse.")
        print("Close the plot window to continue...")
        
        plt.show(block=True)
        
        return fig, ax
        
    except Exception as e:
        print(f"Interactive plotting not available: {e}")
        return None, None


def main():
    print("=== 3D Dataset Visualization Example ===\n")
    
    # File path
    data_file = "/home/sfoucher/DEV/pluginsteplight/twoCylinders.txt"
    
    if not os.path.exists(data_file):
        print(f"Error: File not found: {data_file}")
        return
    
    # Load the data
    print("1. Loading 3D dataset...")
    try:
        points, normals = read_point_cloud_file(data_file)
        print(f"   Loaded {len(points)} points")
        print(f"   Dataset bounds: {points.min(axis=0)} to {points.max(axis=0)}")
        
        # Validate and normalize normals
        validate_point_cloud(points, normals)
        normals = normalize_normals(normals)
        print("   ✓ Data validation passed")
        
    except Exception as e:
        print(f"   Error loading data: {e}")
        return
    
    # Create visualizations
    print("\n2. Creating visualizations...")
    
    # 2.1 Basic 3D scatter plot
    print("   2.1 Creating 3D scatter plot...")
    fig1, ax1 = create_3d_scatter_plot(points, normals, "Two Cylinders Dataset - 3D View")
    plt.savefig("3d_scatter_plot.png", dpi=150, bbox_inches='tight')
    print("   ✓ Saved: 3d_scatter_plot.png")
    
    # 2.2 Density projections
    print("   2.2 Creating density projections...")
    fig2 = create_density_projections(points, "Two Cylinders Dataset - Density Projections")
    plt.savefig("density_projections.png", dpi=150, bbox_inches='tight')
    print("   ✓ Saved: density_projections.png")
    
    # 2.3 Statistical visualizations
    print("   2.3 Creating statistical visualizations...")
    fig3 = create_statistical_visualizations(points, normals)
    plt.savefig("statistical_analysis.png", dpi=150, bbox_inches='tight')
    print("   ✓ Saved: statistical_analysis.png")
    
    # 2.4 Cluster visualization
    print("   2.4 Creating cluster visualization...")
    fig4, centroids, labels = create_cluster_visualization(points, n_clusters=2)
    plt.savefig("cluster_visualization.png", dpi=150, bbox_inches='tight')
    print("   ✓ Saved: cluster_visualization.png")
    
    # 2.5 Interactive 3D plot
    print("   2.5 Creating interactive 3D plot...")
    fig5, ax5 = create_interactive_3d_plot(points, normals)
    if fig5 is not None:
        plt.savefig("interactive_3d_plot.png", dpi=150, bbox_inches='tight')
        print("   ✓ Saved: interactive_3d_plot.png")
    
    # Print dataset statistics
    print("\n3. Dataset Statistics:")
    print(f"   - Number of points: {len(points)}")
    print(f"   - X range: {points[:, 0].min():.3f} to {points[:, 0].max():.3f}")
    print(f"   - Y range: {points[:, 1].min():.3f} to {points[:, 1].max():.3f}")
    print(f"   - Z range: {points[:, 2].min():.3f} to {points[:, 2].max():.3f}")
    print(f"   - Mean position: {points.mean(axis=0)}")
    print(f"   - Standard deviation: {points.std(axis=0)}")
    
    if normals is not None:
        print(f"   - Normal vector mean: {normals.mean(axis=0)}")
        print(f"   - Normal vector std: {normals.std(axis=0)}")
    
    # Cluster information
    from scipy.cluster.vq import kmeans, vq
    centroids, labels = kmeans(points, 2)
    print(f"   - Cluster 1 centroid: {centroids[0]}")
    print(f"   - Cluster 2 centroid: {centroids[1]}")
    print(f"   - Points in cluster 1: {np.sum(labels == 0)}")
    print(f"   - Points in cluster 2: {np.sum(labels == 1)}")
    
    # Save dataset summary
    print("\n4. Saving dataset summary...")
    with open("dataset_summary.txt", "w") as f:
        f.write("Two Cylinders Dataset Summary\n")
        f.write("=" * 30 + "\n\n")
        f.write(f"File: {data_file}\n")
        f.write(f"Number of points: {len(points)}\n")
        f.write(f"Dataset bounds:\n")
        f.write(f"  X: {points[:, 0].min():.6f} to {points[:, 0].max():.6f}\n")
        f.write(f"  Y: {points[:, 1].min():.6f} to {points[:, 1].max():.6f}\n")
        f.write(f"  Z: {points[:, 2].min():.6f} to {points[:, 2].max():.6f}\n")
        f.write(f"Mean position: {points.mean(axis=0)}\n")
        f.write(f"Standard deviation: {points.std(axis=0)}\n")
        if normals is not None:
            f.write(f"Normal vector mean: {normals.mean(axis=0)}\n")
            f.write(f"Normal vector std: {normals.std(axis=0)}\n")
        f.write(f"Cluster centroids:\n")
        f.write(f"  Cluster 1: {centroids[0]}\n")
        f.write(f"  Cluster 2: {centroids[1]}\n")
    
    print("   ✓ Saved: dataset_summary.txt")
    
    print("\n=== Visualization completed successfully! ===")
    print("Files created:")
    print("  - 3d_scatter_plot.png: 3D scatter plot with normal vectors")
    print("  - density_projections.png: 2D density projections")
    print("  - statistical_analysis.png: Statistical distributions")
    print("  - cluster_visualization.png: K-means clustering results")
    print("  - interactive_3d_plot.png: Interactive 3D view")
    print("  - dataset_summary.txt: Text summary of the dataset")


if __name__ == "__main__":
    main() 