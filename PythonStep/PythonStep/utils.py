"""
Utility functions for PythonStep package
"""

import numpy as np
import os
from typing import Tuple, Optional


def read_point_cloud_file(filename: str, skip_rows: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read point cloud data from a text file in twoCylinders.txt format
    
    Args:
        filename: Path to the text file
        skip_rows: Number of header rows to skip (default: 1 for comment line)
    
    Returns:
        Tuple of (points, normals) where:
        - points: numpy array of shape (n_points, 3) with X, Y, Z coordinates
        - normals: numpy array of shape (n_points, 3) with Nx, Ny, Nz normal vectors
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    
    try:
        data = np.loadtxt(filename, skiprows=skip_rows)
    except Exception as e:
        raise ValueError(f"Error reading file {filename}: {e}")
    
    if data.shape[1] != 6:
        raise ValueError(f"File must have 6 columns (X Y Z Nx Ny Nz), got {data.shape[1]}")
    
    points = data[:, :3].astype(np.float64)
    normals = data[:, 3:].astype(np.float64)
    
    return points, normals


def calculate_bounding_box(points: np.ndarray, padding: float = 0.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate bounding box for a set of points
    
    Args:
        points: numpy array of shape (n_points, 3) with X, Y, Z coordinates
        padding: Additional padding to add to the bounding box
    
    Returns:
        Tuple of (min_coords, max_coords) where each is a numpy array of shape (3,)
    """
    min_coords = points.min(axis=0) - padding
    max_coords = points.max(axis=0) + padding
    
    return min_coords, max_coords


def calculate_grid_dimensions(bounds: Tuple[np.ndarray, np.ndarray], 
                            resolution: float) -> Tuple[int, int, int]:
    """
    Calculate grid dimensions from bounding box and resolution
    
    Args:
        bounds: Tuple of (min_coords, max_coords)
        resolution: Grid resolution
    
    Returns:
        Tuple of (nx, ny, nz) grid dimensions
    """
    min_coords, max_coords = bounds
    dimensions = max_coords - min_coords
    
    nx = int(np.ceil(dimensions[0] / resolution))
    ny = int(np.ceil(dimensions[1] / resolution))
    nz = int(np.ceil(dimensions[2] / resolution))
    
    return nx, ny, nz


def validate_point_cloud(points: np.ndarray, normals: Optional[np.ndarray] = None) -> None:
    """
    Validate point cloud data
    
    Args:
        points: numpy array of shape (n_points, 3) with X, Y, Z coordinates
        normals: optional numpy array of shape (n_points, 3) with normal vectors
    
    Raises:
        ValueError: If data is invalid
    """
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("Points must be a 2D array with shape (n_points, 3)")
    
    if points.shape[0] == 0:
        raise ValueError("Points array cannot be empty")
    
    if normals is not None:
        if normals.shape != points.shape:
            raise ValueError("Normals must have the same shape as points")
        
        # Check if normals are unit vectors (approximately)
        norms = np.linalg.norm(normals, axis=1)
        if not np.allclose(norms, 1.0, atol=1e-6):
            print("Warning: Normals are not unit vectors. Consider normalizing.")


def normalize_normals(normals: np.ndarray) -> np.ndarray:
    """
    Normalize normal vectors to unit length
    
    Args:
        normals: numpy array of shape (n_points, 3) with normal vectors
    
    Returns:
        Normalized normals array
    """
    norms = np.linalg.norm(normals, axis=1, keepdims=True)
    # Avoid division by zero
    norms = np.where(norms == 0, 1.0, norms)
    return normals / norms


def sample_points_on_grid(points: np.ndarray, resolution: float, 
                         bounds: Optional[Tuple[np.ndarray, np.ndarray]] = None) -> np.ndarray:
    """
    Sample points to grid coordinates
    
    Args:
        points: numpy array of shape (n_points, 3) with X, Y, Z coordinates
        resolution: Grid resolution
        bounds: Optional bounding box, if None will be calculated from points
    
    Returns:
        numpy array of shape (n_points, 3) with grid coordinates
    """
    if bounds is None:
        bounds = calculate_bounding_box(points)
    
    min_coords, _ = bounds
    
    # Convert to grid coordinates
    grid_coords = np.floor((points - min_coords) / resolution).astype(int)
    
    return grid_coords


def create_grid_from_points(points: np.ndarray, normals: Optional[np.ndarray] = None,
                           resolution: float = 0.1, padding: float = 0.0,
                           na_value: int = -9999) -> 'PySTLGrid3D':
    """
    Create STL_Grid3D from point cloud data
    
    Args:
        points: numpy array of shape (n_points, 3) with X, Y, Z coordinates
        normals: optional numpy array of shape (n_points, 3) with normal vectors
        resolution: Grid resolution
        padding: Additional padding around the bounding box
        na_value: Value to use for NA/missing data
    
    Returns:
        PySTLGrid3D: Initialized grid object
    """
    from .stl_grid3d_wrapper import PySTLGrid3D
    
    # Validate input
    validate_point_cloud(points, normals)
    
    # Calculate bounding box
    bounds = calculate_bounding_box(points, padding)
    min_coords, max_coords = bounds
    
    # Create grid
    grid = PySTLGrid3D(resolution)
    grid.create_from_bounds(
        min_coords[0], min_coords[1], min_coords[2],
        max_coords[0], max_coords[1], max_coords[2],
        resolution, na_value, 0
    )
    
    # TODO: Set point cloud and normals if provided
    # This would require implementing the point cloud wrapper
    
    return grid 