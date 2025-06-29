"""
PythonStep - Python wrapper for STL_Grid3D using Cython

This package provides Python bindings for the STL_Grid3D class from the pluginsteplight project.
"""

from .stl_grid3d_wrapper import PySTLGrid3D, create_grid_from_txt
from .utils import (
    read_point_cloud_file,
    calculate_bounding_box,
    calculate_grid_dimensions,
    validate_point_cloud,
    normalize_normals,
    sample_points_on_grid,
    create_grid_from_points
)

__version__ = "0.1.0"
__author__ = "Your Name"

__all__ = [
    "PySTLGrid3D",
    "create_grid_from_txt", 
    "read_point_cloud_file",
    "calculate_bounding_box",
    "calculate_grid_dimensions", 
    "validate_point_cloud",
    "normalize_normals",
    "sample_points_on_grid",
    "create_grid_from_points",
] 