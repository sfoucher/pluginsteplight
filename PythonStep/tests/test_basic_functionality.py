#!/usr/bin/env python3
"""
Basic unit tests for PythonStep package
"""

import unittest
import numpy as np
import os
import sys
import tempfile

# Add the parent directory to the path to import PythonStep
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

try:
    import PythonStep as ps
except ImportError as e:
    print(f"Error importing PythonStep: {e}")
    print("Make sure the package is properly installed and all dependencies are available.")
    sys.exit(1)


class TestPythonStep(unittest.TestCase):
    """Test cases for PythonStep package"""
    
    def setUp(self):
        """Set up test fixtures"""
        # Create a simple test point cloud
        self.test_points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ], dtype=np.float64)
        
        self.test_normals = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
        ], dtype=np.float64)
        
        # Create a temporary test file
        self.temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
        self.temp_file.write("//X Y Z Nx Ny Nz\n")
        for i in range(len(self.test_points)):
            point = self.test_points[i]
            normal = self.test_normals[i]
            self.temp_file.write(f"{point[0]} {point[1]} {point[2]} {normal[0]} {normal[1]} {normal[2]}\n")
        self.temp_file.close()
    
    def tearDown(self):
        """Clean up test fixtures"""
        if hasattr(self, 'temp_file') and os.path.exists(self.temp_file.name):
            os.unlink(self.temp_file.name)
    
    def test_read_point_cloud_file(self):
        """Test reading point cloud from file"""
        points, normals = ps.read_point_cloud_file(self.temp_file.name)
        
        self.assertEqual(points.shape, (5, 3))
        self.assertEqual(normals.shape, (5, 3))
        np.testing.assert_array_almost_equal(points, self.test_points)
        np.testing.assert_array_almost_equal(normals, self.test_normals)
    
    def test_calculate_bounding_box(self):
        """Test bounding box calculation"""
        min_coords, max_coords = ps.calculate_bounding_box(self.test_points)
        
        expected_min = np.array([0.0, 0.0, 0.0])
        expected_max = np.array([1.0, 1.0, 1.0])
        
        np.testing.assert_array_almost_equal(min_coords, expected_min)
        np.testing.assert_array_almost_equal(max_coords, expected_max)
    
    def test_calculate_bounding_box_with_padding(self):
        """Test bounding box calculation with padding"""
        min_coords, max_coords = ps.calculate_bounding_box(self.test_points, padding=0.5)
        
        expected_min = np.array([-0.5, -0.5, -0.5])
        expected_max = np.array([1.5, 1.5, 1.5])
        
        np.testing.assert_array_almost_equal(min_coords, expected_min)
        np.testing.assert_array_almost_equal(max_coords, expected_max)
    
    def test_calculate_grid_dimensions(self):
        """Test grid dimension calculation"""
        bounds = ps.calculate_bounding_box(self.test_points)
        nx, ny, nz = ps.calculate_grid_dimensions(bounds, resolution=0.5)
        
        # With resolution 0.5, we expect 2x2x2 grid
        self.assertEqual(nx, 2)
        self.assertEqual(ny, 2)
        self.assertEqual(nz, 2)
    
    def test_validate_point_cloud(self):
        """Test point cloud validation"""
        # Valid point cloud should not raise exception
        ps.validate_point_cloud(self.test_points, self.test_normals)
        
        # Invalid shapes should raise exception
        with self.assertRaises(ValueError):
            ps.validate_point_cloud(self.test_points, self.test_normals[:3])
        
        # Empty point cloud should raise exception
        with self.assertRaises(ValueError):
            ps.validate_point_cloud(np.empty((0, 3)))
    
    def test_normalize_normals(self):
        """Test normal vector normalization"""
        # Create non-normalized normals
        non_normalized = np.array([
            [2.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 4.0],
        ])
        
        normalized = ps.normalize_normals(non_normalized)
        
        # Check that norms are approximately 1
        norms = np.linalg.norm(normalized, axis=1)
        np.testing.assert_array_almost_equal(norms, np.ones(3))
    
    def test_sample_points_on_grid(self):
        """Test point sampling to grid coordinates"""
        bounds = ps.calculate_bounding_box(self.test_points)
        grid_coords = ps.sample_points_on_grid(self.test_points, resolution=0.5, bounds=bounds)
        
        # Check that grid coordinates are integers
        self.assertTrue(np.all(grid_coords == grid_coords.astype(int)))
        
        # Check that coordinates are within expected range
        self.assertTrue(np.all(grid_coords >= 0))
        self.assertTrue(np.all(grid_coords < 2))  # 2x2x2 grid
    
    def test_pystl_grid3d_creation(self):
        """Test PySTLGrid3D creation"""
        grid = ps.PySTLGrid3D(resolution=0.1)
        
        # Test creation from bounds
        grid.create_from_bounds(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.1, -9999, 0)
        
        self.assertEqual(grid.dimensions, (10, 10, 10))  # 1.0 / 0.1 = 10
        self.assertEqual(grid.resolution, 0.1)
        self.assertEqual(grid.bounds, (0.0, 0.0, 0.0, 1.0, 1.0, 1.0))
    
    def test_pystl_grid3d_data_access(self):
        """Test PySTLGrid3D data access"""
        grid = ps.PySTLGrid3D(resolution=0.5)
        grid.create_from_bounds(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.5, -9999, 0)
        
        # Test getting data
        data = grid.get_data()
        self.assertEqual(data.shape, (2, 2, 2))
        self.assertEqual(data.dtype, np.int32)
        
        # Test getting values
        value = grid.get_value(0, 0, 0)
        self.assertEqual(value, 0)  # Initial value
    
    def test_pystl_grid3d_coordinate_conversion(self):
        """Test PySTLGrid3D coordinate conversion"""
        grid = ps.PySTLGrid3D(resolution=0.5)
        grid.create_from_bounds(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.5, -9999, 0)
        
        # Test pixel to cartesian
        cartesian = grid.pixel_to_cartesian(0, 0, 0)
        expected = (0.25, 0.25, 0.25)  # Center of first voxel
        np.testing.assert_array_almost_equal(cartesian, expected)
        
        # Test cartesian to pixel
        pixel = grid.cartesian_to_pixel(0.25, 0.25, 0.25)
        self.assertEqual(pixel, (0, 0, 0))
    
    def test_create_grid_from_txt(self):
        """Test creating grid from text file"""
        grid = ps.create_grid_from_txt(self.temp_file.name, resolution=0.5)
        
        self.assertIsInstance(grid, ps.PySTLGrid3D)
        self.assertEqual(grid.dimensions, (2, 2, 2))  # 1.0 / 0.5 = 2
    
    def test_create_grid_from_points(self):
        """Test creating grid from points"""
        grid = ps.create_grid_from_points(
            points=self.test_points,
            normals=self.test_normals,
            resolution=0.5,
            padding=0.1
        )
        
        self.assertIsInstance(grid, ps.PySTLGrid3D)
        # With padding 0.1, bounds should be [-0.1, 1.1] in each dimension
        # With resolution 0.5, dimensions should be 3x3x3
        self.assertEqual(grid.dimensions, (3, 3, 3))


if __name__ == '__main__':
    unittest.main() 