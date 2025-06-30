"""
Tests for the active contours workflow
"""

import numpy as np
import pytest
from PythonStep.workflows import OpenActiveContours, ActiveContoursWorkflow, ContourPoint


class TestContourPoint:
    """Test ContourPoint class"""
    
    def test_contour_point_creation(self):
        """Test creating a contour point"""
        point = ContourPoint(1.0, 2.0, 3.0, 5.0)
        assert point.x == 1.0
        assert point.y == 2.0
        assert point.z == 3.0
        assert point.value == 5.0
    
    def test_contour_point_to_array(self):
        """Test converting contour point to array"""
        point = ContourPoint(1.0, 2.0, 3.0)
        arr = point.to_array()
        assert np.array_equal(arr, np.array([1.0, 2.0, 3.0]))
    
    def test_contour_point_from_array(self):
        """Test creating contour point from array"""
        arr = np.array([1.0, 2.0, 3.0])
        point = ContourPoint.from_array(arr, 5.0)
        assert point.x == 1.0
        assert point.y == 2.0
        assert point.z == 3.0
        assert point.value == 5.0


class TestOpenActiveContours:
    """Test OpenActiveContours class"""
    
    def setup_method(self):
        """Set up test data"""
        # Create a simple 3D grid with some high values
        self.grid_data = np.zeros((10, 10, 10))
        
        # Add some high values in a line pattern (simulating a curve)
        for i in range(3, 7):
            self.grid_data[i, i, i] = 50
            self.grid_data[i, i, i+1] = 40
            self.grid_data[i+1, i, i] = 30
        
        self.resolution = 1.0
        self.contour = OpenActiveContours(
            self.grid_data, self.resolution, min_value=10.0, search_cone_size=3
        )
    
    def test_initialization(self):
        """Test contour initialization"""
        assert self.contour.grid_data.shape == (10, 10, 10)
        assert self.contour.resolution == 1.0
        assert self.contour.min_value == 10.0
        assert len(self.contour.points) == 0
    
    def test_initialize_from_point_success(self):
        """Test successful initialization from a valid point"""
        start_pixel = (4, 4, 4)  # Point with high value
        success = self.contour.initialize_from_point(start_pixel)
        assert success
        assert len(self.contour.points) > 0
    
    def test_initialize_from_point_failure(self):
        """Test failed initialization from invalid point"""
        start_pixel = (0, 0, 0)  # Point with zero value
        success = self.contour.initialize_from_point(start_pixel)
        assert not success
        assert len(self.contour.points) == 0
    
    def test_initialize_from_point_out_of_bounds(self):
        """Test initialization from out-of-bounds point"""
        start_pixel = (15, 15, 15)  # Out of bounds
        success = self.contour.initialize_from_point(start_pixel)
        assert not success
    
    def test_length_3d_empty(self):
        """Test 3D length calculation for empty contour"""
        length = self.contour.length_3d()
        assert length == 0.0
    
    def test_length_3d_with_points(self):
        """Test 3D length calculation with points"""
        # Add some points manually
        self.contour.points = [
            ContourPoint(0.0, 0.0, 0.0),
            ContourPoint(1.0, 0.0, 0.0),
            ContourPoint(1.0, 1.0, 0.0)
        ]
        length = self.contour.length_3d()
        assert length > 0.0
        assert abs(length - 2.0) < 0.1  # Should be approximately 2.0
    
    def test_resample(self):
        """Test contour resampling"""
        # Add some points manually
        self.contour.points = [
            ContourPoint(0.0, 0.0, 0.0),
            ContourPoint(1.0, 0.0, 0.0),
            ContourPoint(2.0, 0.0, 0.0)
        ]
        
        original_length = self.contour.length_3d()
        original_points = len(self.contour.points)
        
        # Resample with smaller resolution
        self.contour.resample(0.5)
        
        # Should have more points after resampling
        assert len(self.contour.points) > original_points
        # Length should be approximately the same
        assert abs(self.contour.length_3d() - original_length) < 0.1
    
    def test_get_tangent_at_point(self):
        """Test tangent calculation"""
        # Add some points manually
        self.contour.points = [
            ContourPoint(0.0, 0.0, 0.0),
            ContourPoint(1.0, 0.0, 0.0),
            ContourPoint(2.0, 0.0, 0.0)
        ]
        
        # Test tangent at middle point
        tangent = self.contour._get_tangent_at_point(1, normalize=True)
        assert np.linalg.norm(tangent) > 0.0
        
        # Test tangent at head
        head_tangent = self.contour._get_tangent_at_head(True)
        assert np.linalg.norm(head_tangent) > 0.0
        
        # Test tangent at tail
        tail_tangent = self.contour._get_tangent_at_tail(True)
        assert np.linalg.norm(tail_tangent) > 0.0
    
    def test_coordinate_conversion(self):
        """Test coordinate conversion methods"""
        # Test pixel to cartesian
        cartesian = self.contour._pixel_to_cartesian(2, 3, 4)
        assert np.array_equal(cartesian, np.array([2.0, 3.0, 4.0]))
        
        # Test cartesian to pixel
        pixel = self.contour._cartesian_to_pixel(np.array([2.5, 3.7, 4.2]))
        assert pixel == (2, 3, 4)
    
    def test_valid_pixel_check(self):
        """Test pixel validity checking"""
        # Valid pixels
        assert self.contour._is_valid_pixel(0, 0, 0)
        assert self.contour._is_valid_pixel(9, 9, 9)
        
        # Invalid pixels
        assert not self.contour._is_valid_pixel(-1, 0, 0)
        assert not self.contour._is_valid_pixel(10, 0, 0)
        assert not self.contour._is_valid_pixel(0, -1, 0)
        assert not self.contour._is_valid_pixel(0, 10, 0)
        assert not self.contour._is_valid_pixel(0, 0, -1)
        assert not self.contour._is_valid_pixel(0, 0, 10)


class TestActiveContoursWorkflow:
    """Test ActiveContoursWorkflow class"""
    
    def setup_method(self):
        """Set up test data"""
        # Create a simple 3D grid with some high values
        self.grid_data = np.zeros((10, 10, 10))
        
        # Add some high values in a line pattern (simulating a curve)
        for i in range(3, 7):
            self.grid_data[i, i, i] = 50
            self.grid_data[i, i, i+1] = 40
            self.grid_data[i+1, i, i] = 30
        
        self.resolution = 1.0
        self.workflow = ActiveContoursWorkflow(self.grid_data, self.resolution)
    
    def test_workflow_initialization(self):
        """Test workflow initialization"""
        assert self.workflow.grid_data.shape == (10, 10, 10)
        assert self.workflow.resolution == 1.0
        assert self.workflow.n_snakes_max == 1
        assert self.workflow.longueur_min == 1.0
    
    def test_get_local_maximas(self):
        """Test local maxima detection"""
        local_maximas = self.workflow._get_local_maximas()
        assert len(local_maximas) > 0
        
        # Check that all maxima are within height constraints
        for x, y, z in local_maximas:
            height = z * self.resolution
            assert (self.workflow.min_height_for_maximum_search <= height <= 
                   self.workflow.max_height_for_maximum_search)
    
    def test_extract_curves(self):
        """Test curve extraction"""
        curves = self.workflow.extract_curves()
        
        # Should return a list
        assert isinstance(curves, list)
        
        # If curves are found, they should be OpenActiveContours objects
        for curve in curves:
            assert isinstance(curve, OpenActiveContours)
            assert curve.length_3d() > self.workflow.longueur_min
    
    def test_mark_repulsion(self):
        """Test repulsion marking"""
        # Create a simple contour
        contour = OpenActiveContours(self.grid_data, self.resolution)
        contour.points = [
            ContourPoint(1.0, 1.0, 1.0),
            ContourPoint(2.0, 1.0, 1.0)
        ]
        
        repulsion_grid = np.zeros_like(self.grid_data, dtype=bool)
        
        # Mark repulsion
        self.workflow._mark_repulsion(contour, repulsion_grid, 1.5)
        
        # Check that some areas are marked as repulsive
        assert np.any(repulsion_grid)


def test_integration():
    """Integration test for the complete workflow"""
    # Create a more complex test grid
    grid_data = np.zeros((20, 20, 20))
    
    # Create a curved path with high values
    for t in np.linspace(0, 2*np.pi, 20):
        x = int(10 + 5*np.cos(t))
        y = int(10 + 5*np.sin(t))
        z = int(10 + 2*np.sin(2*t))
        
        # Add high values along the curve
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                for dz in range(-1, 2):
                    nx, ny, nz = x + dx, y + dy, z + dz
                    if (0 <= nx < 20 and 0 <= ny < 20 and 0 <= nz < 20):
                        grid_data[nx, ny, nz] = 30 + 20*np.exp(-(dx**2 + dy**2 + dz**2)/2)
    
    # Create workflow
    workflow = ActiveContoursWorkflow(grid_data, resolution=1.0)
    workflow.n_snakes_max = 2  # Allow multiple curves
    workflow.longueur_min = 0.5  # Lower minimum length for testing
    
    # Extract curves
    curves = workflow.extract_curves()
    
    # Should find at least one curve
    assert len(curves) > 0
    
    # Each curve should have reasonable properties
    for curve in curves:
        assert curve.length_3d() > workflow.longueur_min
        assert curve.get_n_points() > 2
        assert len(curve.get_points_as_array()) > 0


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"]) 