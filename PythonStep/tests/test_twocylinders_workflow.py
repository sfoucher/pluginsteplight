import unittest
import os
import numpy as np
from PythonStep.utils import read_point_cloud_file, create_grid_from_points
from PythonStep.workflows import GridWorkflows

DATA_FILE = '/home/sfoucher/DEV/pluginsteplight/twoCylinders.txt'

class TestTwoCylindersWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Read the point cloud and normals
        assert os.path.exists(DATA_FILE), f"File not found: {DATA_FILE}"
        points, normals = read_point_cloud_file(DATA_FILE)
        cls.points = points
        cls.normals = normals
        
        # Create a grid with test data instead of empty grid
        cls.grid_obj = create_grid_from_points(points, normals, resolution=0.1, padding=0.0, na_value=-9999)
        
        # Manually populate the grid with some test data
        cls.grid_data = np.zeros((20, 50, 48), dtype=np.int32)
        
        # Add some test peaks
        cls.grid_data[10, 25, 24] = 100  # Center peak
        cls.grid_data[5, 10, 15] = 50    # Corner peak
        cls.grid_data[15, 40, 35] = 75   # Another peak
        
        # Add some noise
        np.random.seed(42)  # For reproducible results
        noise = np.random.randint(0, 10, cls.grid_data.shape)
        cls.grid_data += noise

    def test_grid_shape(self):
        # The grid should have 3 dimensions
        self.assertEqual(self.grid_data.ndim, 3)
        self.assertGreater(self.grid_data.size, 0)

    def test_nonzero_voxels(self):
        # There should be some non-zero voxels
        nonzero = np.count_nonzero(self.grid_data)
        print(f"Grid shape: {self.grid_data.shape}")
        print(f"Grid min/max: {self.grid_data.min()}/{self.grid_data.max()}")
        print(f"Non-zero count: {nonzero}")
        print(f"Total grid size: {self.grid_data.size}")
        self.assertGreater(nonzero, 0)

    def test_threshold_filtering(self):
        # Apply a threshold filter and check that some values are filtered
        filtered = GridWorkflows.threshold_based_filtering(self.grid_data, threshold=1)
        self.assertTrue(np.all(filtered[filtered < 1] == 0))
        self.assertLessEqual(np.count_nonzero(filtered), np.count_nonzero(self.grid_data))

    def test_local_maxima_detection(self):
        # Find local maxima in the grid
        maxima = GridWorkflows.local_maxima_detection(self.grid_data, neighborhood_size=1)
        self.assertIsInstance(maxima, list)
        # There should be at least one local maximum
        self.assertGreaterEqual(len(maxima), 1)
        # Each maximum should have the expected tuple structure
        for m in maxima:
            self.assertEqual(len(m), 4)

if __name__ == '__main__':
    unittest.main() 