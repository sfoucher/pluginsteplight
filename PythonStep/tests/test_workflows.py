import unittest
import numpy as np
from PythonStep.workflows import GridWorkflows

class TestGridWorkflows(unittest.TestCase):
    def setUp(self):
        # Create a simple 3D grid for testing
        self.grid = np.zeros((5, 5, 5), dtype=float)
        self.grid[2, 2, 2] = 10
        self.grid[1, 1, 1] = 5
        self.grid[3, 3, 3] = 7
        self.grid[4, 4, 4] = 2
        self.grid[0, 0, 0] = 1

    def test_neighbor_based_filtering(self):
        filtered = GridWorkflows.neighbor_based_filtering(self.grid, neighbors=1)
        # Only the global maximum should remain
        self.assertEqual(filtered[2, 2, 2], 10)
        self.assertTrue(np.all(filtered[(self.grid < 10)] == 0))

    def test_threshold_based_filtering(self):
        filtered = GridWorkflows.threshold_based_filtering(self.grid, threshold=6)
        self.assertEqual(filtered[2, 2, 2], 10)
        self.assertEqual(filtered[3, 3, 3], 7)
        self.assertEqual(filtered[1, 1, 1], 0)
        self.assertEqual(filtered[0, 0, 0], 0)

    def test_local_maxima_detection(self):
        maxima = GridWorkflows.local_maxima_detection(self.grid, neighborhood_size=1)
        # Should find the three local maxima (10, 7, 5)
        values = [v for _, _, _, v in maxima]
        self.assertIn(10, values)
        self.assertIn(7, values)
        self.assertIn(5, values)
        self.assertNotIn(2, values)
        self.assertNotIn(1, values)

    def test_local_maxima_in_bbox(self):
        maxima = GridWorkflows.local_maxima_in_bbox(self.grid, (1, 1, 1), (4, 4, 4), neighborhood_size=1)
        # Only the maxima within the bounding box (should include 10 and 7)
        values = [v for _, _, _, v in maxima]
        self.assertIn(10, values)
        self.assertIn(7, values)
        self.assertNotIn(5, values)

    def test_percentile_based_filtering(self):
        filtered = GridWorkflows.percentile_based_filtering(self.grid, percentile=80)
        # Only the top 20% values should remain (10 and 7)
        self.assertEqual(filtered[2, 2, 2], 10)
        self.assertEqual(filtered[3, 3, 3], 7)
        self.assertEqual(filtered[1, 1, 1], 0)
        self.assertEqual(filtered[0, 0, 0], 0)

if __name__ == '__main__':
    unittest.main() 