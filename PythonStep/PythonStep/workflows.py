"""
Python implementations of STL_Grid3D processing workflows using standard libraries
"""

import numpy as np
from typing import Tuple, List, Optional, Union
from scipy import ndimage
from scipy.spatial import cKDTree
import warnings


class GridWorkflows:
    """
    Python implementations of STL_Grid3D processing workflows
    """
    
    @staticmethod
    def neighbor_based_filtering(grid_data: np.ndarray, neighbors: int = 1, 
                                preserve_zeros: bool = True) -> np.ndarray:
        """
        Filter grid by keeping only local maxima (neighbor-based filtering)
        
        Args:
            grid_data: 3D numpy array representing the grid
            neighbors: Size of neighborhood to check (radius)
            preserve_zeros: Whether to preserve zero values (skip them)
        
        Returns:
            Filtered grid array with same shape as input
        """
        if grid_data.ndim != 3:
            raise ValueError("Grid data must be a 3D array")
        
        filtered_grid = grid_data.copy()
        
        # Create a mask for non-zero values if preserving zeros
        if preserve_zeros:
            non_zero_mask = grid_data > 0
        else:
            non_zero_mask = np.ones_like(grid_data, dtype=bool)
        
        # Get coordinates of non-zero points
        non_zero_coords = np.where(non_zero_mask)
        
        for i, j, k in zip(*non_zero_coords):
            current_value = grid_data[i, j, k]
            
            # Define neighborhood bounds
            i_min, i_max = max(0, i - neighbors), min(grid_data.shape[0], i + neighbors + 1)
            j_min, j_max = max(0, j - neighbors), min(grid_data.shape[1], j + neighbors + 1)
            k_min, k_max = max(0, k - neighbors), min(grid_data.shape[2], k + neighbors + 1)
            
            # Extract neighborhood
            neighborhood = grid_data[i_min:i_max, j_min:j_max, k_min:k_max]
            
            # Check if current value is maximum in neighborhood
            if current_value < neighborhood.max():
                filtered_grid[i, j, k] = 0
        
        return filtered_grid
    
    @staticmethod
    def threshold_based_filtering(grid_data: np.ndarray, threshold: float) -> np.ndarray:
        """
        Filter grid using a fixed threshold value
        
        Args:
            grid_data: 3D numpy array representing the grid
            threshold: Threshold value - values below this are set to 0
        
        Returns:
            Filtered grid array with same shape as input
        """
        filtered_grid = grid_data.copy()
        filtered_grid[filtered_grid < threshold] = 0
        return filtered_grid
    
    @staticmethod
    def percentile_based_filtering(grid_data: np.ndarray, percentile: float) -> np.ndarray:
        """
        Filter grid using percentile-based threshold
        
        Args:
            grid_data: 3D numpy array representing the grid
            percentile: Percentile threshold (0-100)
        
        Returns:
            Filtered grid array with same shape as input
        """
        if not 0 <= percentile <= 100:
            raise ValueError("Percentile must be between 0 and 100")
        
        # Calculate threshold from percentile
        threshold = np.percentile(grid_data[grid_data > 0], percentile)
        return GridWorkflows.threshold_based_filtering(grid_data, threshold)
    
    @staticmethod
    def local_maxima_detection(grid_data: np.ndarray, neighborhood_size: int = 1,
                              min_value: float = 0.0, sort_descending: bool = True) -> List[Tuple[int, int, int, float]]:
        """
        Detect local maxima in the grid
        
        Args:
            grid_data: 3D numpy array representing the grid
            neighborhood_size: Size of neighborhood to check
            min_value: Minimum value to consider as potential maxima
            sort_descending: Whether to sort results by value in descending order
        
        Returns:
            List of tuples (x, y, z, value) representing local maxima positions and values
        """
        if grid_data.ndim != 3:
            raise ValueError("Grid data must be a 3D array")
        
        local_maxima = []
        
        # Create a mask for values above minimum
        value_mask = grid_data > min_value
        
        # Get coordinates of points above minimum
        coords = np.where(value_mask)
        
        for i, j, k in zip(*coords):
            current_value = grid_data[i, j, k]
            
            # Define neighborhood bounds
            i_min, i_max = max(0, i - neighborhood_size), min(grid_data.shape[0], i + neighborhood_size + 1)
            j_min, j_max = max(0, j - neighborhood_size), min(grid_data.shape[1], j + neighborhood_size + 1)
            k_min, k_max = max(0, k - neighborhood_size), min(grid_data.shape[2], k + neighborhood_size + 1)
            
            # Extract neighborhood
            neighborhood = grid_data[i_min:i_max, j_min:j_max, k_min:k_max]
            
            # Check if current value is maximum in neighborhood
            if current_value >= neighborhood.max():
                local_maxima.append((i, j, k, current_value))
        
        # Sort by value if requested
        if sort_descending:
            local_maxima.sort(key=lambda x: x[3], reverse=True)
        
        return local_maxima
    
    @staticmethod
    def local_maxima_in_bbox(grid_data: np.ndarray, bbox_min: Tuple[int, int, int], 
                            bbox_max: Tuple[int, int, int], neighborhood_size: int = 1,
                            min_value: float = 0.0, sort_descending: bool = True) -> List[Tuple[int, int, int, float]]:
        """
        Detect local maxima within a bounding box
        
        Args:
            grid_data: 3D numpy array representing the grid
            bbox_min: Minimum coordinates (x, y, z) of bounding box
            bbox_max: Maximum coordinates (x, y, z) of bounding box
            neighborhood_size: Size of neighborhood to check
            min_value: Minimum value to consider as potential maxima
            sort_descending: Whether to sort results by value in descending order
        
        Returns:
            List of tuples (x, y, z, value) representing local maxima positions and values
        """
        # Validate bounding box
        for i in range(3):
            if bbox_min[i] < 0 or bbox_max[i] > grid_data.shape[i]:
                raise ValueError(f"Bounding box {bbox_min} to {bbox_max} is outside grid bounds {grid_data.shape}")
            if bbox_min[i] >= bbox_max[i]:
                raise ValueError(f"Invalid bounding box: min {bbox_min} >= max {bbox_max}")
        
        # Extract subgrid
        subgrid = grid_data[bbox_min[0]:bbox_max[0], 
                           bbox_min[1]:bbox_max[1], 
                           bbox_min[2]:bbox_max[2]]
        
        # Find local maxima in subgrid
        local_maxima = GridWorkflows.local_maxima_detection(
            subgrid, neighborhood_size, min_value, sort_descending
        )
        
        # Adjust coordinates back to original grid
        adjusted_maxima = []
        for x, y, z, value in local_maxima:
            adjusted_maxima.append((
                x + bbox_min[0],
                y + bbox_min[1], 
                z + bbox_min[2],
                value
            ))
        
        return adjusted_maxima
    
    @staticmethod
    def gaussian_filtering(grid_data: np.ndarray, sigma: float = 1.0) -> np.ndarray:
        """
        Apply Gaussian smoothing to the grid
        
        Args:
            grid_data: 3D numpy array representing the grid
            sigma: Standard deviation for Gaussian kernel
        
        Returns:
            Smoothed grid array with same shape as input
        """
        return ndimage.gaussian_filter(grid_data, sigma=sigma)
    
    @staticmethod
    def morphological_filtering(grid_data: np.ndarray, operation: str = 'opening', 
                              size: int = 1) -> np.ndarray:
        """
        Apply morphological operations to the grid
        
        Args:
            grid_data: 3D numpy array representing the grid
            operation: Morphological operation ('opening', 'closing', 'erosion', 'dilation')
            size: Size of the structuring element
        
        Returns:
            Filtered grid array with same shape as input
        """
        # Create binary mask for non-zero values
        binary_mask = grid_data > 0
        
        if operation == 'opening':
            filtered_mask = ndimage.binary_opening(binary_mask, 
                                                 structure=np.ones((size, size, size)))
        elif operation == 'closing':
            filtered_mask = ndimage.binary_closing(binary_mask, 
                                                 structure=np.ones((size, size, size)))
        elif operation == 'erosion':
            filtered_mask = ndimage.binary_erosion(binary_mask, 
                                                 structure=np.ones((size, size, size)))
        elif operation == 'dilation':
            filtered_mask = ndimage.binary_dilation(binary_mask, 
                                                  structure=np.ones((size, size, size)))
        else:
            raise ValueError(f"Unknown operation: {operation}")
        
        # Apply mask to original data
        filtered_grid = grid_data.copy()
        filtered_grid[~filtered_mask] = 0
        
        return filtered_grid
    
    @staticmethod
    def connected_components_analysis(grid_data: np.ndarray, min_size: int = 10) -> Tuple[np.ndarray, int]:
        """
        Perform connected components analysis on the grid
        
        Args:
            grid_data: 3D numpy array representing the grid
            min_size: Minimum size of connected components to keep
        
        Returns:
            Tuple of (labeled_grid, num_components) where:
            - labeled_grid: Grid with each component labeled with a unique integer
            - num_components: Number of connected components found
        """
        # Create binary mask
        binary_mask = grid_data > 0
        
        # Label connected components
        labeled_grid, num_components = ndimage.label(binary_mask)
        
        # Remove small components
        if min_size > 1:
            component_sizes = ndimage.sum(binary_mask, labeled_grid, range(1, num_components + 1))
            small_components = np.where(component_sizes < min_size)[0] + 1
            
            for component_id in small_components:
                labeled_grid[labeled_grid == component_id] = 0
            
            # Relabel remaining components
            labeled_grid, num_components = ndimage.label(labeled_grid > 0)
        
        return labeled_grid, num_components
    
    @staticmethod
    def distance_transform(grid_data: np.ndarray, metric: str = 'euclidean') -> np.ndarray:
        """
        Compute distance transform of the grid
        
        Args:
            grid_data: 3D numpy array representing the grid
            metric: Distance metric ('euclidean', 'manhattan', 'chessboard')
        
        Returns:
            Distance transform array with same shape as input
        """
        # Create binary mask
        binary_mask = grid_data > 0
        
        if metric == 'euclidean':
            return ndimage.distance_transform_edt(~binary_mask)
        elif metric == 'manhattan':
            return ndimage.distance_transform_cdt(~binary_mask, metric='taxicab')
        elif metric == 'chessboard':
            return ndimage.distance_transform_cdt(~binary_mask, metric='chessboard')
        else:
            raise ValueError(f"Unknown metric: {metric}")
    
    @staticmethod
    def watershed_segmentation(grid_data: np.ndarray, markers: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Perform watershed segmentation on the grid
        
        Args:
            grid_data: 3D numpy array representing the grid
            markers: Optional marker array for watershed seeds
        
        Returns:
            Segmented grid array with same shape as input
        """
        from scipy import ndimage
        
        # Invert the data (watershed works on minima)
        inverted_data = -grid_data
        
        if markers is None:
            # Use local minima as markers
            markers = ndimage.label(ndimage.minimum_filter(inverted_data, size=3) == inverted_data)[0]
        
        # Apply watershed
        segmented = ndimage.watershed_ift(inverted_data.astype(np.uint8), markers)
        
        return segmented
    
    @staticmethod
    def grid_statistics(grid_data: np.ndarray) -> dict:
        """
        Compute comprehensive statistics for the grid
        
        Args:
            grid_data: 3D numpy array representing the grid
        
        Returns:
            Dictionary containing various statistics
        """
        non_zero_data = grid_data[grid_data > 0]
        
        stats = {
            'shape': grid_data.shape,
            'total_cells': grid_data.size,
            'non_zero_cells': len(non_zero_data),
            'sparsity': 1.0 - len(non_zero_data) / grid_data.size,
            'min_value': float(grid_data.min()),
            'max_value': float(grid_data.max()),
            'mean_value': float(grid_data.mean()),
            'std_value': float(grid_data.std()),
        }
        
        if len(non_zero_data) > 0:
            stats.update({
                'non_zero_min': float(non_zero_data.min()),
                'non_zero_max': float(non_zero_data.max()),
                'non_zero_mean': float(non_zero_data.mean()),
                'non_zero_std': float(non_zero_data.std()),
                'non_zero_median': float(np.median(non_zero_data)),
                'non_zero_percentiles': {
                    '25': float(np.percentile(non_zero_data, 25)),
                    '50': float(np.percentile(non_zero_data, 50)),
                    '75': float(np.percentile(non_zero_data, 75)),
                    '90': float(np.percentile(non_zero_data, 90)),
                    '95': float(np.percentile(non_zero_data, 95)),
                    '99': float(np.percentile(non_zero_data, 99)),
                }
            })
        
        return stats


class GridVisualization:
    """
    Visualization utilities for grid data
    """
    
    @staticmethod
    def create_2d_slice(grid_data: np.ndarray, axis: int = 2, slice_index: Optional[int] = None) -> np.ndarray:
        """
        Create a 2D slice from 3D grid data
        
        Args:
            grid_data: 3D numpy array
            axis: Axis along which to slice (0, 1, or 2)
            slice_index: Index for slicing, if None uses middle
        
        Returns:
            2D numpy array representing the slice
        """
        if slice_index is None:
            slice_index = grid_data.shape[axis] // 2
        
        if axis == 0:
            return grid_data[slice_index, :, :]
        elif axis == 1:
            return grid_data[:, slice_index, :]
        elif axis == 2:
            return grid_data[:, :, slice_index]
        else:
            raise ValueError("Axis must be 0, 1, or 2")
    
    @staticmethod
    def create_2d_projection(grid_data: np.ndarray, axis: int = 2, 
                           projection_type: str = 'max') -> np.ndarray:
        """
        Create a 2D projection from 3D grid data
        
        Args:
            grid_data: 3D numpy array
            axis: Axis along which to project (0, 1, or 2)
            projection_type: Type of projection ('max', 'min', 'mean', 'sum')
        
        Returns:
            2D numpy array representing the projection
        """
        if projection_type == 'max':
            return np.max(grid_data, axis=axis)
        elif projection_type == 'min':
            return np.min(grid_data, axis=axis)
        elif projection_type == 'mean':
            return np.mean(grid_data, axis=axis)
        elif projection_type == 'sum':
            return np.sum(grid_data, axis=axis)
        else:
            raise ValueError(f"Unknown projection type: {projection_type}") 