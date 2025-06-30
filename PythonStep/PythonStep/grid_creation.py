"""
Python implementation of STL_STEPCreateGrid3D::compute() functionality

This module provides a pure Python implementation of the 3D grid creation
and ray tracing algorithm from the C++ STL_STEPCreateGrid3D class.
"""

import numpy as np
from typing import Tuple, Optional, List, Generator
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import time


class RayTracingGrid3D:
    """
    Python implementation of 3D grid creation with ray tracing
    
    This class replicates the functionality of STL_STEPCreateGrid3D::compute()
    using pure Python and numpy for ray tracing and grid filling.
    """
    
    def __init__(self, resolution: float = 0.2):
        self.resolution = resolution
        self.grid_data = None
        self.ray_length_grid = None
        self.bounds = None
        self.dimensions = None
        
    def create_grid_from_point_cloud(self, 
                                   points: np.ndarray, 
                                   normals: np.ndarray,
                                   padding: float = 0.0,
                                   num_threads: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Create and fill 3D grid from point cloud using ray tracing
        
        This is the Python equivalent of STL_STEPCreateGrid3D::compute()
        
        Args:
            points: numpy array of shape (n_points, 3) with X, Y, Z coordinates
            normals: numpy array of shape (n_points, 3) with normal vectors
            padding: Additional padding around the bounding box
            num_threads: Number of threads to use (None for auto-detection)
            
        Returns:
            Tuple of (grid_data, ray_length_grid) where both are 3D numpy arrays
        """
        if num_threads is None:
            num_threads = mp.cpu_count()
            
        print(f"Creating 3D grid with resolution {self.resolution}")
        print(f"Processing {len(points)} points using {num_threads} threads")
        
        # Step 1: Calculate bounding box (like C++ version)
        self.bounds = self._calculate_bounding_box(points, padding)
        min_coords, max_coords = self.bounds
        
        # Step 2: Calculate grid dimensions
        self.dimensions = self._calculate_grid_dimensions(min_coords, max_coords)
        nx, ny, nz = self.dimensions
        
        print(f"Grid dimensions: {nx} x {ny} x {nz}")
        print(f"Grid bounds: {min_coords} to {max_coords}")
        
        # Step 3: Initialize empty grids (like C++ version)
        self.grid_data = np.zeros((nx, ny, nz), dtype=np.int32)
        self.ray_length_grid = np.zeros((nx, ny, nz), dtype=np.float32)
        
        # Step 4: Multi-threaded ray tracing (like C++ version)
        self._process_points_multithreaded(points, normals, num_threads)
        
        # Step 5: Compute statistics (like C++ version)
        self._compute_statistics()
        
        return self.grid_data, self.ray_length_grid
    
    def _calculate_bounding_box(self, points: np.ndarray, padding: float) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate bounding box with padding"""
        min_coords = points.min(axis=0) - padding
        max_coords = points.max(axis=0) + padding
        return min_coords, max_coords
    
    def _calculate_grid_dimensions(self, min_coords: np.ndarray, max_coords: np.ndarray) -> Tuple[int, int, int]:
        """Calculate grid dimensions from bounding box and resolution"""
        dimensions = max_coords - min_coords
        nx = int(np.ceil(dimensions[0] / self.resolution))
        ny = int(np.ceil(dimensions[1] / self.resolution))
        nz = int(np.ceil(dimensions[2] / self.resolution))
        return nx, ny, nz
    
    def _process_points_multithreaded(self, points: np.ndarray, normals: np.ndarray, num_threads: int):
        """Process points using multiple threads (like C++ version)"""
        n_points = len(points)
        points_per_thread = n_points // num_threads
        
        # Create thread pool
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            futures = []
            
            for i in range(num_threads):
                start_idx = i * points_per_thread
                end_idx = start_idx + points_per_thread if i < num_threads - 1 else n_points
                
                future = executor.submit(
                    self._process_point_batch,
                    points[start_idx:end_idx],
                    normals[start_idx:end_idx],
                    i
                )
                futures.append(future)
            
            # Collect results and merge grids
            thread_results = [future.result() for future in futures]
            self._merge_thread_results(thread_results)
    
    def _process_point_batch(self, points: np.ndarray, normals: np.ndarray, thread_id: int) -> Tuple[np.ndarray, np.ndarray]:
        """Process a batch of points (single thread)"""
        nx, ny, nz = self.dimensions
        min_coords, max_coords = self.bounds
        
        # Create local grids for this thread
        local_grid = np.zeros((nx, ny, nz), dtype=np.int32)
        local_ray_length = np.zeros((nx, ny, nz), dtype=np.float32)
        
        print(f"Thread {thread_id}: Processing {len(points)} points")
        
        for i, (point, normal) in enumerate(zip(points, normals)):
            if i % 1000 == 0:
                print(f"Thread {thread_id}: Processed {i}/{len(points)} points")
            
            # Normalize normal vector
            normal_norm = np.linalg.norm(normal)
            if normal_norm == 0:
                continue
                
            normalized_normal = normal / normal_norm
            
            # Create two beams (like C++ version)
            beam_length = 1.5  # Same as C++ version
            end_point_1 = point + normalized_normal * beam_length
            end_point_2 = point - normalized_normal * beam_length
            
            # Trace rays and fill grid
            self._trace_ray_and_fill_grid(point, end_point_1, local_grid, local_ray_length, min_coords)
            self._trace_ray_and_fill_grid(point, end_point_2, local_grid, local_ray_length, min_coords)
        
        return local_grid, local_ray_length
    
    def _trace_ray_and_fill_grid(self, start_point: np.ndarray, end_point: np.ndarray, 
                                grid: np.ndarray, ray_length_grid: np.ndarray, 
                                min_coords: np.ndarray):
        """
        Trace a ray through the grid and fill cells (Python implementation of Woo traversal)
        
        This is a simplified version of STL_Grid3DWooTraversalAlgorithm::compute()
        """
        # Convert to grid coordinates
        start_grid = self._cartesian_to_grid_coords(start_point, min_coords)
        end_grid = self._cartesian_to_grid_coords(end_point, min_coords)
        
        # Use Bresenham-like 3D line algorithm
        grid_cells = self._get_line_cells_3d(start_grid, end_grid)
        
        # Fill each cell along the ray
        for cell in grid_cells:
            x, y, z = cell
            if 0 <= x < grid.shape[0] and 0 <= y < grid.shape[1] and 0 <= z < grid.shape[2]:
                # Increment grid value (like STL_Grid3DBeamVisitor)
                grid[x, y, z] += 1
                
                # Calculate ray length (like C++ version)
                cell_center = self._grid_coords_to_cartesian(cell, min_coords)
                ray_length = np.linalg.norm(cell_center - start_point)
                ray_length_grid[x, y, z] += ray_length
    
    def _cartesian_to_grid_coords(self, point: np.ndarray, min_coords: np.ndarray) -> Tuple[int, int, int]:
        """Convert cartesian coordinates to grid coordinates"""
        grid_coords = np.floor((point - min_coords) / self.resolution).astype(int)
        return tuple(grid_coords)
    
    def _grid_coords_to_cartesian(self, grid_coords: Tuple[int, int, int], min_coords: np.ndarray) -> np.ndarray:
        """Convert grid coordinates to cartesian coordinates (cell center)"""
        return min_coords + (np.array(grid_coords) + 0.5) * self.resolution
    
    def _get_line_cells_3d(self, start: Tuple[int, int, int], end: Tuple[int, int, int]) -> List[Tuple[int, int, int]]:
        """
        Get all grid cells along a 3D line using Bresenham-like algorithm
        
        This is a simplified version of the Woo traversal algorithm
        """
        x0, y0, z0 = start
        x1, y1, z1 = end
        
        dx = abs(x1 - x0)
        dy = abs(y1 - y0)
        dz = abs(z1 - z0)
        
        sx = 1 if x0 < x1 else -1
        sy = 1 if y0 < y1 else -1
        sz = 1 if z0 < z1 else -1
        
        # Find the dominant axis
        if dx >= dy and dx >= dz:
            err_1 = 2 * dy - dx
            err_2 = 2 * dz - dx
            for i in range(dx):
                yield (x0, y0, z0)
                if err_1 > 0:
                    y0 += sy
                    err_1 -= 2 * dx
                if err_2 > 0:
                    z0 += sz
                    err_2 -= 2 * dx
                err_1 += 2 * dy
                err_2 += 2 * dz
                x0 += sx
        elif dy >= dx and dy >= dz:
            err_1 = 2 * dx - dy
            err_2 = 2 * dz - dy
            for i in range(dy):
                yield (x0, y0, z0)
                if err_1 > 0:
                    x0 += sx
                    err_1 -= 2 * dy
                if err_2 > 0:
                    z0 += sz
                    err_2 -= 2 * dy
                err_1 += 2 * dx
                err_2 += 2 * dz
                y0 += sy
        else:
            err_1 = 2 * dy - dz
            err_2 = 2 * dx - dz
            for i in range(dz):
                yield (x0, y0, z0)
                if err_1 > 0:
                    y0 += sy
                    err_1 -= 2 * dz
                if err_2 > 0:
                    x0 += sx
                    err_2 -= 2 * dz
                err_1 += 2 * dy
                err_2 += 2 * dx
                z0 += sz
        
        yield (x1, y1, z1)
    
    def _merge_thread_results(self, thread_results: List[Tuple[np.ndarray, np.ndarray]]):
        """Merge results from all threads (like C++ version)"""
        for local_grid, local_ray_length in thread_results:
            self.grid_data += local_grid
            self.ray_length_grid += local_ray_length
    
    def _compute_statistics(self):
        """Compute and display grid statistics (like C++ version)"""
        if self.grid_data is not None:
            min_val = self.grid_data.min()
            max_val = self.grid_data.max()
            print(f"Grid statistics:")
            print(f"  Min value: {min_val}")
            print(f"  Max value: {max_val}")
            print(f"  Non-zero cells: {np.count_nonzero(self.grid_data)}")
            
        if self.ray_length_grid is not None:
            min_ray = self.ray_length_grid.min()
            max_ray = self.ray_length_grid.max()
            print(f"Ray length statistics:")
            print(f"  Min ray length: {min_ray}")
            print(f"  Max ray length: {max_ray}")


def bresenham_3d(start: np.ndarray, end: np.ndarray) -> Generator[Tuple[int, int, int], None, None]:
    """
    3D Bresenham's line algorithm: yields all voxel indices along a line from start to end (inclusive).
    Args:
        start: (3,) int array, grid coordinates
        end: (3,) int array, grid coordinates
    Yields:
        (x, y, z) tuples
    """
    x0, y0, z0 = start.astype(int)
    x1, y1, z1 = end.astype(int)
    dx = abs(x1 - x0)
    dy = abs(y1 - y0)
    dz = abs(z1 - z0)
    sx = 1 if x1 > x0 else -1
    sy = 1 if y1 > y0 else -1
    sz = 1 if z1 > z0 else -1
    x, y, z = x0, y0, z0
    if dx >= dy and dx >= dz:
        yd = dy - dx // 2
        zd = dz - dx // 2
        for _ in range(dx + 1):
            yield (x, y, z)
            if yd >= 0:
                y += sy
                yd -= dx
            if zd >= 0:
                z += sz
                zd -= dx
            x += sx
            yd += dy
            zd += dz
    elif dy >= dx and dy >= dz:
        xd = dx - dy // 2
        zd = dz - dy // 2
        for _ in range(dy + 1):
            yield (x, y, z)
            if xd >= 0:
                x += sx
                xd -= dy
            if zd >= 0:
                z += sz
                zd -= dy
            y += sy
            xd += dx
            zd += dz
    else:
        xd = dx - dz // 2
        yd = dy - dz // 2
        for _ in range(dz + 1):
            yield (x, y, z)
            if xd >= 0:
                x += sx
                xd -= dz
            if yd >= 0:
                y += sy
                yd -= dz
            z += sz
            xd += dx
            yd += dy


def _trace_ray_and_fill_grid(
    start_point: np.ndarray,
    end_point: np.ndarray,
    grid: np.ndarray,
    ray_length_grid: np.ndarray,
    min_coords: np.ndarray,
    resolution: float
):
    """
    Trace a ray through the grid and fill cells (Python implementation of Woo traversal)
    """
    start_grid = np.floor((start_point - min_coords) / resolution).astype(int)
    end_grid = np.floor((end_point - min_coords) / resolution).astype(int)
    for cell in bresenham_3d(start_grid, end_grid):
        x, y, z = cell
        if 0 <= x < grid.shape[0] and 0 <= y < grid.shape[1] and 0 <= z < grid.shape[2]:
            grid[x, y, z] += 1
            cell_center = min_coords + (np.array([x, y, z]) + 0.5) * resolution
            ray_length = np.linalg.norm(cell_center - start_point)
            ray_length_grid[x, y, z] += ray_length


def _process_point_batch(
    points: np.ndarray,
    normals: np.ndarray,
    grid_shape: Tuple[int, int, int],
    min_coords: np.ndarray,
    resolution: float,
    thread_id: int = 0
) -> Tuple[np.ndarray, np.ndarray]:
    local_grid = np.zeros(grid_shape, dtype=np.int32)
    local_ray_length = np.zeros(grid_shape, dtype=np.float32)
    for i, (point, normal) in enumerate(zip(points, normals)):
        if i % 1000 == 0 and i > 0:
            print(f"Thread {thread_id}: Processed {i}/{len(points)} points")
        normal_norm = np.linalg.norm(normal)
        if normal_norm == 0:
            continue
        normalized_normal = normal / normal_norm
        beam_length = 1.5
        end_point_1 = point + normalized_normal * beam_length
        end_point_2 = point - normalized_normal * beam_length
        _trace_ray_and_fill_grid(point, end_point_1, local_grid, local_ray_length, min_coords, resolution)
        _trace_ray_and_fill_grid(point, end_point_2, local_grid, local_ray_length, min_coords, resolution)
    return local_grid, local_ray_length


def create_grid_from_points_with_ray_tracing(
    points: np.ndarray,
    normals: np.ndarray,
    resolution: float = 0.2,
    padding: float = 0.0,
    num_threads: Optional[int] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    High-level function to create 3D grid from point cloud using ray tracing (multi-threaded).
    Args:
        points: (N, 3) array
        normals: (N, 3) array
        resolution: grid cell size
        padding: extra space around bounding box
        num_threads: number of threads (default: all CPUs)
    Returns:
        (grid_data, ray_length_grid): both are 3D numpy arrays
    """
    if num_threads is None:
        num_threads = mp.cpu_count()
    min_coords = points.min(axis=0) - padding
    max_coords = points.max(axis=0) + padding
    dims = np.ceil((max_coords - min_coords) / resolution).astype(int)
    grid_shape = tuple(dims)
    print(f"Grid shape: {grid_shape}, min: {min_coords}, max: {max_coords}")
    grid_data = np.zeros(grid_shape, dtype=np.int32)
    ray_length_grid = np.zeros(grid_shape, dtype=np.float32)
    n_points = len(points)
    points_per_thread = n_points // num_threads
    batches = []
    for i in range(num_threads):
        start = i * points_per_thread
        end = (i + 1) * points_per_thread if i < num_threads - 1 else n_points
        batches.append((points[start:end], normals[start:end], grid_shape, min_coords, resolution, i))
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(lambda args: _process_point_batch(*args), batches))
    for local_grid, local_ray_length in results:
        grid_data += local_grid
        ray_length_grid += local_ray_length
    print(f"Grid statistics: min={grid_data.min()}, max={grid_data.max()}, nonzero={np.count_nonzero(grid_data)}")
    print(f"Ray length statistics: min={ray_length_grid.min()}, max={ray_length_grid.max()}")
    return grid_data, ray_length_grid 