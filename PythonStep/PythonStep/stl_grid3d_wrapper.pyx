# cython: language_level=3
# distutils: language=c++

"""
Cython wrapper for STL_Grid3D class - Simplified version for testing
"""

import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from cpython.ref cimport PyObject

# C++ declarations are in the .pxd file

# Python wrapper class
cdef class PySTLGrid3D:
    """
    Python wrapper for STL_Grid3D class
    """
    # Attributes _grid and _grid_allocated are declared in the .pxd file
    
    def __cinit__(self, double resolution=0.1):
        self._grid = NULL
        self._grid_allocated = False
    
    def __dealloc__(self):
        if self._grid_allocated and self._grid != NULL:
            del self._grid
    
    def create_from_bounds(self, double xmin, double ymin, double zmin,
                          double xmax, double ymax, double zmax,
                          double resolution, int na_value=-9999, int init_value=0):
        """
        Create grid from bounding box coordinates
        """
        if self._grid_allocated and self._grid != NULL:
            del self._grid
            
        # Use a temporary instance to call createGrid3DFromXYZCoords
        tmp = STL_Grid3D_int()
        self._grid = tmp.createGrid3DFromXYZCoords(
            xmin, ymin, zmin, xmax, ymax, zmax, 
            resolution, na_value, init_value, True)
        self._grid_allocated = True
        return self
    
    def create_from_dimensions(self, double xmin, double ymin, double zmin,
                              size_t dimx, size_t dimy, size_t dimz,
                              double resolution, int na_value=-9999, int init_value=0):
        """
        Create grid from dimensions
        """
        if self._grid_allocated and self._grid != NULL:
            del self._grid
            
        self._grid = new STL_Grid3D_int(xmin, ymin, zmin, dimx, dimy, dimz,
                                       resolution, na_value, init_value)
        self._grid_allocated = True
        return self
    
    def get_data(self):
        """
        Get grid data as numpy array
        """
        if self._grid == NULL:
            raise RuntimeError("Grid not initialized")
        
        cdef int* data_ptr = self._grid.get_data()
        cdef int nx = self._grid.xdim()
        cdef int ny = self._grid.ydim()
        cdef int nz = self._grid.zdim()
        
        # Create numpy array from C++ data
        cdef np.ndarray[np.int32_t, ndim=3] data = np.zeros((nx, ny, nz), dtype=np.int32)
        
        cdef int i, j, k
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    data[i, j, k] = data_ptr[i * ny * nz + j * nz + k]
        
        return data
    
    def get_value(self, int x, int y, int z):
        """
        Get value at specific grid coordinates
        """
        if self._grid == NULL:
            raise RuntimeError("Grid not initialized")
        
        if x < 0 or x >= self._grid.xdim() or \
           y < 0 or y >= self._grid.ydim() or \
           z < 0 or z >= self._grid.zdim():
            raise IndexError("Coordinates out of bounds")
        
        return self._grid.value(x, y, z)
    
    def set_value(self, int x, int y, int z, int value):
        """
        Set value at specific grid coordinates
        """
        if self._grid == NULL:
            raise RuntimeError("Grid not initialized")
        
        if x < 0 or x >= self._grid.xdim() or \
           y < 0 or y >= self._grid.ydim() or \
           z < 0 or z >= self._grid.zdim():
            raise IndexError("Coordinates out of bounds")
        
        # Access the data directly
        cdef int* data_ptr = self._grid.get_data()
        cdef int ny = self._grid.ydim()
        cdef int nz = self._grid.zdim()
        data_ptr[x * ny * nz + y * nz + z] = value
    
    def pixel_to_cartesian(self, int x, int y, int z):
        """
        Convert pixel coordinates to cartesian coordinates
        """
        if self._grid == NULL:
            raise RuntimeError("Grid not initialized")
        
        return (self._grid.pixelToCartesianX(x),
                self._grid.pixelToCartesianY(y),
                self._grid.pixelToCartesianZ(z))
    
    def cartesian_to_pixel(self, double x, double y, double z):
        """
        Convert cartesian coordinates to pixel coordinates
        """
        if self._grid == NULL:
            raise RuntimeError("Grid not initialized")
        
        return (self._grid.cartesianToPixelX(x),
                self._grid.cartesianToPixelY(y),
                self._grid.cartesianToPixelZ(z))
    
    @property
    def dimensions(self):
        """Get grid dimensions"""
        if self._grid == NULL:
            raise RuntimeError("Grid not initialized")
        return (self._grid.xdim(), self._grid.ydim(), self._grid.zdim())
    
    @property
    def resolution(self):
        """Get grid resolution"""
        if self._grid == NULL:
            raise RuntimeError("Grid not initialized")
        return self._grid.resolution()
    
    @property
    def bounds(self):
        """Get grid bounds"""
        if self._grid == NULL:
            raise RuntimeError("Grid not initialized")
        return (self._grid.minX(), self._grid.minY(), self._grid.minZ(),
                self._grid.maxX(), self._grid.maxY(), self._grid.maxZ())

# Utility function to create grid from text file
def create_grid_from_txt(str filename, double resolution=0.1, int na_value=-9999):
    """
    Create STL_Grid3D from twoCylinders.txt format file
    
    Args:
        filename: Path to the text file
        resolution: Grid resolution
        na_value: Value to use for NA/missing data
    
    Returns:
        PySTLGrid3D: Initialized grid object
    """
    # Read the file
    data = np.loadtxt(filename, skiprows=1)  # Skip comment line
    
    if data.shape[1] != 6:
        raise ValueError("File must have 6 columns: X Y Z Nx Ny Nz")
    
    # Extract points and normals
    points = data[:, :3]  # X Y Z
    normals = data[:, 3:]  # Nx Ny Nz
    
    # Calculate bounding box
    xmin, ymin, zmin = points.min(axis=0)
    xmax, ymax, zmax = points.max(axis=0)
    
    # Create grid
    grid = PySTLGrid3D(resolution)
    grid.create_from_bounds(xmin, ymin, zmin, xmax, ymax, zmax, 
                           resolution, na_value, 0)
    
    # TODO: Set point cloud and normals
    # This would require implementing the point cloud wrapper
    
    return grid 