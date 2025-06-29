# cython: language_level=3

"""
Cython declarations for STL_Grid3D wrapper
"""

import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

# Simplified C++ declarations for testing without COMPUTREE dependencies
cdef extern from *:
    """
    // Simplified STL_Grid3D class for testing
    #include <vector>
    #include <cmath>
    
    template<typename T>
    class STL_Grid3D {
    private:
        std::vector<T> _data;
        int _dimx, _dimy, _dimz;
        double _res, _minX, _minY, _minZ, _maxX, _maxY, _maxZ;
        T _na_value;
    
    public:
        STL_Grid3D() : _dimx(0), _dimy(0), _dimz(0), _res(0.1), 
                      _minX(0), _minY(0), _minZ(0), _maxX(1), _maxY(1), _maxZ(1), _na_value(-9999) {}
        
        STL_Grid3D(double xmin, double ymin, double zmin, 
                   size_t dimx, size_t dimy, size_t dimz,
                   double res, T na, T initValue) 
            : _dimx(dimx), _dimy(dimy), _dimz(dimz), _res(res),
              _minX(xmin), _minY(ymin), _minZ(zmin), 
              _maxX(xmin + dimx * res), _maxY(ymin + dimy * res), _maxZ(zmin + dimz * res),
              _na_value(na) {
            _data.resize(dimx * dimy * dimz, initValue);
        }
        
        STL_Grid3D(double xmin, double ymin, double zmin,
                   double xmax, double ymax, double zmax,
                   double res, T na, T initValue)
            : _res(res), _minX(xmin), _minY(ymin), _minZ(zmin), 
              _maxX(xmax), _maxY(ymax), _maxZ(zmax), _na_value(na) {
            _dimx = (int)std::ceil((xmax - xmin) / res);
            _dimy = (int)std::ceil((ymax - ymin) / res);
            _dimz = (int)std::ceil((zmax - zmin) / res);
            _data.resize(_dimx * _dimy * _dimz, initValue);
        }
        
        static STL_Grid3D<T>* createGrid3DFromXYZCoords(
            double xmin, double ymin, double zmin,
            double xmax, double ymax, double zmax,
            double resolution, T na, T initValue, bool extends) {
            return new STL_Grid3D<T>(xmin, ymin, zmin, xmax, ymax, zmax, resolution, na, initValue);
        }
        
        void setPointCloudPtr(void* point_cloud_ptr, void* normal_cloud_ptr) {
            // Placeholder for point cloud setting
        }
        
        T* get_data() { return _data.data(); }
        
        T value(int x, int y, int z) const {
            if (x < 0 || x >= _dimx || y < 0 || y >= _dimy || z < 0 || z >= _dimz)
                return _na_value;
            return _data[x * _dimy * _dimz + y * _dimz + z];
        }
        
        double pixelToCartesianX(int x) const { return _minX + (x + 0.5) * _res; }
        double pixelToCartesianY(int y) const { return _minY + (y + 0.5) * _res; }
        double pixelToCartesianZ(int z) const { return _minZ + (z + 0.5) * _res; }
        
        int cartesianToPixelX(double x) const {
            if (x == _maxX) return _dimx - 1;
            return (int)std::floor((x - _minX) / _res);
        }
        
        int cartesianToPixelY(double y) const {
            if (y == _maxY) return _dimy - 1;
            return (int)std::floor((y - _minY) / _res);
        }
        
        int cartesianToPixelZ(double z) const {
            if (z == _maxZ) return _dimz - 1;
            return (int)std::floor((z - _minZ) / _res);
        }
        
        int xdim() const { return _dimx; }
        int ydim() const { return _dimy; }
        int zdim() const { return _dimz; }
        double resolution() const { return _res; }
        double minX() const { return _minX; }
        double minY() const { return _minY; }
        double minZ() const { return _minZ; }
        double maxX() const { return _maxX; }
        double maxY() const { return _maxY; }
        double maxZ() const { return _maxZ; }
    };
    """
    
    cdef cppclass STL_Grid3D_int "STL_Grid3D<int>":
        STL_Grid3D_int()
        STL_Grid3D_int(double xmin, double ymin, double zmin, 
                      size_t dimx, size_t dimy, size_t dimz,
                      double res, int na, int initValue)
        STL_Grid3D_int(double xmin, double ymin, double zmin,
                      double xmax, double ymax, double zmax,
                      double res, int na, int initValue)
        
        # No 'static' keyword here for Cython compatibility
        STL_Grid3D_int* createGrid3DFromXYZCoords(
            double xmin, double ymin, double zmin,
            double xmax, double ymax, double zmax,
            double resolution, int na, int initValue, bool extends)
        
        # Methods
        void setPointCloudPtr(void* point_cloud_ptr, void* normal_cloud_ptr)
        int* get_data()
        int value(int x, int y, int z) const
        double pixelToCartesianX(int x) const
        double pixelToCartesianY(int y) const  
        double pixelToCartesianZ(int z) const
        int cartesianToPixelX(double x) const
        int cartesianToPixelY(double y) const
        int cartesianToPixelZ(double z) const
        int xdim() const
        int ydim() const
        int zdim() const
        double resolution() const
        double minX() const
        double minY() const
        double minZ() const
        double maxX() const
        double maxY() const
        double maxZ() const

# Note: Eigen and OpenCV declarations removed for simplified testing version

# Python wrapper class declaration
cdef class PySTLGrid3D:
    cdef STL_Grid3D_int* _grid
    cdef bool _grid_allocated 