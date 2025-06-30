# PythonStep Examples

This directory contains example scripts demonstrating how to use the PythonStep package for 3D point cloud and grid processing workflows. All examples use the sample file `twoCylinders.txt`.

## Example Scripts

### 1. `basic_workflow_example.py`
- **Purpose:** Demonstrates a complete workflow: loading a point cloud, creating a grid, applying filtering and analysis, and visualizing results.
- **Features:**
  - Loads `twoCylinders.txt`
  - Creates a synthetic 3D grid
  - Applies threshold filtering, neighbor-based filtering, local maxima detection, smoothing, and connected components analysis
  - Visualizes results with matplotlib (6-panel figure)
  - **Output:** `workflow_results.png`

![Workflow Results](workflow_results.png)
*Complete workflow visualization showing original grid, filtered results, and analysis outputs*

### 2. `simple_workflow_example.py`
- **Purpose:** Lightweight version of the workflow for environments without matplotlib.
- **Features:**
  - Same workflow as above, but prints results and saves to text files
  - **Output:**
    - `local_maxima_results.txt`: Coordinates and values of detected peaks
    - `grid_statistics.txt`: Comprehensive grid statistics

### 3. `visualize_3d_dataset.py`
- **Purpose:** Visualizes the 3D dataset in `twoCylinders.txt` using various techniques.
- **Features:**
  - 3D scatter plot of points (with color by Z)
  - 3D scatter plot with normal vectors
  - 2D density projections (XY, XZ, YZ)
  - Statistical visualizations (histograms, scatter plots)
  - K-means clustering visualization
  - Interactive 3D plot (if supported)
  - **Output:**
    - `3d_scatter_plot.png`
    - `density_projections.png`
    - `statistical_analysis.png`
    - `cluster_visualization.png`
    - `interactive_3d_plot.png` (if interactive backend is available)
    - `dataset_summary.txt`

### 4. `active_contours_example.py`
- **Purpose:** Demonstrates the Active Contours workflow for extracting characteristic curves from 3D grid data.
- **Features:**
  - Creates synthetic 3D grid data with curve-like structures (straight line, helix, S-curve)
  - Applies the active contours algorithm to extract curves
  - Visualizes the grid data and extracted curves in 3D and 2D projections
  - Analyzes curve properties (length, number of points, bounding box)
  - **Output:**
    - `active_contours_visualization.png` - 3D and 2D visualizations of grid and curves
    - `extracted_curves.txt` - Coordinates of extracted curve points
    - `grid_statistics.txt` - Grid analysis statistics

![Active Contours Visualization](active_contours_visualization.png)
*Active contours workflow showing synthetic grid data and extracted characteristic curves*

### 5. `active_contours_two_cylinders_example.py`
**Purpose:** Demonstrates the Active Contours workflow using real point cloud data from `twoCylinders.txt`.
**Features:**
  - Loads point cloud data from `twoCylinders.txt`
  - Creates 3D grid from real point cloud data
  - Applies active contours algorithm to extract curves from real data
  - Enhanced visualization with multiple 2D projections (XY, XZ, YZ)
  - Comprehensive analysis report with point cloud statistics
  - **Output:**
    - `two_cylinders_active_contours_visualization.png` - 3D and 2D visualizations of grid and curves
    - `two_cylinders_extracted_curves.txt` - Coordinates of extracted curve points
    - `two_cylinders_analysis_report.txt` - Comprehensive analysis report with point cloud and grid statistics

### 6. `grid_creation_two_cylinders_example.py`
- **Purpose:** Demonstrates the Python implementation of STL_STEPCreateGrid3D::compute() using ray tracing to fill a 3D grid from point cloud data.
- **Features:**
  - Loads point cloud and normals from `twoCylinders.txt`
  - Uses multi-threaded ray tracing algorithm (3D Bresenham/Woo traversal)
  - Creates filled 3D grid by tracing rays along normal vectors
  - Visualizes 2D slices of the resulting grid at different Z levels
  - Saves grid data and ray length information
  - **Output:**
    - `output_grid_creation/` directory containing:
      - `grid_slice_2_*.png` - 2D slice visualizations at different Z levels
      - `grid_data.npy` - 3D grid data as numpy array
      - `ray_length_grid.npy` - Ray length information for each grid cell

![Grid Creation Visualization](output_grid_creation/grid_slice_2_0.png)
*2D slice visualization showing the filled grid from ray tracing*

#### Visualization Examples

**3D Scatter Plot with Normal Vectors:**
![3D Scatter Plot](3d_scatter_plot.png)
*3D visualization of the two cylinders dataset with normal vectors*

**Density Projections:**
![Density Projections](density_projections.png)
*2D density projections showing the spatial distribution of points*

**Statistical Analysis:**
![Statistical Analysis](statistical_analysis.png)
*Statistical distributions and coordinate relationships*

**Cluster Visualization:**
![Cluster Visualization](cluster_visualization.png)
*K-means clustering showing the two distinct cylinder clusters*

## Requirements
- Python 3.7+
- numpy
- scipy
- matplotlib

Install requirements with:
```bash
pip install numpy scipy matplotlib
```

## Usage
Run any example from this directory:

```bash
python3 basic_workflow_example.py
python3 simple_workflow_example.py
python3 visualize_3d_dataset.py
python3 active_contours_example.py
python3 active_contours_two_cylinders_example.py
python3 grid_creation_two_cylinders_example.py
```

## Data File
All examples require the file `twoCylinders.txt` to be present at the path `/home/sfoucher/DEV/pluginsteplight/twoCylinders.txt`.

## Output
Each script will generate output files (images, text) in this directory. See the script comments and the above descriptions for details.

### Sample Output Files
After running the examples, you should see files like:
- `workflow_results.png` - Complete workflow visualization
- `3d_scatter_plot.png` - 3D point cloud visualization
- `density_projections.png` - 2D density projections
- `statistical_analysis.png` - Statistical distributions
- `cluster_visualization.png` - Clustering results
- `active_contours_visualization.png` - Active contours curve extraction
- `two_cylinders_active_contours_visualization.png` - Active contours on real data
- `local_maxima_results.txt` - Detected peak coordinates
- `extracted_curves.txt` - Extracted curve coordinates
- `grid_statistics.txt` - Grid analysis statistics
- `dataset_summary.txt` - Dataset overview
- `output_grid_creation/` - Grid creation results with slice images and .npy files

---
For more information, see the main project README or the script docstrings. 