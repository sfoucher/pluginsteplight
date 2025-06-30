import numpy as np
import matplotlib.pyplot as plt
import os
from PythonStep.grid_creation import create_grid_from_points_with_ray_tracing

# Path to the twoCylinders.txt file (adjust if needed)
DATA_PATH = os.path.join(os.path.dirname(__file__), '../..', 'twoCylinders.txt')

#DATA_PATH = os.path.join(os.path.dirname(__file__), '../..', 'placetteNorm1Arbre.txt')

# Output directory
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), 'output_grid_creation')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 1. Load the point cloud and normals from twoCylinders.txt
def load_two_cylinders_txt(path):
    """
    Assumes twoCylinders.txt has columns: x y z nx ny nz (space or tab separated),
    and may have a header or comment lines starting with // or #.
    Returns:
        points: (N, 3) numpy array
        normals: (N, 3) numpy array
    """
    # Try to skip comment lines and headers
    data = np.loadtxt(path, comments='/', dtype=float)
    points = data[:, :3]
    normals = data[:, 3:6]
    return points, normals

print(f"Loading data from {DATA_PATH}")
points, normals = load_two_cylinders_txt(DATA_PATH)
print(f"Loaded {points.shape[0]} points.")

# 2. Create the 3D grid using ray tracing
grid_resolution = 0.2  # meters (adjust as needed)
padding = 0.0
num_threads = 4  # or None for auto

print("Creating 3D grid from point cloud using ray tracing...")
grid_data, ray_length_grid = create_grid_from_points_with_ray_tracing(
    points, normals, resolution=grid_resolution, padding=padding, num_threads=num_threads)

# 3. Print grid statistics
print(f"Grid shape: {grid_data.shape}")
print(f"Nonzero cells: {np.count_nonzero(grid_data)}")
print(f"Grid min: {grid_data.min()}, max: {grid_data.max()}")

# 4. Visualize a few 2D slices (XY at different Z)
def plot_grid_slices(grid, output_dir, axis=2, num_slices=5, title_prefix="Grid slice"):
    """Save a few 2D slices of the 3D grid as images."""
    nz = grid.shape[axis]
    slice_indices = np.linspace(0, nz-1, num_slices, dtype=int)
    for i, idx in enumerate(slice_indices):
        if axis == 2:
            img = grid[:, :, idx]
        elif axis == 1:
            img = grid[:, idx, :]
        else:
            img = grid[idx, :, :]
        plt.figure(figsize=(6, 5))
        plt.imshow(img.T, origin='lower', cmap='viridis')
        plt.colorbar(label='Votes')
        plt.title(f"{title_prefix} {i+1} (axis={axis}, idx={idx})")
        plt.xlabel('X')
        plt.ylabel('Y')
        out_path = os.path.join(output_dir, f"grid_slice_{axis}_{idx}.png")
        plt.savefig(out_path)
        plt.close()
        print(f"Saved slice image: {out_path}")

plot_grid_slices(grid_data, OUTPUT_DIR, axis=2, num_slices=5, title_prefix="Grid XY slice")

# 5. Save the grid as .npy files
np.save(os.path.join(OUTPUT_DIR, 'grid_data.npy'), grid_data)
np.save(os.path.join(OUTPUT_DIR, 'ray_length_grid.npy'), ray_length_grid)
print(f"Saved grid data and ray length grid as .npy files in {OUTPUT_DIR}") 