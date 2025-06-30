"""
Python implementation of Open Active Contours workflow for STL_Grid3D
Based on the C++ STL_OpenActiveContours implementation
"""

import numpy as np
from typing import List, Tuple, Optional, Union
from scipy.spatial import cKDTree
from scipy.linalg import inv
from scipy.ndimage import gaussian_filter
import warnings
from dataclasses import dataclass
from enum import Enum


class GrowthDirection(Enum):
    """Enumeration for growth directions"""
    HEAD = "head"
    TAIL = "tail"
    BOTH = "both"
    NONE = "none"


@dataclass
class ContourPoint:
    """Represents a point in the active contour"""
    x: float
    y: float
    z: float
    value: float = 0.0
    
    def to_array(self) -> np.ndarray:
        """Convert to numpy array"""
        return np.array([self.x, self.y, self.z])
    
    @classmethod
    def from_array(cls, arr: np.ndarray, value: float = 0.0) -> 'ContourPoint':
        """Create from numpy array"""
        return cls(arr[0], arr[1], arr[2], value)


class OpenActiveContours:
    """
    Python implementation of Open Active Contours for 3D grid processing
    
    This class implements the core functionality of the C++ STL_OpenActiveContours
    for extracting characteristic curves from 3D grid data.
    """
    
    def __init__(self, grid_data: np.ndarray, resolution: float = 1.0,
                 min_value: float = 10.0, search_cone_size: int = 20):
        """
        Initialize the open active contours
        
        Args:
            grid_data: 3D numpy array representing the grid
            resolution: Grid resolution in units per pixel
            min_value: Minimum value to consider for growth
            search_cone_size: Size of search cone for direction finding
        """
        self.grid_data = grid_data
        self.resolution = resolution
        self.min_value = min_value
        self.search_cone_size = search_cone_size
        
        # Initialize contour points
        self.points: List[ContourPoint] = []
        
        # Repulsion grid to prevent overlap
        self.repulsion_grid = np.zeros_like(grid_data, dtype=bool)
        
        # Growth parameters
        self.grow_coeff = 1.0
        self.cone_angle_max_degrees = 20.0
        self.cone_size_max_pixels = 3
        self.seuil_sigma_l1 = 0.75
        
        # Relaxation parameters
        self.alpha = 1.5  # Elasticity
        self.beta = 1.5   # Bending stiffness
        self.gamma = 1.3  # Image attraction
        self.global_weight = 1.0
        self.time_step = 0.001
        self.thresh_average_movement = 0.001
        
    def initialize_from_point(self, start_pixel: Tuple[int, int, int]) -> bool:
        """
        Initialize the contour from a starting pixel
        
        Args:
            start_pixel: Starting pixel coordinates (x, y, z)
            
        Returns:
            True if initialization successful, False otherwise
        """
        x, y, z = start_pixel
        
        # Check if starting point is valid
        if not self._is_valid_pixel(x, y, z):
            return False
            
        if self.grid_data[x, y, z] < self.min_value:
            return False
            
        # Convert pixel to cartesian coordinates
        start_point = self._pixel_to_cartesian(x, y, z)
        
        # Find initial growth direction using PCA
        growth_direction = self._find_initial_direction(start_pixel, start_point)
        
        if growth_direction is None:
            return False
            
        # Initialize contour with two points
        head_point = start_point - growth_direction * self.resolution
        tail_point = start_point + growth_direction * self.resolution
        
        self.points = [
            ContourPoint(head_point[0], head_point[1], head_point[2]),
            ContourPoint(tail_point[0], tail_point[1], tail_point[2])
        ]
        
        # Resample to get more points
        self.resample(self.length_3d() / 3.0)
        
        return True
    
    def _find_initial_direction(self, center_pixel: Tuple[int, int, int], 
                               reference_point: np.ndarray) -> Optional[np.ndarray]:
        """
        Find initial growth direction using PCA on potential directions
        
        Args:
            center_pixel: Center pixel coordinates
            reference_point: Reference point in cartesian coordinates
            
        Returns:
            Normalized growth direction vector or None if failed
        """
        # Get potential directions and scores
        directions_and_scores = self._get_directions_and_scores_in_bbox(
            center_pixel, reference_point, self.search_cone_size, 0, 40,
            self.min_value, True
        )
        
        if len(directions_and_scores) < 3:
            return None
            
        # Extract directions and weights
        directions = np.array([d for d, _ in directions_and_scores])
        weights = np.array([s for _, s in directions_and_scores])
        
        # Perform weighted PCA
        weighted_mean = np.average(directions, weights=weights, axis=0)
        centered_directions = directions - weighted_mean
        
        # Compute covariance matrix
        cov_matrix = np.cov(centered_directions.T, aweight=weights)
        
        # Get principal components
        eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
        
        # Return the direction of maximum variance
        principal_direction = eigenvectors[:, -1]
        
        # Normalize and scale by resolution
        norm = np.linalg.norm(principal_direction)
        if norm == 0:
            return None
            
        return principal_direction * self.resolution / norm
    
    def _get_directions_and_scores_in_bbox(self, center_pixel: Tuple[int, int, int],
                                         reference_point: np.ndarray,
                                         search_bbox_size: int, search_bbox_radius: int,
                                         n_max_dir_to_keep: int, min_value: float,
                                         sort_crescent_order: bool) -> List[Tuple[np.ndarray, float]]:
        """
        Get potential growth directions and scores within a bounding box
        
        Args:
            center_pixel: Center pixel coordinates
            reference_point: Reference point in cartesian coordinates
            search_bbox_size: Size of search bounding box
            search_bbox_radius: Radius for search (unused in this implementation)
            n_max_dir_to_keep: Maximum number of directions to keep
            min_value: Minimum value to consider
            sort_crescent_order: Whether to sort in ascending order
            
        Returns:
            List of (direction, score) tuples
        """
        cx, cy, cz = center_pixel
        directions_and_scores = []
        
        # Define bounding box
        bbox_min = (max(0, cx - search_bbox_size), 
                   max(0, cy - search_bbox_size), 
                   max(0, cz - search_bbox_size))
        bbox_max = (min(self.grid_data.shape[0], cx + search_bbox_size + 1),
                   min(self.grid_data.shape[1], cy + search_bbox_size + 1),
                   min(self.grid_data.shape[2], cz + search_bbox_size + 1))
        
        # Search within bounding box
        for i in range(bbox_min[0], bbox_max[0]):
            for j in range(bbox_min[1], bbox_max[1]):
                for k in range(bbox_min[2], bbox_max[2]):
                    
                    if (i, j, k) == center_pixel:
                        continue
                        
                    pixel_value = self.grid_data[i, j, k]
                    
                    if pixel_value > min_value:
                        # Convert to cartesian coordinates
                        curr_point = self._pixel_to_cartesian(i, j, k)
                        
                        # Calculate potential growth direction
                        potential_dir = curr_point - reference_point
                        dir_norm_3d = np.linalg.norm(potential_dir)
                        
                        # Weight the direction
                        weighted_direction = potential_dir * pixel_value * (
                            (search_bbox_size * self.resolution) - dir_norm_3d
                        )
                        
                        directions_and_scores.append((weighted_direction, pixel_value))
        
        # Sort if requested
        if sort_crescent_order:
            directions_and_scores.sort(key=lambda x: x[1])
        
        # Limit number of directions
        if len(directions_and_scores) > n_max_dir_to_keep:
            directions_and_scores = directions_and_scores[:n_max_dir_to_keep]
            
        return directions_and_scores
    
    def grow(self, n_iter_max: int = 100) -> int:
        """
        Grow the contour by adding points at head and tail
        
        Args:
            n_iter_max: Maximum number of growth iterations
            
        Returns:
            Number of iterations performed
        """
        self._update_points()
        
        sigma_l1_head = float('inf')
        sigma_l1_back = float('inf')
        has_grown_head = True
        has_grown_back = True
        grow_dir_head = np.zeros(3)
        grow_dir_back = np.zeros(3)
        
        for i in range(n_iter_max):
            if (sigma_l1_head <= self.seuil_sigma_l1 and 
                sigma_l1_back <= self.seuil_sigma_l1) or \
               (not has_grown_head and not has_grown_back):
                break
                
            # Get growing directions
            self._get_growing_directions(
                self.cone_angle_max_degrees,
                self.cone_size_max_pixels,
                self.seuil_sigma_l1,
                sigma_l1_head, sigma_l1_back,
                has_grown_head, has_grown_back,
                grow_dir_head, grow_dir_back
            )
            
            # Update growth flags
            if sigma_l1_head < self.seuil_sigma_l1:
                has_grown_head = False
            if sigma_l1_back < self.seuil_sigma_l1:
                has_grown_back = False
            
            # Add new points
            if has_grown_head and has_grown_back:
                # Add points at both ends
                new_head = self.points[0].to_array() - self.grow_coeff * grow_dir_head
                new_tail = self.points[-1].to_array() - self.grow_coeff * grow_dir_back
                
                self.points.insert(0, ContourPoint.from_array(new_head))
                self.points.append(ContourPoint.from_array(new_tail))
                
            elif has_grown_head:
                # Add point at head
                new_head = self.points[0].to_array() - self.grow_coeff * grow_dir_head
                self.points.insert(0, ContourPoint.from_array(new_head))
                
            elif has_grown_back:
                # Add point at tail
                new_tail = self.points[-1].to_array() - self.grow_coeff * grow_dir_back
                self.points.append(ContourPoint.from_array(new_tail))
            
            self._update_points()
        
        return i
    
    def _get_growing_directions(self, cone_angle_degrees: float, cone_size_pixels: int,
                               seuil_sigma_l1: float, sigma_l1_head: float, sigma_l1_back: float,
                               has_grown_head: bool, has_grown_back: bool,
                               grow_dir_head: np.ndarray, grow_dir_back: np.ndarray) -> None:
        """
        Calculate growing directions at head and tail
        
        Args:
            cone_angle_degrees: Maximum cone angle for growth
            cone_size_pixels: Size of search cone in pixels
            seuil_sigma_l1: Reliability threshold
            sigma_l1_head: Head reliability measure (output)
            sigma_l1_back: Back reliability measure (output)
            has_grown_head: Whether head grew last iteration (output)
            has_grown_back: Whether back grew last iteration (output)
            grow_dir_head: Head growth direction (output)
            grow_dir_back: Back growth direction (output)
        """
        # Reset directions
        grow_dir_head[:] = 0
        grow_dir_back[:] = 0
        
        # Calculate head growth if needed
        if sigma_l1_head > seuil_sigma_l1 and has_grown_head:
            head_point = self.points[0].to_array()
            head_pixel = self._cartesian_to_pixel(head_point)
            head_tangent = self._get_tangent_at_head(True)
            head_tangent *= -1  # Reverse for head
            
            grow_dir_head[:] = self._get_gradient_energy_grow_proportional_resolution(
                head_point, head_pixel, head_tangent, cone_angle_degrees,
                cone_size_pixels, sigma_l1_head, has_grown_head
            )
        
        # Calculate back growth if needed
        if sigma_l1_back > seuil_sigma_l1 and has_grown_back:
            tail_point = self.points[-1].to_array()
            tail_pixel = self._cartesian_to_pixel(tail_point)
            tail_tangent = self._get_tangent_at_tail(True)
            
            grow_dir_back[:] = self._get_gradient_energy_grow_proportional_resolution(
                tail_point, tail_pixel, tail_tangent, cone_angle_degrees,
                cone_size_pixels, sigma_l1_back, has_grown_back
            )
    
    def _get_gradient_energy_grow_proportional_resolution(self, snake_extremity_point: np.ndarray,
                                                        snake_extremity_pix: Tuple[int, int, int],
                                                        snake_tangent_at_extremity: np.ndarray,
                                                        cone_angle_degrees: float,
                                                        cone_size_pixels: int,
                                                        out_sigma_l1: float,
                                                        out_has_grown: bool) -> np.ndarray:
        """
        Calculate gradient energy for growth with proportional resolution
        
        Args:
            snake_extremity_point: Point at snake extremity
            snake_extremity_pix: Pixel coordinates of extremity
            snake_tangent_at_extremity: Tangent vector at extremity
            cone_angle_degrees: Maximum cone angle
            cone_size_pixels: Size of search cone
            out_sigma_l1: Output reliability measure
            out_has_grown: Output growth flag
            
        Returns:
            Growth direction vector
        """
        # Get directions and scores within cone
        directions_and_scores = self._get_directions_and_scores_in_cone(
            snake_extremity_pix, snake_extremity_point, snake_tangent_at_extremity,
            cone_angle_degrees, cone_size_pixels, 0, 40, self.min_value, True
        )
        
        if len(directions_and_scores) == 0:
            out_sigma_l1 = 0.0
            out_has_grown = False
            return np.zeros(3)
        
        # Calculate weighted average direction
        directions = np.array([d for d, _ in directions_and_scores])
        weights = np.array([s for _, s in directions_and_scores])
        
        weighted_direction = np.average(directions, weights=weights, axis=0)
        
        # Calculate reliability measure (simplified)
        total_weight = np.sum(weights)
        out_sigma_l1 = total_weight / len(weights) if len(weights) > 0 else 0.0
        out_has_grown = total_weight > 0
        
        return weighted_direction
    
    def _get_directions_and_scores_in_cone(self, center_pixel: Tuple[int, int, int],
                                         reference_point: np.ndarray,
                                         cone_direction: np.ndarray,
                                         cone_angle_deg: float,
                                         search_bbox_size: int, search_bbox_radius: int,
                                         n_max_dir_to_keep: int, min_value: float,
                                         sort_crescent_order: bool) -> List[Tuple[np.ndarray, float]]:
        """
        Get potential growth directions within a cone
        
        Args:
            center_pixel: Center pixel coordinates
            reference_point: Reference point in cartesian coordinates
            cone_direction: Direction of cone
            cone_angle_deg: Maximum angle of cone
            search_bbox_size: Size of search bounding box
            search_bbox_radius: Radius for search (unused)
            n_max_dir_to_keep: Maximum number of directions to keep
            min_value: Minimum value to consider
            sort_crescent_order: Whether to sort in ascending order
            
        Returns:
            List of (direction, score) tuples
        """
        cx, cy, cz = center_pixel
        cone_direction_spatial = cone_direction / np.linalg.norm(cone_direction)
        directions_and_scores = []
        
        # Define bounding box
        bbox_min = (max(0, cx - search_bbox_size), 
                   max(0, cy - search_bbox_size), 
                   max(0, cz - search_bbox_size))
        bbox_max = (min(self.grid_data.shape[0], cx + search_bbox_size + 1),
                   min(self.grid_data.shape[1], cy + search_bbox_size + 1),
                   min(self.grid_data.shape[2], cz + search_bbox_size + 1))
        
        # Search within bounding box
        for i in range(bbox_min[0], bbox_max[0]):
            for j in range(bbox_min[1], bbox_max[1]):
                for k in range(bbox_min[2], bbox_max[2]):
                    
                    if (i, j, k) == center_pixel:
                        continue
                        
                    pixel_value = self.grid_data[i, j, k]
                    
                    if pixel_value > min_value:
                        # Convert to cartesian coordinates
                        curr_point = self._pixel_to_cartesian(i, j, k)
                        
                        # Calculate potential growth direction
                        potential_dir = curr_point - reference_point
                        potential_dir_spatial = potential_dir
                        potential_dir_norm_3d = np.linalg.norm(potential_dir_spatial)
                        
                        if potential_dir_norm_3d > 0:
                            # Test angle in 3D
                            scalar_product = np.dot(
                                potential_dir_spatial / potential_dir_norm_3d,
                                cone_direction_spatial
                            )
                            
                            # Check if vectors are in same direction
                            if scalar_product >= 0:
                                angle = np.arccos(scalar_product) * 180.0 / np.pi
                                if angle <= cone_angle_deg:
                                    weighted_direction = potential_dir * pixel_value
                                    directions_and_scores.append((weighted_direction, pixel_value))
        
        # Sort if requested
        if sort_crescent_order:
            directions_and_scores.sort(key=lambda x: x[1])
        
        # Limit number of directions
        if len(directions_and_scores) > n_max_dir_to_keep:
            directions_and_scores = directions_and_scores[:n_max_dir_to_keep]
            
        return directions_and_scores
    
    def relax(self, n_iter_max: int = 1000) -> int:
        """
        Relax the contour to optimize its shape
        
        Args:
            n_iter_max: Maximum number of relaxation iterations
            
        Returns:
            Number of iterations performed
        """
        if len(self.points) < 3:
            return 0
        
        # Get geometric matrix inverse
        a_plus_rho_i_inverse = self._get_a_plus_rho_i_inverse(
            self.alpha, self.beta, self.time_step
        )
        
        n_pts = len(self.points)
        rho = 1.0 / self.time_step
        average_movement = float('inf')
        
        for relax_step in range(n_iter_max):
            if average_movement <= self.thresh_average_movement:
                break
                
            # Save current points
            save_points = [p.to_array() for p in self.points]
            
            # Calculate energy gradients
            gradient_energy_image = self._get_gradient_energy_orthogonal_image_matrix(
                self.global_weight
            )
            
            second_diff_on_tangent = self._get_second_differential_on_tangent_divided_by_tangent_norm_multiplied_by_image_energy(
                self.global_weight
            )
            
            # Calculate next iteration positions
            current_points = np.array([p.to_array() for p in self.points])
            points_next_iteration = a_plus_rho_i_inverse @ (
                (rho * current_points) - 
                (self.gamma * gradient_energy_image) + 
                second_diff_on_tangent
            )
            
            # Update points and calculate movement
            average_movement = self._update_points_and_get_average_movement(points_next_iteration)
        
        return relax_step
    
    def _get_a_plus_rho_i_inverse(self, alpha: float, beta: float, time_step: float) -> np.ndarray:
        """
        Calculate (A + ρI)^(-1) matrix for relaxation
        
        Args:
            alpha: Elasticity parameter
            beta: Bending stiffness parameter
            time_step: Time step for relaxation
            
        Returns:
            Inverse matrix
        """
        n_pts = len(self.points)
        rho = 1.0 / time_step
        
        # Identity matrix
        identity = np.eye(n_pts)
        
        # Geometric matrix A
        matrix_a = self._get_geometric_matrix(alpha, beta)
        
        # Calculate inverse
        a_plus_rho_i = matrix_a + (rho * identity)
        
        return inv(a_plus_rho_i)
    
    def _get_geometric_matrix(self, alpha: float, beta: float) -> np.ndarray:
        """
        Calculate geometric matrix A
        
        Args:
            alpha: Elasticity parameter
            beta: Bending stiffness parameter
            
        Returns:
            Geometric matrix
        """
        n_pts = len(self.points)
        matrix_a = np.zeros((n_pts, n_pts))
        
        # Fill matrix (excluding boundary conditions)
        for i in range(n_pts):
            if 2 <= i < n_pts - 2:
                if i - 2 >= 0:
                    matrix_a[i, i-2] = beta
                if i - 1 >= 0:
                    matrix_a[i, i-1] = -alpha - 3*beta
                
                matrix_a[i, i] = 2*alpha + 6*beta
                
                if i + 1 <= n_pts - 1:
                    matrix_a[i, i+1] = -alpha - 3*beta
                if i + 2 <= n_pts - 1:
                    matrix_a[i, i+2] = beta
        
        return matrix_a
    
    def _get_gradient_energy_orthogonal_image_matrix(self, global_weight: float) -> np.ndarray:
        """
        Calculate gradient energy orthogonal to image matrix
        
        Args:
            global_weight: Global weight for energy calculation
            
        Returns:
            Gradient energy matrix
        """
        n_pts = len(self.points)
        result = np.zeros((n_pts, 3))
        
        for i in range(n_pts):
            # Calculate gradient energy at point
            curr_grad = self._gradient_energy_image_interpolated_at_point(i, global_weight)
            
            # Calculate contour direction
            direction_contour = self._direction_contours_at_point(i)
            direction_contour_norm = np.linalg.norm(direction_contour)
            
            if direction_contour_norm > 0:
                direction_contour_normalized = direction_contour / direction_contour_norm
                
                # Calculate orthogonal component
                curr_grad = curr_grad - (direction_contour_normalized * 
                                       np.dot(curr_grad, direction_contour_normalized))
            
            result[i, :] = curr_grad
        
        return result
    
    def _get_second_differential_on_tangent_divided_by_tangent_norm_multiplied_by_image_energy(
        self, global_weight: float) -> np.ndarray:
        """
        Calculate second differential energy matrix
        
        Args:
            global_weight: Global weight for energy calculation
            
        Returns:
            Second differential energy matrix
        """
        n_pts = len(self.points)
        result = np.zeros((n_pts, 3))
        
        for i in range(n_pts):
            point = self.points[i].to_array()
            cur_tangent = self._get_tangent_at_point(i, True)
            cur_second_diff = self._get_second_differential_at_point(i)
            tangent_norm = np.linalg.norm(cur_tangent)
            
            if tangent_norm > 0:
                scal_prod = np.dot(cur_tangent, cur_second_diff)
                scal_prod_normalized = scal_prod / (tangent_norm * tangent_norm)
                
                pixel_at_i = self._cartesian_to_pixel(point)
                energy_at_pixel = self._energy_image_at_pixel(pixel_at_i, global_weight)
                
                for c in range(3):
                    result[i, c] = ((cur_second_diff[c] - 
                                   (scal_prod_normalized * cur_tangent[c])) / 
                                   tangent_norm) * energy_at_pixel
        
        return result
    
    def _update_points_and_get_average_movement(self, new_points: np.ndarray) -> float:
        """
        Update points and calculate average movement
        
        Args:
            new_points: New point positions
            
        Returns:
            Average movement distance
        """
        n_pts = len(self.points)
        movement_sum = 0.0
        
        for i in range(n_pts):
            old_point = self.points[i].to_array()
            new_point = new_points[i, :]
            
            # Calculate movement
            movement = np.linalg.norm(old_point - new_point)
            movement_sum += movement
            
            # Update point position
            self.points[i] = ContourPoint.from_array(new_point)
        
        return movement_sum / n_pts if n_pts > 0 else 0.0
    
    def resample(self, sample_res: float) -> None:
        """
        Resample the contour to maintain consistent point spacing
        
        Args:
            sample_res: Target sampling resolution
        """
        if len(self.points) < 2:
            return
            
        total_length = self.length_3d()
        
        if total_length <= sample_res:
            return
        
        # Calculate number of points needed
        n_points_resampled = int(np.floor(total_length / sample_res)) + 1
        
        if n_points_resampled < 2:
            return
        
        # Create new points array
        old_points = [p.to_array() for p in self.points]
        new_points = []
        
        # Add first point
        new_points.append(ContourPoint.from_array(old_points[0]))
        
        # Calculate cumulative distances
        cumulative_distances = [0.0]
        for i in range(1, len(old_points)):
            dist = np.linalg.norm(old_points[i] - old_points[i-1])
            cumulative_distances.append(cumulative_distances[-1] + dist)
        
        # Add intermediate points
        for i in range(1, n_points_resampled - 1):
            target_distance = i * sample_res
            
            # Find segment containing target distance
            segment_idx = 0
            for j in range(len(cumulative_distances) - 1):
                if (cumulative_distances[j] <= target_distance < 
                    cumulative_distances[j + 1]):
                    segment_idx = j
                    break
            
            # Interpolate point
            if segment_idx < len(old_points) - 1:
                p1 = old_points[segment_idx]
                p2 = old_points[segment_idx + 1]
                d1 = cumulative_distances[segment_idx]
                d2 = cumulative_distances[segment_idx + 1]
                
                if d2 > d1:
                    t = (target_distance - d1) / (d2 - d1)
                    interpolated_point = p1 + t * (p2 - p1)
                    new_points.append(ContourPoint.from_array(interpolated_point))
        
        # Add last point
        new_points.append(ContourPoint.from_array(old_points[-1]))
        
        self.points = new_points
    
    def length_3d(self) -> float:
        """
        Calculate the 3D length of the contour
        
        Returns:
            Total length of the contour
        """
        if len(self.points) < 2:
            return 0.0
        
        total_length = 0.0
        for i in range(len(self.points) - 1):
            p1 = self.points[i].to_array()
            p2 = self.points[i + 1].to_array()
            total_length += np.linalg.norm(p2 - p1)
        
        return total_length
    
    def _get_tangent_at_point(self, i: int, normalize: bool = True) -> np.ndarray:
        """
        Calculate tangent vector at a point
        
        Args:
            i: Point index
            normalize: Whether to normalize the tangent
            
        Returns:
            Tangent vector
        """
        if i == 0:
            # Forward difference
            curr = self.points[i].to_array()
            next_point = self.points[i + 1].to_array()
            tangent = next_point - curr
        elif i == len(self.points) - 1:
            # Backward difference
            curr = self.points[i].to_array()
            prev = self.points[i - 1].to_array()
            tangent = curr - prev
        else:
            # Central difference
            prev = self.points[i - 1].to_array()
            next_point = self.points[i + 1].to_array()
            tangent = next_point - prev
        
        if normalize:
            norm = np.linalg.norm(tangent)
            if norm > 0:
                tangent = tangent / norm
        
        return tangent
    
    def _get_tangent_at_head(self, normalize: bool = True) -> np.ndarray:
        """Get tangent at head (first point)"""
        return self._get_tangent_at_point(0, normalize)
    
    def _get_tangent_at_tail(self, normalize: bool = True) -> np.ndarray:
        """Get tangent at tail (last point)"""
        return self._get_tangent_at_point(len(self.points) - 1, normalize)
    
    def _get_second_differential_at_point(self, i: int) -> np.ndarray:
        """
        Calculate second differential at a point
        
        Args:
            i: Point index
            
        Returns:
            Second differential vector
        """
        if i == 0 or i == len(self.points) - 1:
            return np.zeros(3)
        
        prev = self.points[i - 1].to_array()
        curr = self.points[i].to_array()
        next_point = self.points[i + 1].to_array()
        
        return prev - 2 * curr + next_point
    
    def _direction_contours_at_point(self, i: int) -> np.ndarray:
        """
        Calculate contour direction at a point
        
        Args:
            i: Point index
            
        Returns:
            Contour direction vector
        """
        if i == 0:
            return self._get_tangent_at_head(True)
        elif i == len(self.points) - 1:
            return self._get_tangent_at_tail(True)
        else:
            prev = self.points[i - 1].to_array()
            next_point = self.points[i + 1].to_array()
            return next_point - prev
    
    def _gradient_energy_image_interpolated_at_point(self, i: int, global_weight: float) -> np.ndarray:
        """
        Calculate interpolated gradient energy at a point
        
        Args:
            i: Point index
            global_weight: Global weight for energy calculation
            
        Returns:
            Gradient energy vector
        """
        curr_point = self.points[i].to_array()
        curr_pixel = self._cartesian_to_pixel(curr_point)
        curr_pixel_center = self._pixel_to_cartesian(*curr_pixel)
        
        grad_interpolated = np.zeros(3)
        
        for coord in range(3):
            # Get gradient at current pixel
            curr_pixel_grad = self._partial_derivative_energy_image_at_pixel(
                curr_pixel, global_weight, coord
            )
            
            # Interpolate based on position within pixel
            if curr_point[coord] < curr_pixel_center[coord]:
                # Interpolate with previous pixel
                prev_pixel = list(curr_pixel)
                prev_pixel[coord] -= 1
                prev_pixel = tuple(prev_pixel)
                
                if self._is_valid_pixel(*prev_pixel):
                    prev_pixel_center = self._pixel_to_cartesian(*prev_pixel)
                    prev_pixel_grad = self._partial_derivative_energy_image_at_pixel(
                        prev_pixel, global_weight, coord
                    )
                    
                    # Linear interpolation
                    t = (curr_point[coord] - prev_pixel_center[coord]) / (
                        curr_pixel_center[coord] - prev_pixel_center[coord]
                    )
                    grad_interpolated[coord] = prev_pixel_grad + t * (curr_pixel_grad - prev_pixel_grad)
                else:
                    grad_interpolated[coord] = curr_pixel_grad
            else:
                # Interpolate with next pixel
                next_pixel = list(curr_pixel)
                next_pixel[coord] += 1
                next_pixel = tuple(next_pixel)
                
                if self._is_valid_pixel(*next_pixel):
                    next_pixel_center = self._pixel_to_cartesian(*next_pixel)
                    next_pixel_grad = self._partial_derivative_energy_image_at_pixel(
                        next_pixel, global_weight, coord
                    )
                    
                    # Linear interpolation
                    t = (curr_point[coord] - curr_pixel_center[coord]) / (
                        next_pixel_center[coord] - curr_pixel_center[coord]
                    )
                    grad_interpolated[coord] = curr_pixel_grad + t * (next_pixel_grad - curr_pixel_grad)
                else:
                    grad_interpolated[coord] = curr_pixel_grad
        
        return grad_interpolated
    
    def _partial_derivative_energy_image_at_pixel(self, p: Tuple[int, int, int], 
                                                global_weight: float, coord: int) -> float:
        """
        Calculate partial derivative of energy image at pixel
        
        Args:
            p: Pixel coordinates
            global_weight: Global weight for energy calculation
            coord: Coordinate to derive (0, 1, or 2)
            
        Returns:
            Partial derivative value
        """
        # Get previous and next pixels in direction
        prev_pixel = list(p)
        prev_pixel[coord] -= 1
        prev_pixel = tuple(prev_pixel)
        
        next_pixel = list(p)
        next_pixel[coord] += 1
        next_pixel = tuple(next_pixel)
        
        is_prev_in_image = self._is_valid_pixel(*prev_pixel)
        is_next_in_image = self._is_valid_pixel(*next_pixel)
        
        if not is_prev_in_image and not is_next_in_image:
            return 0.0
        elif is_prev_in_image and is_next_in_image:
            # Central difference
            energy_prev = self._energy_image_at_pixel(prev_pixel, global_weight)
            energy_next = self._energy_image_at_pixel(next_pixel, global_weight)
            return (energy_next - energy_prev) / 2.0
        elif is_prev_in_image and not is_next_in_image:
            # Forward difference
            energy_prev = self._energy_image_at_pixel(prev_pixel, global_weight)
            energy_curr = self._energy_image_at_pixel(p, global_weight)
            return energy_curr - energy_prev
        else:
            # Backward difference
            energy_curr = self._energy_image_at_pixel(p, global_weight)
            energy_next = self._energy_image_at_pixel(next_pixel, global_weight)
            return energy_next - energy_curr
    
    def _energy_image_at_pixel(self, p: Tuple[int, int, int], global_weight: float) -> float:
        """
        Calculate energy image value at pixel
        
        Args:
            p: Pixel coordinates
            global_weight: Global weight for energy calculation
            
        Returns:
            Energy value
        """
        return (global_weight * self._energy_image_global_at_pixel(p) + 
                (1.0 - global_weight) * self._energy_image_local_at_pixel(p))
    
    def _energy_image_global_at_pixel(self, p: Tuple[int, int, int]) -> float:
        """Calculate global energy image value at pixel"""
        if not self._is_valid_pixel(*p):
            return 0.0
        
        pixel_value = self.grid_data[p]
        data_min = np.min(self.grid_data[self.grid_data > 0])
        data_max = np.max(self.grid_data)
        
        if data_max <= data_min:
            return 0.0
        
        return (data_min - pixel_value) / (data_max - data_min)
    
    def _energy_image_local_at_pixel(self, p: Tuple[int, int, int]) -> float:
        """Calculate local energy image value at pixel (simplified)"""
        if not self._is_valid_pixel(*p):
            return 0.0
        
        # Simplified local energy - could be enhanced with gradient magnitude
        return 0.0
    
    def _update_points(self) -> None:
        """Update points to ensure they stay within grid bounds"""
        bbox_bot = np.array([0, 0, 0])
        bbox_top = np.array(self.grid_data.shape) * self.resolution
        
        for point in self.points:
            for coord in range(3):
                if point.to_array()[coord] < bbox_bot[coord]:
                    point.to_array()[coord] = bbox_bot[coord] + 0.00001
                if point.to_array()[coord] > bbox_top[coord]:
                    point.to_array()[coord] = bbox_top[coord] - 0.00001
    
    def _pixel_to_cartesian(self, x: int, y: int, z: int) -> np.ndarray:
        """Convert pixel coordinates to cartesian coordinates"""
        return np.array([x, y, z]) * self.resolution
    
    def _cartesian_to_pixel(self, point: np.ndarray) -> Tuple[int, int, int]:
        """Convert cartesian coordinates to pixel coordinates"""
        return tuple(np.floor(point / self.resolution).astype(int))
    
    def _is_valid_pixel(self, x: int, y: int, z: int) -> bool:
        """Check if pixel coordinates are valid"""
        return (0 <= x < self.grid_data.shape[0] and
                0 <= y < self.grid_data.shape[1] and
                0 <= z < self.grid_data.shape[2])
    
    def get_points_as_array(self) -> np.ndarray:
        """Get contour points as numpy array"""
        return np.array([p.to_array() for p in self.points])
    
    def get_n_points(self) -> int:
        """Get number of points in contour"""
        return len(self.points)


class ActiveContoursWorkflow:
    """
    High-level workflow for extracting characteristic curves using active contours
    """
    
    def __init__(self, grid_data: np.ndarray, resolution: float = 1.0):
        """
        Initialize the active contours workflow
        
        Args:
            grid_data: 3D numpy array representing the grid
            resolution: Grid resolution in units per pixel
        """
        self.grid_data = grid_data
        self.resolution = resolution
        
        # Workflow parameters
        self.min_value_for_maxima = 10
        self.maxima_neighborhood_size = 1
        self.min_height_for_maximum_search = 1.1
        self.max_height_for_maximum_search = 2.0
        self.min_value_for_maximum_search = 5
        self.tree_height_maximum = 30.0
        self.grow_coeff = 1.0
        self.angle_cone_recherche = 20.0
        self.taille_cone_recherche_cm = 10.0
        self.seuil_sigma_l1 = 0.75
        self.n_iter_max_optim = 1000
        self.time_step = 0.001
        self.longueur_min = 1.0
        self.thresh_grad_move = 0.001
        self.n_snakes_max = 1
        
    def extract_curves(self) -> List[OpenActiveContours]:
        """
        Extract characteristic curves from the grid
        
        Returns:
            List of OpenActiveContours objects representing the extracted curves
        """
        # Calculate cone size in pixels
        taille_cone_recherche = int(np.ceil((self.taille_cone_recherche_cm / 100.0) / self.resolution))
        n_iter_max_grow = max(int(self.tree_height_maximum / self.resolution), 1)
        
        # Get local maxima
        local_maximas = self._get_local_maximas()
        
        if not local_maximas:
            warnings.warn("No local maxima found in the grid")
            return []
        
        # Create repulsion grid
        repulsion_grid = np.zeros_like(self.grid_data, dtype=bool)
        
        # Extract curves
        curves = []
        for local_maxima in local_maximas:
            if len(curves) >= self.n_snakes_max:
                break
                
            if not repulsion_grid[local_maxima]:
                # Initialize active contour
                contour = OpenActiveContours(
                    self.grid_data, self.resolution,
                    self.min_value_for_maximum_search, 20
                )
                
                if contour.initialize_from_point(local_maxima):
                    # Grow the contour
                    contour.grow(n_iter_max_grow)
                    
                    # Check if contour is long enough
                    if contour.length_3d() > self.longueur_min:
                        # Resample contour
                        contour.resample(0.1)
                        
                        # Mark as repulsive
                        self._mark_repulsion(contour, repulsion_grid, 1.5)
                        
                        # Relax/optimize contour
                        contour.relax(self.n_iter_max_optim)
                        
                        # Final resample
                        contour.resample(0.1)
                        
                        curves.append(contour)
        
        return curves
    
    def _get_local_maximas(self) -> List[Tuple[int, int, int]]:
        """
        Get local maxima in the grid
        
        Returns:
            List of local maxima coordinates
        """
        from .workflows import GridWorkflows
        
        local_maxima_data = GridWorkflows.local_maxima_detection(
            self.grid_data,
            self.maxima_neighborhood_size,
            self.min_value_for_maxima,
            True
        )
        
        # Filter by height constraints
        filtered_maximas = []
        for x, y, z, value in local_maxima_data:
            height = z * self.resolution
            if (self.min_height_for_maximum_search <= height <= 
                self.max_height_for_maximum_search):
                filtered_maximas.append((x, y, z))
        
        return filtered_maximas
    
    def _mark_repulsion(self, contour: OpenActiveContours, 
                       repulsion_grid: np.ndarray, repulse_factor: float) -> None:
        """
        Mark contour area as repulsive in the repulsion grid
        
        Args:
            contour: Active contour to mark
            repulsion_grid: Repulsion grid to update
            repulse_factor: Repulsion factor
        """
        # Simplified repulsion marking
        # In a full implementation, this would mark the area around each contour point
        points = contour.get_points_as_array()
        
        for point in points:
            pixel = contour._cartesian_to_pixel(point)
            if contour._is_valid_pixel(*pixel):
                # Mark surrounding area
                for dx in range(-2, 3):
                    for dy in range(-2, 3):
                        for dz in range(-2, 3):
                            nx, ny, nz = pixel[0] + dx, pixel[1] + dy, pixel[2] + dz
                            if contour._is_valid_pixel(nx, ny, nz):
                                repulsion_grid[nx, ny, nz] = True 