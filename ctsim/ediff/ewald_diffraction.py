import typing
from typing import List, Tuple, Union, Literal
import numpy as np
from itertools import product


from celltools.linalg import Vector
from celltools import Cell


def calculate_structure_factor(cell: Cell, k_point: Vector) -> complex:
    """
    calculates structure factor for 
    Parameters
    ----------
    cell : Cell
        _description_
    k_point : Vector
        _description_

    Returns
    -------
    complex
        _description_
    """
    pass

def ewald_diffraction(k_vector: Vector, cell: Cell, grid: (int, int, int), *,
                      ds_cutoff: float = 1e-6, 
                      distance_correction: Union[None, Literal["exponential"], Literal["gaussian"]] = None) -> [List[Vector], List[complex], List[float]]:
    """_summary_

    Parameters
    ----------
    k_vector : Vector
        k-vector of incoming wave. This defines the radius of the Ewald sphere. Units are expected to be the same as the reciprocal lattice of the input cell
    cell : Cell
        input unit cell
    grid : (int, int, int)
        defines the size of the reciprocal grid through which the cut with the Ewald sphere is performed. Input is a tuple (nx, ny, nz) which defines an isotropic grid size around the reciprocal space origin (0,0,0) with +-ni points for each direction. For example (10,10,10) gives (-10, 0, 0), (-9,0,0), ... ,(0,0,0), ..., (10,0,0); (-10,-10,0),(-9,-10,0),... and so on.  
    ds_cutoff : float, optional
        In reality, the scattering condition is not only fulfilled if the Gamma point lies on the Ewald sphere but in a region around the Gamma point due to mosaicity, crystal size effects etc. This value gives the cutoff ratio in which range the scattering condition is viewed as fulfilled, i.e. if |(distance to Ewald sphere)/(radius of Ewald sphere)| < ds_cutoff reciprocal lattice point fulfills scattering condition. by default 1e-6
    distance_correction : None, "exponential", "gaussian", optional
        Describes the intensity modulation for points fulfilling the scattering condition with distance from the Ewald sphere. ONLY None IMPLEMENTED so far. by default None

    Returns
    -------
    refelctions, structure_factor, ds
    
    reflections: List[Vector]
        list of reciprocal lattice vectors fulfilling the scattering condition
    
    structure_factor: List[complex]
        list of structure factors of points fulfilling the scattering condition

    ds: List[float]
        if distance_correction = None
            distances to Ewald sphere of points fulfilling the scattering condition
        else:
            weights of scattering intensities modulated with distance from the Ewald sphere according to distance_correction
    """

    # make k-grid
    _nx, _ny, _nz = [-grid[0] + ix for ix in range(2*grid[0]+1)],  [-grid[1] + iy for iy in range(2*grid[1]+1)],  [-grid[2] + iz for iz in range(2*grid[2]+1)]

    k_grid = [Vector(hkl, cell.reciprocal_lattice) for hkl in product(_nx,_ny,_nz)]
    k_origin = Vector((0,0,0), cell.reciprocal_lattice)
    # find center of Ewald sphere 
    sphere_center = -1*k_vector
    # build distance check sphere center and k grid
    def distance_to_sphere_check(k_point):
        distance_to_sphere = (k_point - sphere_center).abs_global
        if distance_to_sphere / k_vector.abs_global < ds_cutoff:
            return True, distance_to_sphere
        else:
            return False, distance_to_sphere

    # run distance check over k-grid 
    reflections, structure_factor, ds = [], [], []

    for k_point in k_grid:
        condition_fulfilled, distance_to_sphere = distance_to_sphere_check(k_point)
        if condition_fulfilled:
            reflections.append(k_point)
            structure_factor.append(calculate_structure_factor(k_point))
            ds.append(distance_to_sphere)
    
    return reflections, structure_factor, ds
