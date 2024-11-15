import typing
from typing import List, Tuple, Union, Literal
import numpy as np
from itertools import product
from skued import electron_wavelength as skued_electron_wavelength
from skued.simulation import affe


from celltools.linalg import Vector, Plane, BasisTransformation
from celltools import Cell


def calculate_structure_factor(cell: Cell, k_point: Vector) -> complex:
    """
    calculates structure factor for reflection at reciprocal lattice point given by K_point at 0K (DW NOT INCLUDED).

    Parameters
    ----------
    cell : Cell
        unit cell
    k_point : Vector
        reciprocal lattice vector

    Returns
    -------
    complex
        complex structure factor. Scattering intensity is proportional to absolute square.
    """
    cell_atoms = []
    for atm in cell.base[0]:
        cell_atoms.append(atm)
    for molc in cell.base[1]:
        for atm in molc.atoms:
            cell_atoms.append(atm)

    ff, r = [], []

    for atm in cell_atoms:
        ff.append(affe(atm.element, k_point.abs_global))
        r.append(atm.coords.global_coord)

    ff = np.array(ff)
    r = np.array(r)

    structure_factor = np.sum(
        ff * np.exp(-1j * np.sum(k_point.global_coord * r, axis=1))
    )

    ### Can basically be copied from diffraction_from_supercell

    return structure_factor


def ewald_diffraction(
    k_vector: Vector,
    cell: Cell,
    grid: (int, int, int),
    *,
    ds_cutoff: float = 1e-3,
    distance_correction: Union[
        None, Literal["exponential"], Literal["gaussian"]
    ] = None,
) -> [List[Vector], List[complex], List[float]]:
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
        In reality, the scattering condition is not only fulfilled if the Gamma point lies on the Ewald sphere but in a region around the Gamma point due to mosaicity, crystal size effects etc. This value gives the cutoff ratio in which range the scattering condition is viewed as fulfilled, i.e. if |(distance to Ewald sphere)/(radius of Ewald sphere)| < ds_cutoff reciprocal lattice point fulfills scattering condition. by default 1e-3
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
    _nx, _ny, _nz = (
        [-grid[0] + ix for ix in range(2 * grid[0] + 1)],
        [-grid[1] + iy for iy in range(2 * grid[1] + 1)],
        [-grid[2] + iz for iz in range(2 * grid[2] + 1)],
    )

    k_grid = [Vector(hkl, cell.reciprocal_lattice) for hkl in product(_nx, _ny, _nz)]
    k_origin = Vector((0, 0, 0), cell.reciprocal_lattice)
    # find center of Ewald sphere
    sphere_center = -1 * k_vector

    # build distance check sphere center and k grid
    def distance_to_sphere_check(k_point):
        distance_to_sphere = (k_point - sphere_center).abs_global - k_vector.abs_global
        if abs(distance_to_sphere / k_vector.abs_global) < ds_cutoff:
            return True, distance_to_sphere
        else:
            return False, distance_to_sphere

    # run distance check over k-grid
    reflections, structure_factor, ds = [], [], []

    for k_point in k_grid:
        condition_fulfilled, distance_to_sphere = distance_to_sphere_check(k_point)
        if condition_fulfilled:
            reflections.append(k_point)
            structure_factor.append(calculate_structure_factor(cell, k_point))
            ds.append(distance_to_sphere)

    return reflections, structure_factor, ds


class DiffractionExperiment:
    def __init__(
        self,
        cell: Cell,
        direction: Vector,
        electron_energy: float,
        grid: (int, int, int) = (5, 5, 5),
        *,
        ds_cutoff: float = 1e-3,
        distance_correction: Union[
            None, Literal["exponential"], Literal["gaussian"]
        ] = None,
    ):
        """
        This class represents a diffraction experiment and calculates the diffraction using the Ewald sphere method on projects the reciprocal lattice points fulfilling the scattering condition on a plane normal to the incoming electron beam as would be the case for a planar detector.
        NOTE: The projection does not take any distortions due to the curvature of the Ewald sphere into account but assumes an "Ewald plane".
        Parameters
        ----------
        cell : Cell
            input unit cell
        direction : Vector
            direction of the incoming electron beam in reciprocal lattice coordinates of cell
        electron_energy : float
            electron energy in kV
        grid : (int, int, int)
        defines the size of the reciprocal grid through which the cut with the Ewald sphere is performed. Input is a tuple (nx, ny, nz) which defines an isotropic grid size around the reciprocal space origin (0,0,0) with +-ni points for each direction. For example (10,10,10) gives (-10, 0, 0), (-9,0,0), ... ,(0,0,0), ..., (10,0,0); (-10,-10,0),(-9,-10,0),... and so on. by default (5,5,5)
        ds_cutoff : float, optional
            In reality, the scattering condition is not only fulfilled if the Gamma point lies on the Ewald sphere but in a region around the Gamma point due to mosaicity, crystal size effects etc. This value gives the cutoff ratio in which range the scattering condition is viewed as fulfilled, i.e. if |(distance to Ewald sphere)/(radius of Ewald sphere)| < ds_cutoff reciprocal lattice point fulfills scattering condition. by default 1e-3,
        distance_correction :  None, "exponential", "gaussian", optional
            Describes the intensity modulation for points fulfilling the scattering condition with distance from the Ewald sphere. ONLY None IMPLEMENTED so far. by default None

        Returns
        -------
        DiffractionExperiment
        """

        self._cell = cell
        self._direction = direction
        self._electron_energy = electron_energy
        self._grid = grid
        self._ewald_diffraction_kwargs = [ds_cutoff, distance_correction]

        self.hkl = []
        self.reflections = []
        self.structure_factor = None
        self.ds = None

        self._experiment_setup()

    def _experiment_setup(self):
        """
        set up diffraction experiment
        """
        # make k_vector incoming k-vector
        self._direction = (
            (2 * np.pi / skued_electron_wavelength(self.electron_energy))
            * (1 / self._direction.abs_global)
            * self._direction
        )

        # Ewald diffraction
        reflections, structure_factor, ds = ewald_diffraction(
            self.direction,
            self.cell,
            self.grid,
            ds_cutoff=self._ewald_diffraction_kwargs[0],
            distance_correction=self._ewald_diffraction_kwargs[1],
        )

        # make "Ewald plane"
        # incoming k-vector is normal vector
        self.diffraction_plane = Plane(
            Vector([0, 0, 0], self.cell.reciprocal_lattice), self.direction
        )

        # project Bragg reflections onto plane
        # save old coordinates as (hkl) labels

        self.structure_factor = structure_factor
        self.ds = ds

        to_plane_coordinates = BasisTransformation(
            self.cell.reciprocal_lattice, self.diffraction_plane.basis
        )

        for reflection in reflections:
            self.hkl.append(f"({reflection[0]} {reflection[1]} {reflection[2]})")
            _reflection_in_plane_coordinates = to_plane_coordinates.transform(
                reflection
            )
            _reflection_in_plane_coordinates = Vector(
                (
                    _reflection_in_plane_coordinates[0],
                    _reflection_in_plane_coordinates[1],
                    0,
                ),
                _reflection_in_plane_coordinates.basis,
            )
            self.reflections.append(_reflection_in_plane_coordinates)

    @property
    def cell(self):
        return self._cell

    @property
    def direction(self):
        return self._direction

    @property
    def electron_energy(self):
        return self._electron_energy

    @property
    def grid(self):
        return self._grid

    @direction.setter
    def direction(self, direction_new):
        self._direction = direction_new
        self._experiment_setup()

    @electron_energy.setter
    def electron_energy(self, electron_energy_new):
        self._electron_energy = electron_energy_new
        self._experiment_setup()

    @grid.setter
    def grid(self, grid_new):
        self._grid = grid_new
        self._experiment_setup()
