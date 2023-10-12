import typing
from typing import List, Tuple

import numpy as np
from skued.simulation import structure_factor

from . import IndexLike
from celltools import Cell, SuperCell
from celltools.cell.generate import cell_to_crystal

def diffraction_from_cell(hkl: List[IndexLike], cell: Cell) -> Tuple[List[np.ndarray], List[complex]]:
    """
    calculates the scattering vectors and respective structure factor as the complex scattering amplitude for a list of
    hkl indices from a unit cell. Uses structure_factor method from scikit-ued to calculate the structure factor.
    Parameters
    ----------
    hkl: iterable of (3,)-tuples of length N
        list of hkl indices to simulate the structure factor from
    cell: :class:`Cell`
        cell to simulate the diffraction from

    Returns
    -------
        scatt_vector: list of nd.arrays size (3,) of length N
            scattering vectors q of respective hkl index in reciprocal space
        structure_factor: list of complex numbers of length N
            structure_factor of reflection of respective hkl index
    """
    crystal = cell_to_crystal(cell)
    _scatt_vector = []
    _structure_factor = []
    for _hkl in hkl:
        _scatt_vector.append(crystal.scattering_vector(_hkl))
        _structure_factor.append((structure_factor(crystal, *_hkl)))
    return _scatt_vector, _structure_factor
