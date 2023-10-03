import typing
from typing import List, Tuple

import numpy as np
from skued.simulation import structure_factor
import celltools.cell.contents as cc
from celltools.cell.generate import cell_to_crystal



def diffraction_from_cell(hkl: List[Tuple[int, int ,int]], cell: cc.cell) -> Tuple[List[float, ...], List[float, ...]]:
    """

    Parameters
    ----------
    hkl: iterable of (3,)-tuples
        list of hkl indices to simulate the structure factor from
    cell: :class:`cell`
        cell to simulate the diffraction from

    Returns
    -------

    """
    crystal = cell_to_crystal(cell)
    _scatt_vector = []
    _structure_factor = []
    for _hkl in hkl:
        _scatt_vector.append(crystal.scattering_vector(_hkl))
        _structure_factor.append((structure_factor(crystal, _hkl)))
    return _scatt_vector, _structure_factor