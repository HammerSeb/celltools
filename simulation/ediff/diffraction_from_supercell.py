import typing
from typing import List, Tuple

import numpy as np
from numpy import exp, sqrt, pi
from skued.simulation import structure_factor
import celltools.cell.tools as ct
from crystals.lattice import Lattice
from skued.simulation import affe

def diffraction_from_supercell(hkl: List[Tuple[int, int ,int]], scell: ct.super_cell, norm: bool = True):
    """

    Parameters
    ----------
    hkl: iterable of (3,)-tuples
        list of hkl indices to simulate the structure factor from
    cell: :class:`cell`
        cell to simulate the diffraction from
    norm: bool (default true)
        normalize to number of unit cells

    Returns
    -------
    """
    if norm:
        N = scell.size[0]*scell.size[1]*scell.size[2]
    else:
        N=1

    lattice = Lattice(scell.lattice.basis)
    _scatt_vector = []
    _amplitude = []
    for _hkl in hkl:
        _scatt_vector.append(lattice.scattering_vector(_hkl))

        _atms = []
        for atm in scell.base[0]:
            _atms.append(atm)
        for molc in scell.base[1]:
            for atm in molc.atoms:
                _atms.append((atm))

        _formfactor = []
        _coords = []
        for atm in _atms:
            _formfactor.append(affe(atm.element, np.sqrt(np.sum(np.square(_scatt_vector[-1])))))
            _coords.append(atm.coords.global_coord)

        ff = np.array(_formfactor)
        r = np.array(_coords)

        _amplitude.append(
            (1/N) * np.sum( ff * exp(-1j *  np.sum(_scatt_vector[-1]*r, axis=1) ))
        )

    return _scatt_vector, _amplitude
