import typing
from typing import List, Tuple

import numpy as np
from numpy import exp, sqrt, pi

from . import IndexLike
from celltools import SuperCell
from crystals.lattice import Lattice
from skued.simulation import affe

def diffraction_from_supercell(hkl: List[IndexLike], supercell: SuperCell, norm: bool = True):
    """
    calculates the scattering vectors and respective complex scattering amplitude for a list of hkl indices from a
    super cell. Uses atomic form factors form scikit-ued and crystals.lattice.Lattice to calculate the scattering
    vectors from the base unit cell.
    Parameters
    ----------
    hkl: iterable of (3,)-tuples of length N
        list of hkl indices to simulate the structure factor from
    supercell: :class:`SuperCell`
        super cell to simulate the diffraction from

    Returns
    -------
        scatt_vector: list of nd.arrays size (3,) of length N
            scattering vectors q of respective hkl index in reciprocal space
        amplitude: list of complex numbers of length N
            amplitude of reflection of respective hkl index
    """
    if norm:
        N = supercell.size[0]*supercell.size[1]*supercell.size[2]
    else:
        N=1

    lattice = Lattice(supercell.lattice.basis)
    _scatt_vector = []
    _amplitude = []
    for _hkl in hkl:
        _scatt_vector.append(lattice.scattering_vector(_hkl))

        _atms = []
        for atm in supercell.base[0]:
            _atms.append(atm)
        for molc in supercell.base[1]:
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
