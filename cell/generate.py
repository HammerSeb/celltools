import re
import numpy as np
from numpy import deg2rad, cos, sin, sqrt
import cell.contents as cc
from linalg.basis import basis, vector

# TODO: atom list from cif, lattice from cif, unit cell from cif
# TODO: function to generate unit_cell from Crystals object


def lattice_from_cell_parameters(a, b, c, alpha, beta, gamma):
    """
    generates a :class:`lattice` from given cell parameters
    Parameters
    ----------
    a: float
        lattice distance a
    b: float
        lattice distance b
    c: float
        lattice distance c
    alpha: float
        lattice angle alpha in degree
    beta: float
        lattice angle beta in degree
    gamma: float
        lattice angle gamma in degree

    Returns
    -------
        :class:`lattice`
    """
    vec_a = np.array([a, 0, 0])
    vec_b = np.array([b * cos(deg2rad(gamma)),
                      b * sin(deg2rad(gamma)),
                      0 ])
    _c1 = c * cos(deg2rad(beta))
    _c2 = c * ( cos(deg2rad(alpha)) - cos(deg2rad(beta)) * cos(deg2rad(gamma)) ) / sin(deg2rad(gamma))
    _c3 = sqrt( c**2 - _c1**2 - _c2**2 )
    vec_c = np.array([_c1, _c2, _c3])
    return cc.lattice(basis(vec_a, vec_b, vec_c))

def cell_from_crystal(cryst):
    _basis = basis(*cryst.lattice_vectors)
    latt = cc.lattice(_basis)
    atms = []
    for atm in cryst.atoms:
        atms.append(cc.atom(atm.element, vector(atm.coords_fractional, latt)))
    return cc.cell(latt,atms)

def cell_from_cif(file, type="file"):
    """
    Generatures a :class:`cell` from a given cif (crystallographic information framework). Supported data type is cif1!
    Parameters
    ----------
    file: cif file
        input file with crystallographic information

    type: str
        "file" (default): atom list is generated from file without considering symmetries
        "sym": generate all atoms from listed symmetries - NOT IMPLEMENTED

    Returns: :class:`cell`
        crystal structure
    -------
    """
    with open(file, 'r') as cif:
        _latt_params = []
        for line in cif.readlines():
            #generates lattice from file
            line = line.lower()
            if "_cell_length_a" in line:
                _latt_params.append(re.findall("\d+\.\d+"), line)
            elif "_cell_length_b" in line:
                _latt_params.append(re.findall("\d+\.\d+"), line)
            elif "_cell_length_c" in line:
                _latt_params.append(re.findall("\d+\.\d+"), line)
            elif "_cell_length_alpha" in line:
                _latt_params.append(re.findall("\d+\.\d+"), line)
            elif "_cell_length_beta" in line:
                _latt_params.append(re.findall("\d+\.\d+"), line)
            elif "_cell_length_gamma" in line:
                _latt_params.append(re.findall("\d+\.\d+"), line)
        _latt = lattice_from_cell_parameters(*_latt_params)





def export_cell_to_cif(cell):
    pass



