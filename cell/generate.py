import re
from itertools import combinations
import numpy as np
from numpy import deg2rad, cos, sin, sqrt
import cell.contents as cc
import cell.tools as ct
from cell.spacegroup_data import SPACE_GROUP
from linalg.basis import basis, vector

from CifFile import ReadCif

def auto_bonds(atm_list, rmin=1, rmax=2):
    """
    generates list of bonds from given atom list. A bond between two atoms is generated if the distances of the two
    atoms lies between rmin and rmax
    Parameters
    ----------
    atm_list: list of :class:`atom`
        atom list from which bonds are created
    rmin: float
        minimum distance between two atoms to create a bond
    rmax: float
        maximum distance between two atoms to create a bond

    Returns
    -------
        list of :class:`bond`
    """
    bonds = []
    atm_pairs = combinations(atm_list, 2)
    for pair in atm_pairs:
        _dist = (pair[0].coords - pair[1].coords).abs_global
        if rmin < _dist <rmax:
            bonds.append(cc.bond(pair[0], pair[1]))
    return bonds

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
    cif = ReadCif(file)[ReadCif(file).keys()[0]]

    ltt_blocks = ["_cell_length_a", "_cell_length_b", "_cell_length_c",
                   "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"]

    # creating lattice from cif
    _latt_params = []
    for block in ltt_blocks:
        _latt_params.append(float(cif[block]))
    _latt = lattice_from_cell_parameters(*_latt_params)

    _elem = list(map(lambda label : re.findall("\D+" ,label) , cif["_atom_site_label"]))
    _coords = list(map(
        lambda coord: vector([float(coord[0]), float(coord[1]), float(coord[2])], _latt),
        zip(cif["_atom_site_fract_x"], cif["_atom_site_fract_y"], cif["_atom_site_fract_z"])
    ))
    _atms = []
    for el, coord, label in zip(_elem, _coords, cif["_atom_site_label"]):
        _atms.append(cc.atom(el[0],coord,label))

    if cif["_symmetry_int_tables_number"] in SPACE_GROUP.keys():
        _atmssym = []
        for operator in SPACE_GROUP[cif["_symmetry_int_tables_number"]]:
            # generate atms from symmetry element except inversion
            for atm in _atms:
                _atm, sym_out = ct.generate_from_symmetry(atm, operator)
                if sym_out:
                    _atmssym.append(_atm)
        return cc.cell(_latt, _atmssym)
    else:
        return cc.cell(_latt, _atms)

def _export_atom_list_to_cif(atoms, file=None):
    """
    atom list is formatted in cif format and written to stdout or file
    Parameters
    ----------
    atoms: iterable of :class:`cell.contents.atoms`
        list of atoms for cif file.
    file: file path (default: None)
        file path to save output too, if left unspecified the formatted atom list is returned as a list of strings where
        each entry is one line of the cif file.

    Returns
    -------
    list of str
        lines of output file as list

    """

def _export_lattice_to_cif(lattice, file=None):
    #TODO
    pass

def export_cell_to_cif(cell):
    #TODO
    pass




