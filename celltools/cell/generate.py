import re
from os import PathLike
from typing import Union, Literal

import numpy as np
from CifFile import ReadCif
from crystals import Crystal, Lattice, AtomicStructure, Atom
from numpy import deg2rad, cos, sin, sqrt

from celltools.cell.spacegroup_data import SPACE_GROUP
from celltools.cell.tools import generate_from_symmetry
from celltools.linalg import Basis, Vector
from . import Atom, Lattice, Cell


def lattice_from_cell_parameters(a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> Lattice:
    """
    generates a :class:`Lattice` from given cell parameters
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
        :class:`Lattice`
    """
    vec_a = np.array([a, 0, 0])
    vec_b = np.array([b * cos(deg2rad(gamma)),
                      b * sin(deg2rad(gamma)),
                      0])
    _c1 = c * cos(deg2rad(beta))
    _c2 = c * (cos(deg2rad(alpha)) - cos(deg2rad(beta)) * cos(deg2rad(gamma))) / sin(deg2rad(gamma))
    _c3 = sqrt(c ** 2 - _c1 ** 2 - _c2 ** 2)
    vec_c = np.array([_c1, _c2, _c3])
    return Lattice(Basis(vec_a, vec_b, vec_c))


def cell_from_crystal(cryst: Crystal) -> Cell:
    """
    creats a :class:`cell` object from a :class:`crystals.Crystal` object
    Parameters
    ----------
    cryst

    Returns
    -------

    """
    _basis = Basis(*cryst.lattice_vectors)
    latt = Lattice(_basis)
    atms = []
    for atm in cryst.atoms:
        atms.append(Atom(atm.element, Vector(atm.coords_fractional, latt)))
    return Cell(latt, atms)


def cell_to_crystal(cell: Cell) -> Crystal:
    """
    crystals and scikit-ued interface
    transforms a cell instance to a Crystal instance
    Parameters
    ----------
    cell: :class:`Cell`

    Returns
    -------
        :class:`Crystal`
    """
    _lattice = Lattice(cell.lattice.basis)
    _atoms = []
    for atm in cell.atoms:
        _atoms.append(Atom(atm.element, atm.coords.vector, _lattice))
    for molc in cell.molecules:
        for atm in molc:
            _atoms.append(Atom(atm.element, atm.coords.vector, _lattice))
    return Crystal(AtomicStructure(_atoms), _lattice.lattice_vectors)


def cell_from_cif(file: Union[str, PathLike], mode: Union[Literal["file"], Literal["sym"]] = "sym") -> Cell:
    """
    Generatures a :class:`Cell` from a given cif (crystallographic information framework). Supported data type is cif1!
    Parameters
    ----------
    file: cif file
        input file with crystallographic information

    mode: str
        "file": atom list is generated from file without considering symmetries
        "sym" (default): generate all atoms from listed symmetries
            - NOT ALL SYMMETRIES ARE IMPLEMENTED, CHECK cell.spacegroup_data TO CHECK -

    Returns
    -------
        :class:`Cell`
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

    _elem = list(map(lambda label: re.findall("\D+", label), cif["_atom_site_label"]))
    _coords = list(map(
        lambda coord: Vector([float(coord[0]), float(coord[1]), float(coord[2])], _latt),
        zip(cif["_atom_site_fract_x"], cif["_atom_site_fract_y"], cif["_atom_site_fract_z"])
    ))
    _atms = []
    for el, coord, label in zip(_elem, _coords, cif["_atom_site_label"]):
        _atms.append(Atom(el[0], coord, label))

    if cif["_symmetry_int_tables_number"] in SPACE_GROUP.keys() and mode == "sym":
        _atmssym = []
        for operator in SPACE_GROUP[cif["_symmetry_int_tables_number"]]:
            # generate atms from symmetry element except inversion
            for atm in _atms:
                _atm, sym_out = generate_from_symmetry(atm, operator)
                if sym_out:
                    _atmssym.append(_atm)
        return Cell(_latt, _atmssym)
    else:
        return Cell(_latt, _atms)


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
    # TODO


def _export_lattice_to_cif(lattice, file=None):
    # TODO
    pass


def export_cell_to_cif(cell):
    # TODO
    pass
