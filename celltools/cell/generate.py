import re
from copy import copy
from os import PathLike
from typing import Union, Literal, Optional

import numpy as np
from CifFile import ReadCif
from crystals import Crystal, Lattice, AtomicStructure
from numpy import deg2rad, cos, sin, sqrt

from celltools.cell.spacegroup_data import SPACE_GROUP, read_sym_file
from celltools.cell.tools import generate_from_symmetry, create_SymmetryOperator
from celltools.linalg import Basis, Vector, BasisTransformation, standard_basis
from . import Atom, Molecule, Lattice, Cell


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
        _latt_params.append(float(re.findall('\d*\.?\d*',cif[block])[0]))
    _latt = lattice_from_cell_parameters(*_latt_params)

    _elem = list(map(lambda label: re.findall("\D+", label), cif["_atom_site_label"]))
    _coords = list(map(
        lambda coord: Vector([float(re.findall('\-?\d*\.?\d*', coord[0])[0]),
                              float(re.findall('\-?\d*\.?\d*',coord[1])[0]),
                              float(re.findall('\-?\d*\.?\d*',coord[2])[0])],
                              _latt),
        zip(cif["_atom_site_fract_x"], cif["_atom_site_fract_y"], cif["_atom_site_fract_z"])
    ))
    _atms = []
    for el, coord, label in zip(_elem, _coords, cif["_atom_site_label"]):
        _atms.append(Atom(el[0], coord, label))

    if "_symmetry_int_tables_number" in cif:
        sg_number = cif["_symmetry_int_tables_number"]
    elif "_space_group_it_number" in cif:
        sg_number = cif["_space_group_it_number"]
    else:
        sg_number = '0'

    if mode == 'sym':
        _atmssym = []
        symmetries = []
        if "_space_group_symop_operation_xyz" in cif:
            for gen_string in cif["_space_group_symop_operation_xyz"]:
                symmetries.append(create_SymmetryOperator(gen_string))
        elif "_symmetry_equiv_pos_as_xyz" in cif:
            # only for older cifs
            for gen_string in cif["_symmetry_equiv_pos_as_xyz"]:
                symmetries.append(create_SymmetryOperator(gen_string))
        else:
            symmetries = read_sym_file(SPACE_GROUP[sg_number])

        for atm in _atms:
            # generate atms from identity operator (leave positions unchanged)
            _atmssym.append(copy(atm))
        for operator in symmetries:
            # generate atms from symmetry element except identity
            for atm in _atms:
                _atm, sym_out = generate_from_symmetry(atm, operator)
                if sym_out:
                    _atmssym.append(_atm)
        return Cell(_latt, _atmssym)
    else:
        return Cell(_latt, _atms)


def qe_export_molecule(molc: Molecule, file: Optional[Union[str, PathLike]] = None) -> None:
    """
    Exports atom coordinates to ATOMIC_SPECIES and ATOMIC_POSITIONS card  as well as nat and ntyp values.
    Units of Angstrom are assumed for export of atomic positions in the standard basis. Use ibrav=0 or ibrav=1.
    If no file is specified, card is printed into stdout.
    Parameters
    ----------
    molc: class:`Molecule`
    file: filepath (Optional)

    """
    lines = []
    ### SYSTEM and ATOMIC_SPECIES CARD
    _elements = []
    _atomtypes = []
    for atm in molc.atoms:
        if atm.element not in _elements:
            _elements.append(atm.element)
            _atomtypes.append(atm)

    lines.append("&SYSTEM\n")
    lines.append(f"nat = {len(molc.atoms)}\n")
    lines.append(f"ntyp = {len(_elements)}\n")
    lines.append('\n')

    lines.append("ATOMIC_SPECIES\n")
    for atm in _atomtypes:
        lines.append(f"{atm.element} {atm.mass} {atm.element}.pseudo\n")
    lines.append("\n")

    ### ATOMIC_POSITIONS CARD
    _molcbasis = molc.atoms[0].coords.basis
    to_std_basis = BasisTransformation(_molcbasis, standard_basis)

    lines.append("ATOMIC_POSITIONS 'angstrom'\n")
    for atm in molc.atoms:
        _coords = to_std_basis.transform(atm.coords)
        lines.append(f"{atm.label} {_coords[0]:.5f} {_coords[1]:.5f} {_coords[2]:.5f}\n")

    if not file:
        for line in lines:
            print(line.strip("\n"))
    else:
        with open(file, "w") as f:
            f.writelines(lines)


def qe_export_cell(cell: Cell, file: Optional[Union[str, PathLike]] = None) -> None:
    """
    Exports given cell to CELL_PARAMETERS, ATOMIC_SPECIES and ATOMIC_POSITIONS card as well as nat and ntyp values.
    Needs to be used with ibrav=0 (no explotation of symmetry! Handle with care!). Atom positions are given in units of
    lattice (alat), which is assumed to be in angstrom. If no file is specified, card is printed into stdout.
    Parameters
    ----------
    cell: :class:`Cell`
    file: filepath (Optional)
    """
    lines = []
    ### CELL_PARAMETERS CARD
    lines.append("CELL_PARAMETERS 'angstrom'\n")
    for vec in cell.lattice:
        lines.append(f"{vec[0]:.5f} {vec[1]:.5f} {vec[2]:.5f}\n")
    lines.append("\n")

    ### SYSTEM and ATOMIC_SPECIES CARD
    _elements = []
    _atomtypes = []
    nat = 0
    for atm in cell.atoms:
        nat +=1
        if atm.element not in _elements:
            _elements.append(atm.element)
            _atomtypes.append(atm)
    for molc in cell.molecules:
        nat += 1
        for atm in molc.atoms:
            if atm.element not in _elements:
                _elements.append(atm.element)
                _atomtypes.append(atm)

    lines.append("&SYSTEM\n")
    lines.append(f"nat = {nat}\n")
    lines.append(f"ntyp = {len(_elements)}\n")
    lines.append('\n')

    lines.append("ATOMIC_SPECIES\n")
    for atm in _atomtypes:
        lines.append(f"{atm.element} {atm.mass} {atm.element}.pseudo\n")
    lines.append("\n")

    ### ATOMIC_POSITIONS CARD
    lines.append("ATOMIC_POSITIONS 'crystal'\n")
    for atm in cell.atoms:
        lines.append(f"{atm.label} {atm.coords[0]:.5f} {atm.coords[1]:.5f} {atm.coords[2]:.5f}\n")
    for molc in cell.molecules:
        for atm in molc.atoms:
            lines.append(f"{atm.label} {atm.coords[0]:.5f} {atm.coords[1]:.5f} {atm.coords[2]:.5f}\n")

    if not file:
        for line in lines:
            print(line.strip("\n"))
    else:
        with open(file, "w") as f:
            f.writelines(lines)


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
