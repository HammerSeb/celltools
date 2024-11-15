from copy import copy, deepcopy
from itertools import product
from typing import Union, List, Tuple, Literal, Callable

import numpy as np
from numpy import deg2rad

from celltools.linalg import BasisTransformation, Rotation
from celltools.linalg import Vector, Line, Plane, standard_basis
from celltools.linalg.find import average_plane
from . import Atom, Molecule, Cell


def move(obj: Union[Atom, Molecule], vec: Vector):
    """
    takes a cell object, e.g. atom or molecule and moves it in the direction of vector
    Parameters
    ----------
    obj: :class:`atom` or :class:`molecule`
    vec: :class:`vector`
    """
    if isinstance(obj, Atom):
        to_atm_basis = BasisTransformation(obj.coords.basis, vec.basis)
        obj.coords = obj.coords + to_atm_basis.inv_transform(vec)
    elif isinstance(obj, Molecule):
        to_mol_basis = BasisTransformation(obj.atoms[0].coords.basis, vec.basis)
        for atm in obj.atoms:
            atm.coords = atm.coords + to_mol_basis.inv_transform(vec)


def rotate(
    obj: Union[Molecule, List[Atom]],
    axis: Line,
    angle: float,
    mode: Union[Literal["deg"], Literal["rad"]] = "deg",
) -> None:
    """
    rotates object counterclockwise around line about specified angle
    Parameters
    ----------
    obj: list of :class:`atom`, :class:`molecule`
        objects to be rotated
    axis: :class:linalg.line
        axis around object is rotated
    angle: float
        rotation angle in degree or radians
    mode: "deg", "rad"
        specify angle input - "deg": for degree, "rad": for radial. default: "deg"

    """
    if isinstance(obj, list) and isinstance(obj[0], Atom):
        atom_list = obj
    elif isinstance(obj, Molecule):
        atom_list = obj.atoms
    else:
        raise ValueError("obj needs to be list of atoms or molecule")

    if mode == "deg":
        _rotation = Rotation(deg2rad(angle), axis.direction)
    elif mode == "rad":
        _rotation = Rotation(angle, axis.direction)
    else:
        raise ValueError("mode needs to be deg or rad")

    _to_line_basis = BasisTransformation(standard_basis, axis.basis)
    _offset = _to_line_basis.inv_transform(axis.origin)

    for atm in atom_list:
        _to_atm_basis = BasisTransformation(standard_basis, atm.coords.basis)
        atm.coords = _to_atm_basis.inv_transform(atm.coords)
        atm.coords -= _offset
        atm.coords = _rotation.rotate(atm.coords)
        atm.coords += _offset
        atm.coords = _to_atm_basis.transform(atm.coords)


def molecular_plane(molc: Molecule) -> Plane:
    """
    finds the average plane of a planar molecule.

    Parameters
    ----------
    molc: :class:`Molecule

    Returns
    -------
        :class:`Plane`
            average plane through all atoms of a given molecule

    """
    molc_coords = [atm.coords for atm in molc.atoms]
    return average_plane(molc_coords)


class SuperCell(Cell):
    """
    class defines a supercell generated from a unit cell with given size

    Parameters
    ----------
    unit_cell: :class:`Cell`
        unit cell from which supercell is generated
    size: tuple (int, int, int)
        specifies the size of the supercell along the three lattice dimensions of cell
    """

    def __init__(self, unit_cell: Cell, size: Tuple[int, int, int]):
        self._size = size
        self._atoms, self._molecules = [], []
        self.set_lattice(unit_cell.lattice)
        self._basevectors = [
            Vector([1, 0, 0], self.lattice),
            Vector([0, 1, 0], self.lattice),
            Vector([0, 0, 1], self.lattice),
        ]

        self._translation_vector = [
            l * Vector([1, 0, 0], self.lattice)
            + m * Vector([0, 1, 0], self.lattice)
            + n * Vector([0, 0, 1], self.lattice)
            for (l, m, n) in product(
                range(self.size[0]), range(self.size[1]), range(self.size[2])
            )
        ]

        for trans_vec in self._translation_vector:
            uc_atms, uc_molcs = unit_cell.base
            for _atm in uc_atms:
                _atm_new = copy(_atm)
                move(_atm_new, trans_vec)
                self.add_atom(_atm_new)
            for _molc in uc_molcs:
                _molc_new = deepcopy(_molc)
                move(_molc_new, trans_vec)
                self.add_molecule(_molc_new)

    @property
    def size(self):
        return self._size


class SymmetryOperator:
    """
    class specifing a symmetry operation in a crystal lattice to generate full basis from smallest possible basis

    Parameters
    ----------
    x_op, y_op, z_op: float
        parameters specifing symmetry operation on coordinates [x,y,z] in basis [b1, b2, b3] as
        [x_op(x,y,z), y_op(x,y,z), z_op(x,y,z)]
    label: str
        custom name of symmetry operation used in __repr__
        Should be chosen to be Seitz Symbol {R/m axis | translation}
    """

    def __init__(
        self,
        x_op: Callable,
        y_op: Callable,
        z_op: Callable,
        label: str = "Symmetry Operation",
    ) -> None:
        self._x_op = x_op
        self._y_op = y_op
        self._z_op = z_op

        self.label = label

    def __repr__(self):
        return f"< {self.label} >"

    @property
    def x_op(self):
        return self._x_op

    @property
    def y_op(self):
        return self._y_op

    @property
    def z_op(self):
        return self._z_op


def create_SymmetryOperator(generator_string: str) -> SymmetryOperator:
    """
    generates a symmetry operator from a generating string of teh form "operation on x, operation on y, operation on z"
    For example the identity is "x,y,z", a reflection on the xz plane is "x,-y,z"
    Parameters
    ----------
    generator_string: str
        generator string

    Returns SymmetryOperator
        generated symmetry operator
    -------

    """
    x_str, y_str, z_str = generator_string.split(",")

    def x_op(x, y, z):
        return eval(x_str)

    def y_op(x, y, z):
        return eval(y_str)

    def z_op(x, y, z):
        return eval(z_str)

    return SymmetryOperator(x_op, y_op, z_op, label=generator_string)


def generate_from_symmetry(atom: Atom, operator: SymmetryOperator) -> Tuple[Atom, bool]:
    """
    generates a new atom from a symmetry operation. The function returns an additional flag marking if the generated
    atom is truely a new position or if the atom is kept in place by the symmetry operator!
    Use together with an identity operation to generate full unit cell.
    Parameters
    ----------
    atom: :class:`atom`
        original atom
    operator: :class:`SymmetryOperator`
        symmetry operator generating a new set of atomic positions (x,y,z) -> (x',y',z')
    Returns
    -------
    :class:`atom`:
        new atom

    bool:
        return if the generated atom is kept in place (False) or if a new position is generated (True).
    """
    _atm = copy(atom)
    _coords = _atm.coords
    _atm.coords = Vector(
        [
            operator.x_op(_coords[0], _coords[1], _coords[2]),
            operator.y_op(_coords[0], _coords[1], _coords[2]),
            operator.z_op(_coords[0], _coords[1], _coords[2]),
        ],
        _coords.basis,
    )
    if _atm.label:
        _atm.label = _atm.label + "*"
    if _atm.coords == atom.coords:
        return _atm, False
    else:
        return _atm, True
