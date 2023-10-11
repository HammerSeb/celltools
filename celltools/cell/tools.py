from typing import List, Tuple, Literal, Union
from itertools import  product
import numpy as np
from numpy import deg2rad
from copy import copy, deepcopy
from celltools.linalg.basis import standard_basis

from celltools.linalg import BasisTransformation, Rotation
from celltools.linalg import Vector, Line, standard_basis
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


def rotate(obj: Union[Molecule, List[Atom]], axis: Line, angle: float,
           mode: Union[Literal["deg"], Literal["rad"]] = "deg") -> None:
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
        self._basevectors = [Vector([1, 0, 0], self.lattice),
                             Vector([0, 1, 0], self.lattice),
                             Vector([0, 0, 1], self.lattice)]

        self._translation_vector = [l * Vector([1, 0, 0], self.lattice) +
                                    m * Vector([0, 1, 0], self.lattice) +
                                    n * Vector([0, 0, 1], self.lattice)
                                    for (l, m, n) in
                                    product(range(self.size[0]), range(self.size[1]), range(self.size[2]))]

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
    a, b, c, x0, y0, z0: float
        parameters specifing symmetry operation on coordinates [x,y,z] as
        [a*x + x0, b*y + y0, c*z + z0]
    label: str
        custom name of symmetry operation used in __repr__
    id: list
        list of coordinates that are kept in place by symmetry operation, e.g. [[0,0,0]] for inversion [-x, -y, -z].
        leave empty for identity
    """

    def __init__(self, a, b, c, x0, y0, z0, id=[], label=None):
        self._a = a
        self._b = b
        self._c = c
        self._x0 = x0
        self._y0 = y0
        self._z0 = z0
        self._id = id

        self.label = label

    def __repr__(self):
        if self.label:
            return f"< Symm {self.label} >"
        else:
            return f"< Symm {self.a:.0f}x+{self.x0:.2f} | {self.b:.0f}y+{self.y0:.2f} | {self.c:.0f}y+{self.z0:.2f} >"

    def __eq__(self, other):
        if isinstance(other, SymmetryOperator):
            return np.all([
                self.a == other.a, self.b == other.b, self.c == other.c,
                self.x0 == other.x0, self.y0 == other.y0, self.z0 == other.z0
            ])
        else:
            raise TypeError(f"cannot compare symmetry_operator with {type(other)}")

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @property
    def c(self):
        return self._c

    @property
    def x0(self):
        return self._x0

    @property
    def y0(self):
        return self._y0

    @property
    def z0(self):
        return self._z0

    @property
    def id(self):
        return self._id

    def addid(self, coord):
        """
        Parameters
        ----------
        coord: list (3,0
            list of coordinates not affected by symmetry operation
        """
        self._id.append(coord)


def generate_from_symmetry(atom: Atom, operator: SymmetryOperator) -> Tuple[Atom, bool]:
    """
    generates a new atom from a symmetry operation
    Parameters
    ----------
    atom: :class:`atom`
        original atom
    operator: :class:`symmetry_operator

    Returns
    -------
    :classl:`atom`:
        new atom
    bool:
        returns if performed symmetry operation is valid according to operator.id values. Only use atoms for which the
        returned value is True.


    """
    _atm = copy(atom)
    _coords = _atm.coords
    if not operator.id:
        # this escapes the identity
        return _atm, True
    else:
        for id in operator.id:
            if np.all(_atm.coords.vector == id):
                return _atm, False
        _atm.coords = Vector(
            [operator.a * _coords[0] + operator.x0,
             operator.b * _coords[1] + operator.y0,
             operator.c * _coords[2] + operator.z0,
             ], _coords.basis
        )
        if _atm.label:
            _atm.label = _atm.label + "*"
        return _atm, True
