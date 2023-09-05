from copy import copy
import cell.contents as cc
from linalg.basis import basis, vector
from linalg.transformations import basis_transformation
#TODO: function to generate avergage plane through a list of atoms

def move(obj, vec):
    """
    ONLY FUNCTION SKETCH
    takes a cell object, e.g. atom or molecule and moves it in the direction of vector
    Parameters
    ----------
    obj: :class:`atom` or :class:`molecule`
    vec: :class:`vector`
    """
    if isinstance(obj, cc.atom):
        to_atm_basis = basis_transformation(obj.coords.basis, vec.basis)
        obj.coords = obj.coords + to_atm_basis.inv_transform(vec)
    elif isinstance(obj, cc.molecule):
        to_mol_basis = basis_transformation(obj.atoms[0].coords.basis, vec.basis)
        for atm in obj.atoms:
            atm.coords = atm.coords + to_mol_basis.inv_transform(vec)

def generate_from_symmetry(atom, operator):
    """
    generates a new atom from a symmetry operation
    Parameters
    ----------
    atom: :class:`atom`
        original atom
    operator: :class:`symmetry_operator

    Returns
    -------
    :classl:`atom`
        new atom

    """
    _atm = copy(atom)
    _coords =  _atm.coords
    _atm.coords = vector(
        [operator.a*_coords[0]+operator.x0,
         operator.b*_coords[1]+operator.y0,
         operator.c*_coords[2]+operator.z0,
         ], _coords.basis
    )
    if _atm.label:
        _atm.coords = _atm.label + "*"
    return _atm



class symmetry_operator:
    def __init__(self,a ,b ,c ,x0 , y0, z0):
        """

        Parameters
        ----------
        a
        b
        c
        x0
        y0
        z0
        """
        self._a = a
        self._b = b
        self._c = c
        self._x0 = x0
        self._y0 = y0
        self._z0 = z0

    def __repr__(self):
        return f"Symm {self.a:.0f}x+{self.x0:.2f} | {self.b:.0f}y+{self.y0:.2f} | {self.c:.0f}y+{self.z0:.2f} "

    def __eq__(self, other):
        if isinstance(other, symmetry_operator):
            return np.all([
                self.a == other.a, self.b == other.b, self.c == other.c,
                self. x0 == other.x0, self.y0 == other.y0, self.z0 == other.z0
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
