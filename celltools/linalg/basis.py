from copy import copy
from typing import Union, List, Tuple, Callable, Literal

import numpy as np
from numpy.linalg import det


class LinearAlgebraError(Exception):
    """ custom exception for linear algebra errors"""
    pass


def check_linear_independence() -> Callable:
    """
    Decorator function to check if the three vectors passed to the target function are linear independent
    (assumes to be wrapper for basis class init function)
    Raises
    ------
        LinearAlgebraError
            if the vectors are not linear independent
    """

    def check_accepts(function: Callable) -> Callable:
        def new_function(*args, **kwargs) -> Callable:
            basis = np.array([*args[-3:]])
            if det(basis) == 0:
                raise LinearAlgebraError("vectors are not linear independent")

            return function(*args[-4:], **kwargs)

        new_function.__name__ = function.__name__

        return new_function

    return check_accepts


def parallel(vec1, vec2):
    """checks if two vectors are parallel in standard coordinates"""
    if abs(np.dot(vec1.global_coord, vec2.global_coord)) == vec1.abs_global * vec2.abs_global:
        return True
    else:
        return False


VectorLike = Union[Tuple[float, float, float], np.ndarray]
IndexLike = Union[Literal[1], Literal[2], Literal[3]]


class Basis:
    """
    Class for a linear algebra basis.
    Basis object is generatred from three linear independent vectors v1, v2, v3.
    An offset can be specified to shift origin of basis system within a global reference frame

    Parameters
    ----------
    v1: list, tuple or nd.array
        three entries as coordinates v11, v12, v13
    v2: list, tuple or nd.array
        three entries as coordinates v21, v22, v23
    v3: list, tuple or nd.array
        three entries as coordinates v31, v32, v33
    offset: list, tuple or nd.array
        three entries as coordinates in global reference frame. default is [0,0,0]

    Raises
    ------
        LinearAlgebraError
            if input vectors are not linear independent
    """

    @check_linear_independence()
    def __init__(self, v1: VectorLike, v2: VectorLike, v3: VectorLike, offset: VectorLike = np.zeros((3,))):
        self._basis = np.array([v1, v2, v3])

        self._offset = offset

    def __repr__(self) -> str:
        return f'<basis {self.basis[0]}, {self.basis[1]}, {self.basis[2]}>'

    def __eq__(self, other: 'Basis') -> bool:
        """ checks only if all basis vectors are the same, does not check for offset! """
        if isinstance(other, Basis):
            if np.all(self.v1 == other.v1) and np.all(self.v2 == other.v2) and np.all(self.v3 == other.v3):
                return True
            else:
                return False
        else:
            raise TypeError("must be compared to instance of basis")

    def __getitem__(self, index: int) -> np.ndarray:
        if not isinstance(index, int):
            raise ValueError("index needs to be integer")
        elif abs(index) > 2 and index != -3:
            raise IndexError("index out of bounds")
        else:
            if index == 0 or index == -3:
                return self.basis[0]
            elif index == 1 or index == -2:
                return self.basis[1]
            elif index == 2 or index == -1:
                return self.basis[2]

    @property
    def basis(self) -> np.ndarray:
        """ returns basis systems """
        return self._basis

    @property
    def offset(self) -> np.ndarray:
        """ returns offset """
        return self._offset

    @property
    def v1(self) -> np.ndarray:
        """return v1 in global coordinates"""
        return self.basis[0] + self.offset

    @property
    def v2(self) -> np.ndarray:
        """return v2 in global coordinates"""
        return self.basis[1] + self.offset

    @property
    def v3(self) -> np.ndarray:
        """return v3 in global coordinates"""
        return self.basis[2] + self.offset

    @offset.setter
    def offset(self, offset: VectorLike):
        """
        sets offset of basis system
        Parameters
        ----------
        offset:  list, tuple or nd.array
            three entries as coordinates in global reference frame
        """
        if isinstance(offset, np.ndarray):
            self._offset = offset
        else:
            self._offset = np.array(offset)

    def permute(self, order: Tuple[IndexLike, IndexLike, IndexLike]):
        """
        permutes the basis vectors specified by a 3-tuple of (new position v1, new position v2, new position v3) given
        by numbers 1, 2 and 3.
        Parameters
        ----------
        order: tuple of size 3
        """
        if np.sum(order) != 6:
            raise ValueError("order must be given in integers of 1,2 and 3 only")
        else:
            self._basis = np.array([self._basis[order[0] - 1], self._basis[order[1] - 1], self._basis[order[2] - 1]])


standard_basis = Basis([1, 0, 0], [0, 1, 0], [0, 0, 1])


class Vector:
    """
    Class representing vectors.
    Vector is given by three coordinate entries and a basis in which the coordinates are
    defined. Vectors can be added/substracted if they are in the basis. Scalar multiplication is possible as well,
    division can be achieved by *1/divisor.
    Parameters
    ----------
    vector: list, tuple or nd.array
        three entries as coordinates v1, v2, v3
    basis: :class:`Basis`
        basis in which coordinates are specified. default is standard_basis
    """
    def __init__(self, vector: VectorLike, basis: Basis = standard_basis):

        self.vector = np.array([vector[0], vector[1], vector[2]])
        self._basis = basis

    def __repr__(self) -> str:
        return f"{self.vector}"

    def __add__(self, other: 'Vector') -> 'Vector':
        """ addition of 2 vector instances """
        if not self.basis == other.basis:
            raise LinearAlgebraError("basis do not match")
        else:
            return Vector([self.vector[0] + other.vector[0],
                           self.vector[1] + other.vector[1],
                           self.vector[2] + other.vector[2]],
                          self.basis)

    def __sub__(self, other: 'Vector') -> 'Vector':
        """ subtraction of 2 vector instances """
        if not self.basis == other.basis:
            raise LinearAlgebraError("basis do not match")
        else:
            return Vector([self.vector[0] - other.vector[0],
                           self.vector[1] - other.vector[1],
                           self.vector[2] - other.vector[2]],
                          self.basis)

    def __mul__(self, other: float) -> 'Vector':
        """ multiplication with a scalar """
        if isinstance(other, float) or isinstance(other, int):
            return Vector([self.vector[0] * other,
                           self.vector[1] * other,
                           self.vector[2] * other],
                          self.basis)
        else:
            raise ValueError("only scalar multiplication with real numbers allowed")

    def __rmul__(self, other: float) -> 'Vector':
        """ multiplication with a scalar """
        if isinstance(other, float) or isinstance(other, int):
            return Vector([self.vector[0] * other,
                           self.vector[1] * other,
                           self.vector[2] * other],
                          self.basis)
        else:
            raise ValueError("only scalar multiplication with real numbers allowed")

    def __eq__(self, other: 'Vector') -> bool:
        """ checks if vectors are the same in global reference frame """
        if isinstance(other, Vector):
            return np.all(self.global_coord == other.global_coord)
        else:
            raise TypeError("must be compared to vector instance")

    def __getitem__(self, index: int) -> float:
        if not isinstance(index, int):
            raise ValueError("index needs to be integer")
        elif abs(index) > 2 and index != -3:
            raise IndexError("index out of bounds")
        else:
            if index == 0 or index == -3:
                return self.vector[0]
            elif index == 1 or index == -2:
                return self.vector[1]
            elif index == 2 or index == -1:
                return self.vector[2]

    @property
    def abs(self) -> float:
        """
        returns absolute value of vector length in its basis coordinates. This might only be a usefull metric in
        some basis systems. Use with care!
        """
        return np.sqrt(self.vector[0] ** 2 + self.vector[1] ** 2 + self.vector[2] ** 2)

    @property
    def global_coord(self) -> np.ndarray:
        """return vector in global coordinates"""
        return self.vector[0] * self.basis[0] + self.vector[1] * self.basis[1] + self.vector[2] * self.basis[2]

    @property
    def abs_global(self) -> float:
        """ returns absolute length of vector in global reference frame """
        return np.sqrt(self.global_coord[0] ** 2 + self.global_coord[1] ** 2 + self.global_coord[2] ** 2)

    @property
    def basis(self):
        return self._basis

    def parallel(self, other: 'Vector') -> bool:
        """
        checks if vector is parallel to another vector

        Parameters
        ----------
        other: :class:`Vector`

        Returns
        -------
            bool
        """
        if parallel(self, other):
            return True
        else:
            return False


class Line:
    """
    Class representing a line through a point in a specific direction. All vectors must be given in same basis
    Parameters
    ----------
    origin: :class:`Vector`
        defines point on line
    direction: :class:`Vector`
        defines direction of line
    """
    def __init__(self, origin: Vector, direction: Vector):
        if origin.basis == direction.basis:
            self._origin = origin
            self._direction = direction
            self._basis = self.origin.basis
        else:
            raise LinearAlgebraError("input vectors need to be expressed in the same basis")

    def __repr__(self) -> str:
        return f"< line: {self.origin} + k*{self.direction} >"

    def __eq__(self, other: 'Line') -> bool:
        if isinstance(other, Line):
            if self.on_line(other.origin) and self.parallel(other):
                return True
            else:
                return False
        else:
            raise ValueError("must be compared to line instance")

    @property
    def origin(self) -> Vector:
        """ return origin """
        return self._origin

    @property
    def direction(self) -> Vector:
        """ return direction """
        return self._direction

    @property
    def basis(self) -> Basis:
        """ return basis in which line is defined"""
        return self._basis

    def on_line(self, point: Vector) -> bool:
        """
        Check if a given point lies on the line.
        Parameters
        ----------
        point: :class:`Vector`

        Returns
        -------
            bool
        """
        if ((point.global_coord[0] - self.origin.global_coord[0]) / self.direction.global_coord[0]) == \
                ((point.global_coord[1] - self.origin.global_coord[1]) / self.direction.global_coord[1]) and \
                ((point.global_coord[1] - self.origin.global_coord[1]) / self.direction.global_coord[1]) == \
                ((point.global_coord[2] - self.origin.global_coord[2]) / self.direction.global_coord[2]) and \
                ((point.global_coord[2] - self.origin.global_coord[2]) / self.direction.global_coord[2]) == \
                ((point.global_coord[0] - self.origin.global_coord[0]) / self.direction.global_coord[0]):
            return True
        else:
            return False

    def parallel(self, line: 'Line') -> bool:
        """
        checks if line is parralel to another line
        Parameters
        ----------
        line: :class:`line`
            other line
        Returns
        -------
            bool
        """
        if parallel(self.direction, line.direction):
            return True
        else:
            return False

    def distance(self, point: Vector) -> float:
        """
        distance from line to point in standard basis
        formula adapted from https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        Parameters
        ----------
        point: :class:`Vector`

        Returns
        -------
            float
        """
        return ((Vector(self.origin.global_coord - point.global_coord).abs ** 2 * self.direction.abs_global ** 2
                 - np.dot((self.origin.global_coord - point.global_coord), self.direction.global_coord) ** 2)
                / self.direction.abs_global ** 2)

    def point_on_line(self, k: float) -> Vector:
        """
        returns a point which is given by self.origin + k * self.direction
        Parameters
        ----------
        k: float

        Returns
        -------
            :class:`vector`
        """
        return self.origin + k * self.direction


class Plane:
    """
    Class representing a plane defined by a point on the plane and a normal vector.
    The plane is defined in the global reference frame but can be defined by any two vectors and the same basis

    Parameters
    ----------
    origin: :class:`Vector`
    normal: :class:`Vector`
    """
    def __init__(self, origin: Vector, normal: Vector):

        if origin.basis == normal.basis:
            self._origin = origin
            self._normal = normal
            self._define_basis()
        else:
            raise LinearAlgebraError("input vectors need to be expressed in the same basis")

    def __repr__(self) -> str:
        return (f"< plane: {self.parametric_form[0]:.2f}*x + {self.parametric_form[1]:.2f}*y "
                f"+ {self.parametric_form[2]:.2f}*z = {self.parametric_form[3]:.2f} >")

    @property
    def normal(self) -> Vector:
        """ returns normal vector """
        return self._normal

    @property
    def origin(self) -> Vector:
        """ returns origin vecttor """
        return self._origin

    @property
    def basis(self) -> Basis:
        """
        returns a basis of the plane. Usually that means that two basis vectors lie within the plane, and the third one
        is the normal vector defining the plane.
        (This property has a setter, so you can change this basis to whatever basis you like, which means it might not
        be a basis related to the plane any more.)
        Returns
        -------
            Basis
                basis related to the plane
        """
        return self._basis

    @basis.setter
    def basis(self, basis: Basis):
        """
        Setter for the basis property. Be careful when using
        Parameters
        ----------
        basis: :class:`basis`
        """
        if isinstance(basis, Basis):
            self._basis = basis
        else:
            raise ValueError("input needs to be of type basis")

    @property
    def parametric_form(self) -> np.ndarray:
        """
        gives an array containing the for parameters of the normalized parametric form n1*x + n2*y +n3*z = d
        in the standard basis. ni are the coordinates of the normalized normal vector on the plane

        Returns
        -------
            nd.array
                [n1, n2, n3, d]
        """
        return np.array([self.normal.global_coord[0] / self.normal.abs_global,
                         self.normal.global_coord[1] / self.normal.abs_global,
                         self.normal.global_coord[2] / self.normal.abs_global,
                         np.dot(self.normal.global_coord / self.normal.abs_global, self.origin.global_coord)])

    def on_plane(self, point: Vector) -> bool:
        """
        checks if given point lies withing the plane
        Parameters
        ----------
        point: :class:`Vector`

        Returns
        -------
            bool
        """
        if np.dot(self.normal.global_coord, point.global_coord) == self.parametric_form[3]:
            return True
        else:
            return False

    def distance(self, point: Vector):
        """
        returns the distance from the plane to a given point
        formula adapted from: https://www.cuemath.com/geometry/distance-between-point-and-plane/

        Parameters
        ----------
        point: :class:`Vector`

        Returns
        -------
            float
        """
        return (abs((np.sum(self.parametric_form[:3] * point.global_coord) - self.parametric_form[3]))
                / np.sqrt(np.sum(np.square(self.parametric_form[:3]))))

    def _define_basis(self):
        """
        defines orthonormal basis of the plane centered at origin. Two basis vectors lie within the plane, third vector
        is the in direction of the normal vector.
        """
        # parallel to xy, xz, yz plane
        if np.count_nonzero(self.parametric_form[:-1] == 0) == 2:
            self._basis = copy(standard_basis)
            if not self.parametric_form[0] == 0:
                self._basis.permute((2, 3, 1))
            if not self.parametric_form[1] == 0:
                self._basis.permute((1, 3, 2))

            self._basis.offset = self.origin.global_coord
        # parallel to x,y or z axis
        elif np.count_nonzero(self.parametric_form[:-1] == 0) == 1:
            if self.parametric_form[0] == 0:  # x-axis
                self._basis = Basis(
                    [0, np.sqrt(2), np.sqrt(2)], [0, -1 * np.sqrt(2), np.sqrt(2)], [1, 0, 0]
                )

            elif self.parametric_form[1] == 0:  # y-axis
                self._basis = Basis(
                    [np.sqrt(2), 0, np.sqrt(2)], [np.sqrt(2), 0, -1 * np.sqrt(2)], [0, 1, 0]
                )

            elif self.parametric_form[2] == 0:  # z-axis
                self._basis = Basis(
                    [np.sqrt(2), np.sqrt(2), 0], [-1 * np.sqrt(2), np.sqrt(2), 0], [0, 0, 1],
                )

            self._basis.offset = self.origin.global_coord
            # plane with no restrictions
        else:
            _origin = Vector(self.origin.global_coord)
            _x_intsec = Vector([self.parametric_form[3] / self.parametric_form[0], 0, 0])
            _y_intsec = Vector([0, self.parametric_form[3] / self.parametric_form[1], 0])
            _basis3 = self.normal.global_coord / self.normal.abs_global
            if not _origin == _x_intsec:
                _basis1 = ((_x_intsec - _origin) * (1 / Vector(_x_intsec - _origin).abs_global)).global_coord
            else:
                _basis1 = ((_y_intsec - _origin) * (1 / Vector(_y_intsec - _origin).abs_global)).global_coord
            _basis2 = Vector(np.cross(_basis1, _basis3) * (
                        1 / Vector(np.cross(_basis1, _basis3)).abs_global)).global_coord
            self._basis = Basis(_basis1, _basis2, _basis3)
            self.basis.offset = self.origin.global_coord
