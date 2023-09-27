from copy import copy
import numpy as np
from numpy.linalg import det

class LinearAlgebraError(Exception):
    pass

def check_linear_independence():
    """
    Decorator function to check if the three vectors passed to the target function are linear independent
    (assumes to be wrapper for basis class init function)
    :raises:
    LinearAlgebraError if the vectors are not linear independent
    """
    def check_accepts(function):
        def new_function(*args, **kwargs):
            basis = np.array([*args[1:]])
            if det(basis)==0:
                raise LinearAlgebraError("vectors are not linear independent")

            return function(*args, **kwargs)

        new_function.__name__ = function.__name__

        return new_function
    return check_accepts


def parallel(vec1, vec2):
    """checks if two vectors are parallel in standard coordinates"""
    if abs(np.dot(vec1.global_coord, vec2.global_coord)) == vec1.abs_global*vec2.abs_global:
        return True
    else:
        False

class basis:
    @check_linear_independence()
    def __init__(self, v1, v2, v3, offset=np.zeros((3,))):
        """
        Generates basis from three linear independent vectors v1, v2, v3. Offset can be given to
        :param v1: list, tuple or nd.array with three entries as coordinates v11, v12, v13
        :param v2: list, tuple or nd.array with three entries as coordinates v21, v22, v23
        :param v3: list, tuple or nd.array with three entries as coordinates v31, v32, v33
        :param offset:  list, tuple or nd.array with three entries as coordinates
        """
        self._basis = np.array([v1, v2, v3])

        self._offset = offset

    def __repr__(self):
        return f'<basis {self.basis[0]}, {self.basis[1]}, {self.basis[2]}>'

    def __eq__(self, other):
        if isinstance(other, basis):
            if np.all(self.v1 == other.v1) and np.all(self.v2 == other.v2) and np.all(self.v3 == other.v3):
                return True
            else:
                return False
        else:
            raise TypeError("must be compared to instance of basis")

    @property
    def basis(self):
        return self._basis

    @property
    def offset(self):
        return self._offset

    @offset.setter
    def offset(self, offset):
        if isinstance(offset, np.ndarray):
            self._offset = offset
        else:
            self._offset = np.array(offset)

    @property
    def v1(self):
        """return v1 in global coordinates"""
        return self.basis[0] + self.offset

    @property
    def v2(self):
        """return v2 in global coordinates"""
        return self.basis[1] + self.offset

    @property
    def v3(self):
        """return v3 in global coordinates"""
        return self.basis[2] + self.offset

    def permute(self, order):
        """
        permutes the basis vectors specified by a 3-tuple of (new position v1, new position v2, new position v3) given
        by numbers 1, 2 and 3.
        :param order: tuple of size 3
        """
        if np.sum(order) != 6:
            raise ValueError("order must be given in integers of 1,2 and 3 only")
            return
        else:
            self._basis = np.array([self._basis[order[0]-1], self._basis[order[1]-1], self._basis[order[2]-1]])

    def __getitem__(self, index):
        if not isinstance(index, int):
            raise ValueError("index needs to be integer")
            return
        elif abs(index) > 2 and index != -3:
            raise IndexError("index out of bounds")
            return
        else:
            if index == 0 or index == -3:
                return self.basis[0]
            elif index == 1 or index == -2:
                return self.basis[1]
            elif index ==2 or index == -1:
                return self.basis[2]

standard_basis = basis([1,0,0], [0,1,0], [0,0,1])

class vector:
    def __init__(self, vector, basis=standard_basis):
        """
        :param vector: list, tuple, ndarray shape (3,)
        :param basis: instance of basis class
        """
        self.vector = np.array([vector[0], vector[1], vector[2]])
        self.basis = basis

    def __repr__(self):
        return f"{self.vector}"

    def __add__(self, other):
        """addition of 2 vector instances"""
        if not self.basis == other.basis:
            raise LinearAlgebraError("basis do not match")
            return
        else:
            return vector([self.vector[0]+other.vector[0],
                           self.vector[1]+other.vector[1],
                           self.vector[2]+other.vector[2]],
                          self.basis)
        pass
    def __sub__(self, other):
        """substraction of 2 vector instances"""
        if not self.basis == other.basis:
            raise LinearAlgebraError("basis do not match")
            return
        else:
            return vector([self.vector[0]-other.vector[0],
                           self.vector[1]-other.vector[1],
                           self.vector[2]-other.vector[2]],
                          self.basis)
        pass
    def __mul__(self, other):
        """multiplication with a scalar"""
        if isinstance(other, float) or isinstance(other, int):
            return vector([self.vector[0]*other,
                           self.vector[1]*other,
                           self.vector[2]*other],
                          self.basis)
        else:
            raise ValueError("only scalar multiplication with real numbers allowed")
            return
        pass

    def __rmul__(self, other):
        """multiplication with a scalar"""
        if isinstance(other, float) or isinstance(other, int):
            return vector([self.vector[0]*other,
                           self.vector[1]*other,
                           self.vector[2]*other],
                          self.basis)
        else:
            raise ValueError("only scalar multiplication with real numbers allowed")
            return
        pass

    def __eq__(self, other):
        if isinstance(other, vector):
            return np.all(self.global_coord == other.global_coord)
        else:
            raise TypeError("must be compared to vector instance")

    def __getitem__(self, index):
        if not isinstance(index, int):
            raise ValueError("index needs to be integer")
        elif abs(index) > 2 and index != -3:
            raise IndexError("index out of bounds")
        else:
            if index == 0 or index == -3:
                return self.vector[0]
            elif index == 1 or index == -2:
                return self.vector[1]
            elif index ==2 or index == -1:
                return self.vector[2]

    @property
    def abs(self):
        """returns absolute value of vector instance in its basis coordinates"""
        return np.sqrt(self.vector[0]**2 + self.vector[1]**2 + self.vector[2]**2)

    @property
    def global_coord(self):
        """return vector in global coordinates"""
        return self.vector[0]*self.basis[0] + self.vector[1]*self.basis[1] + self.vector[2]*self.basis[2]

    @property
    def abs_global(self):
        return np.sqrt(self.global_coord[0]**2 + self.global_coord[1]**2 + self.global_coord[2]**2)

    def parallel(self, other):
        if parallel(self, other):
            return True
        else:
            False

class line:
    def __init__(self, origin, direction):
        if origin.basis == direction.basis:
            self._origin = origin
            self._direction = direction
            self._basis = self.origin.basis
        else:
            raise LinearAlgebraError("input vectors need to be expressed in the same basis")

    def __repr__(self):
        return f"< line: {self.origin} + k*{self.direction} >"

    def __eq__(self, other):
        if isinstance(other, line):
            if self.on_line(other.origin) and self.parallel(other):
                return True
            else:
                return False
        else:
            raise ValueError("must be compared to line instance")

    @property
    def origin(self):
        return self._origin

    @property
    def direction(self):
        return self._direction

    def on_line(self, point):
        if  ((point.global_coord[0] - self.origin.global_coord[0]) / self.direction.global_coord[0]) == \
            ((point.global_coord[1] - self.origin.global_coord[1]) / self.direction.global_coord[1]) and \
            ((point.global_coord[1] - self.origin.global_coord[1]) / self.direction.global_coord[1]) == \
            ((point.global_coord[2] - self.origin.global_coord[2]) / self.direction.global_coord[2]) and \
            ((point.global_coord[2] - self.origin.global_coord[2]) / self.direction.global_coord[2]) == \
            ((point.global_coord[0] - self.origin.global_coord[0]) / self.direction.global_coord[0]):
            return True
        else:
            return False

    def parallel(self, line):
        if parallel(self.direction, line.direction):
            return True
        else:
            False

    def distance(self, point):
        """
        distance from line to point in standard basis
        formula adapted from https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        Parameters
        ----------
        point: :class:`vector`

        Returns
        -------
            float
        """
        return (( vector(self.origin.global_coord-point.global_coord).abs**2 * self.direction.abs_global**2
                  - np.dot((self.origin.global_coord-point.global_coord), self.direction.global_coord)**2)
                / self.direction.abs_global**2)

    def point_on_line(self,k):
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

class plane:
    def __init__(self, origin, normal):
        if origin.basis == normal.basis:
            self._origin = origin
            self._normal = normal
            self._define_basis()
        else:
            raise LinearAlgebraError("input vectors need to be expressed in the same basis")

    def __repr__(self):
        return (f"< plane: {self.parametric_form[0]:.2f}*x + {self.parametric_form[1]:.2f}*y "
                f"+ {self.parametric_form[2]:.2f}*z = {self.parametric_form[3]} >")

    @property
    def normal(self):
        return self._normal

    @property
    def origin(self):
        return self._origin

    @property
    def basis(self):
        return self._basis

    @basis.setter
    def basis(self, basis):
        if isinstance(basis,basis):
            self._basis = basis
        else:
            raise ValueError("input needs to be of type basis")

    @property
    def parametric_form(self):
        """
        gives an array containing the for parameters of the normalized parametric form n1*x + n2*y +n3*z = d
        in the standard basis. ni are the coordinates of the normalized normal vector on the plane

        Returns
        -------
            nd.array
                [n1, n2, n3, d]
        """
        return np.array([self.normal.global_coord[0]/self.normal.abs_global,
                         self.normal.global_coord[1]/self.normal.abs_global,
                         self.normal.global_coord[2]/self.normal.abs_global,
                         np.dot(self.normal.global_coord//self.normal.abs_global, self.origin.global_coord)])

    def on_plane(self, point):
        if np.dot(self.normal.global_coord, point.global_coord) == self.parametric_form[3]:
            return True
        else:
            return False

    def _define_basis(self):
        """
        defines orthonormal basis of the plane centered at origin
            TODO: normal vector must always be 3rd basis vector
        """
        # parallel to xy, xz, yz plane
        if np.count_nonzero(self.parametric_form[:-1]==0) == 2:
            self._basis = copy(standard_basis)
            if not self.parametric_form[0]==0:
                self._basis.permute((2,3,1))
            if not self.parametric_form[1]==0:
                self._basis.permute((1,3,2))


            self._basis.offset = self.origin.global_coord
        # parallel to x,y or z axis
        elif np.count_nonzero(self.parametric_form[:-1]==0) == 1:
            if self.parametric_form[0]==0: #x-axis
                self._basis = basis(
                 [0,np.sqrt(2),np.sqrt(2)], [0,-1*np.sqrt(2),np.sqrt(2)], [1,0,0]
                )

            elif self.parametric_form[1]==0: #y-axis
                self._basis = basis(
                    [np.sqrt(2),0,np.sqrt(2)], [np.sqrt(2),0,-1*np.sqrt(2)], [0,1,0]
                )

            elif self.parametric_form[2]==0: #z-axis
                self._basis = basis(
                    [np.sqrt(2),np.sqrt(2),0], [-1*np.sqrt(2), np.sqrt(2), 0], [0,0,1],
                )

            self._basis.offset = self.origin.global_coord
            #plane with no restrictions
        else:
            _origin = vector(self.origin.global_coord)
            _x_intsec = vector([self.parametric_form[3]/self.parametric_form[0],0,0])
            _y_intsec = vector([0,self.parametric_form[3]/self.parametric_form[1],0])
            _basis3 = self.normal.global_coord/self.normal.abs_global
            if not _origin == _x_intsec:
                _basis1 = ((_x_intsec - _origin)* (1 / vector(_x_intsec - _origin).abs_global)).global_coord
            else:
                _basis1 = ((_y_intsec - _origin) * (1 / vector(_y_intsec - _origin).abs_global)).global_coord

            _basis2 = (vector(np.cross(_basis1, _basis3))*(1/vector(np.cross(_basis1, _basis3)).abs_global)).global_coord

            self._basis = basis(_basis1, _basis2, _basis3)
            self.basis.offset = self.origin.global_coord