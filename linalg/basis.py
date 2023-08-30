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

class basis:
    #TODO: construction of basis from three vector instances should be possible
    @check_linear_independence()
    def __init__(self, v1, v2, v3, offset=None):
        """
        Generates basis from three linear independent vectors v1, v2, v3. Offset can be given to
        :param v1: list, tuple or nd.array with three entries as coordinates v11, v12, v13
        :param v2: list, tuple or nd.array with three entries as coordinates v21, v22, v23
        :param v3: list, tuple or nd.array with three entries as coordinates v31, v32, v33
        :param offset:  list, tuple or nd.array with three entries as coordinates
        """
        self._basis = np.array([v1, v2, v3])

        if offset == None:
            self._offset = None
        else:
            self.offset(offset)

    @property
    def basis(self):
        return self._basis

    @property
    def offset(self):
        if self._offset:
            return self._offset
        else:
            return np.zeros((3,))

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

    @offset.setter
    def offset(self,offset):
        if isinstance(offset, np.ndarray):
            self._offset = offset
        else:
            self._offset = np.array(offset)

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
            self.basis = np.array([self.basis[order[0]-1], self.basis[order[1]-1], self.basis[order[2]-1]])

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

    def __add__(self, other):
        """addition of 2 vector instances"""
        if not np.all(self.basis == other.basis):
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
        if not np.all(self.basis == other.basis):
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




