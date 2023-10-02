import typing

import numpy as np
from numpy import cos, sin
from numpy.linalg import inv

from .basis import basis, vector, standard_basis, line

class basis_transformation:
    """
    Sets up a basis transformation between basis1 and basis2
    Parameters
    ----------
    basis1: :class:`basis`
    basis2: :class:`basis`
    """
    def __init__(self, basis1, basis2):

        self.basis1 = basis1
        self.basis2 = basis2
        self.t_matrix = np.dot(basis1.basis, inv(basis2.basis))
        self.invt_matrix = inv(self.t_matrix)


    def transform(self, v):
        """
        transforms vector v in basis1 coordinates to coordinates in basis2
        :param v: nd.array of shape (3,)
        :return: nd.array of shape (3,)
        """
        _v = np.dot(v.vector, self.t_matrix)
        return vector(_v, self.basis2)

    def inv_transform(self, v):
        """
        transforms vector v in basis2 coordinates to coordinates in basis1
        :param v: nd.array of shape (3,)
        :return: nd.array of shape (3,)
        """
        _v = np.dot(v.vector, self.invt_matrix)
        return vector(_v, self.basis1)

class rotation:
    def __init__(self, angle: float, axis: vector):
        """
        class defining a counterclockwise rotation around a given angle and axis in the standard basis, expressed by
        a rotation matrix (from https://en.wikipedia.org/wiki/Rotation_matrix)
        Parameters
        ----------
        angle: float
        axis: :class:`vector`
        """
        self._angle = angle
        self._axis = axis.global_coord / axis.abs_global
        self._matrix = self._set_matrix()


    def __repr__(self):
        return (f"< rotation around [{self.axis[0]:.2f}, {self.axis[1]:.2f}, {self.axis[2]:.2f}]"
                f" about {self.angle:.2f} rad >")

    @property
    def angle(self) -> float:
        return self._angle

    @property
    def axis(self) -> np.ndarray:
        return self._axis

    @angle.setter
    def angle(self, angle):
        self._angle = angle
        self._set_matrix()

    def rotate(self, point: vector) -> vector:
        """rotate point around axis about specified angle"""
        to_std_basis = basis_transformation(point.basis, standard_basis)
        _new_point = vector(np.dot(self._matrix, point.global_coord))
        return to_std_basis.inv_transform(_new_point)

    def _set_matrix(self) -> np.ndarray:
        return np.array(
            [
            [cos(self.angle) + self.axis[0]**2 * (1-cos(self.angle)) ,
             self.axis[0]*self.axis[1] * (1-cos(self.angle)) - self.axis[2] * sin(self.angle),
             self.axis[0]*self.axis[2] * (1-cos(self.angle)) + self.axis[1] * sin(self.angle)],
            [self.axis[1]*self.axis[0] * (1-cos(self.angle)) + self.axis[2] * sin(self.angle),
             cos(self.angle) + self.axis[1] ** 2 * (1 - cos(self.angle)),
             self.axis[1] * self.axis[2] * (1 - cos(self.angle)) - self.axis[0] * sin(self.angle)],
            [self.axis[2] * self.axis[0] * (1 - cos(self.angle)) - self.axis[1] * sin(self.angle),
             self.axis[2] * self.axis[1] * (1 - cos(self.angle)) + self.axis[0] * sin(self.angle),
             cos(self.angle) + self.axis[2] ** 2 * (1 - cos(self.angle))]
            ]
        )

