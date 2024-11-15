import typing

import numpy as np
from numpy import cos, sin
from numpy.linalg import inv

from . import Basis, Vector, standard_basis


class BasisTransformation:
    """
    Class representing a basis transformation between two bases, basis1 and basis2. Vectors are transformed using the
    methods self.transform() (basis1 to basis2), and self.inv_transform() (basis2 to basis1)
    Parameters
    ----------
    basis1: :class:`basis`
    basis2: :class:`basis`
    """

    def __init__(self, basis1: Basis, basis2: Basis):

        self.basis1 = basis1
        self.basis2 = basis2
        self.t_matrix = np.dot(basis1.basis, inv(basis2.basis))
        self.invt_matrix = inv(self.t_matrix)

    def transform(self, v: Vector) -> Vector:
        """
        transforms vector v in basis1 coordinates to coordinates in basis2

        Parameters
        ----------
        v: :class:`Vector`

        Returns
        -------
            :class:`Vector`
        """
        _v = np.dot(v.vector, self.t_matrix)
        return Vector(_v, self.basis2)

    def inv_transform(self, v: Vector) -> Vector:
        """
        transforms vector v in basis2 coordinates to coordinates in basis1

        Parameters
        ----------
        v: :class:`Vector`

        Returns
        -------
            :class:`Vector`
        """
        _v = np.dot(v.vector, self.invt_matrix)
        return Vector(_v, self.basis1)


class Rotation:
    """
    class defining a counterclockwise rotation around a given angle and axis in the standard basis, expressed by
    a rotation matrix (from https://en.wikipedia.org/wiki/Rotation_matrix)
    Parameters
    ----------
    angle: float
        angle in rad
    axis: :class:`Vector`
        defining the axis of rotation
    """

    def __init__(self, angle: float, axis: Vector):
        self._angle = angle
        self._axis = axis.global_coord / axis.abs_global
        self._matrix = self._set_matrix()

    def __repr__(self) -> str:
        return (
            f"< rotation around [{self.axis[0]:.2f}, {self.axis[1]:.2f}, {self.axis[2]:.2f}]"
            f" about {self.angle:.2f} rad >"
        )

    @property
    def angle(self) -> float:
        """returns angle"""
        return self._angle

    @property
    def axis(self) -> np.ndarray:
        """returns axis"""
        return self._axis

    @angle.setter
    def angle(self, angle: float):
        """
        setter of rotation angle
        Parameters
        ----------
        angle: float
            angle in rad
        """
        self._angle = angle
        self._set_matrix()

    def rotate(self, point: Vector) -> Vector:
        """
        rotate point around axis about specified angle
        Parameters
        ----------
        point: :class:`Vector`

        Returns
        -------
            :class:`Vector`
                rotated point
        """
        to_std_basis = BasisTransformation(point.basis, standard_basis)
        _new_point = Vector(np.dot(self._matrix, point.global_coord))
        return to_std_basis.inv_transform(_new_point)

    def _set_matrix(self) -> np.ndarray:
        """defines the rotation matrix - auxiliary function"""
        return np.array(
            [
                [
                    cos(self.angle) + self.axis[0] ** 2 * (1 - cos(self.angle)),
                    self.axis[0] * self.axis[1] * (1 - cos(self.angle))
                    - self.axis[2] * sin(self.angle),
                    self.axis[0] * self.axis[2] * (1 - cos(self.angle))
                    + self.axis[1] * sin(self.angle),
                ],
                [
                    self.axis[1] * self.axis[0] * (1 - cos(self.angle))
                    + self.axis[2] * sin(self.angle),
                    cos(self.angle) + self.axis[1] ** 2 * (1 - cos(self.angle)),
                    self.axis[1] * self.axis[2] * (1 - cos(self.angle))
                    - self.axis[0] * sin(self.angle),
                ],
                [
                    self.axis[2] * self.axis[0] * (1 - cos(self.angle))
                    - self.axis[1] * sin(self.angle),
                    self.axis[2] * self.axis[1] * (1 - cos(self.angle))
                    + self.axis[0] * sin(self.angle),
                    cos(self.angle) + self.axis[2] ** 2 * (1 - cos(self.angle)),
                ],
            ]
        )
