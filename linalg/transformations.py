import numpy as np
from numpy.linalg import inv

from .basis import basis, vector, standard_basis

class basis_transformation:
    def __init__(self, basis1, basis2):
        """

        :param start_basis:
        :param target_basis:
        """
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