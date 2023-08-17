import numpy as np
from numpy.linalg import det


class basis():
    def __init__(self, v1, v2, v3):
        """
        Generates basis from three linear independent vectors v1, v2, v3
        :param v1: list, tuple or nd.array with three entries as coordinated v11, v12, v13
        :param v2: list, tuple or nd.array with three entries as coordinated v21, v22, v23
        :param v3: list, tuple or nd.array with three entries as coordinated v31, v32, v33
        """
        _vectors = [v1, v2, v3]
        for _v in _vectors:
            if _v is type(np.ndarray):
                pass
            else:
                _v = np.array(_v)

        self.v1, self.v2, self.v3 = v1, v2, v3

        try:
            if self.__check_linear_independece == False:
                raise LinearAlgebraError
        except LinearAlgebraError:
            print("vectors are not linear independent and cannot form a basis")




    def __check_linear_independece(self):
        """
        checks if vectors v1, v2, v3 are linear independent by calculating the determinant
        """
        M = np.array([self.v1,self.v2,self.v3])
        return det(M) == 0




class LinearAlgebraError(Exception):
    pass
