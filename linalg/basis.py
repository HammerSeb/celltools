import numpy as np
from numpy.linalg import det

class LinearAlgebraError(Exception):
    pass

def check_linear_independence():
    """
    Decorator function to check if the three vectors passed to the target function are linear independent
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
    @check_linear_independence()
    def __init__(self, v1, v2, v3):
        """
        Generates basis from three linear independent vectors v1, v2, v3
        :param v1: list, tuple or nd.array with three entries as coordinated v11, v12, v13
        :param v2: list, tuple or nd.array with three entries as coordinated v21, v22, v23
        :param v3: list, tuple or nd.array with three entries as coordinated v31, v32, v33
        """
        self.basis = np.array([v1, v2, v3])








