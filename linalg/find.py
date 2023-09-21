import numpy as np
from itertools import combinations
from linalg.basis import basis, vector, standard_basis, line, plane



#TODO: function that finds avergage plane through a list of points
#TODO: function that finds avergage line through list of points



# parametric plane n1*x + n2*y + n3*z = d
# there is of course no parametric line equation in 3 dimensions

# I think the easiest way to do it to calculate all possible lines/planes through all combinations of 2/3 points of a
# given point list and then calcualte the average line/plane from that


def _average_vectors(list_of_vectors):
    avg_vec = list_of_vectors[0]
    for vec in list_of_vectors[1:]:
        avg_vec += vec
    return avg_vec*(1/len(list_of_vectors))





def average_line(list_of_points):
    """
    returns average through a given set of points. If all are given in the same basis, the output line will be given in
    the same basis, otherwise the standard cartesian basis is used.
    Parameters
    ----------
    list_of_points: list of vectors
        list of :class:`vector`, not necessarly the same basis

    Returns
    -------
    avg_line: :class:`line`
        average line through points

    """
    point_pairs = list(combinations(list_of_points, 2))
    _same_basis, _origin, _direction = [], [], []
    for pair in point_pairs:
        _same_basis.append(pair[0].basis == pair[1].basis)

    if np.all(_same_basis) == True:
        for pair in point_pairs:
            _origin.append(pair[0])
            _direction.append(pair[1]-pair[0])
        return line(_average_vectors(_origin), _average_vectors(_direction))
    else:
        for pair in point_pairs:
            _origin.append(vector(pair[0].global_coord))
            _direction.append(vector(pair[1].global_coord)-vector(pair[0].global_coord))
        return line(_average_vectors(_origin), _average_vectors(_direction))
