import numpy as np
from itertools import combinations
from linalg.basis import basis, vector, standard_basis, line, plane
from scipy.optimize import minimize



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
    returns average line through a given set of points in the standard basis. The avergage line is determined by a
    least square minimization of the distance of all points to the line.

    Parameters
    ----------
    list_of_points: list of vectors
        list of :class:`vector`, not necessarly the same basis

    Returns
    -------
    avg_line: :class:`line`
        average line through points

    """
    def _minimize(x, pnts):
        """
        minimization function: giving the sum of all distances squared to the line defined by x
        Parameters
        ----------
        x: nd.array
            parameters defining the line's origin and direction [o1, o2, o3, d1, d2, d3]
        pnts: list of vectors
              list of :class:`vector`, not necessarly the same basis

        Returns
        -------
            float
                sum of squared distances to line
        """
        _d = []
        _line = line(vector([x[0], x[1], x[2]]), vector([x[3], x[4], x[5]]))
        for pnt in pnts:
            _d.append(_line.distance(pnt))

        return np.sum(np.square(np.array(_d)))

    # estimating starting values as averages
    point_pairs = list(combinations(list_of_points, 2))
    _origin, _direction = [], []
    for pair in point_pairs:
        _origin.append(vector(pair[0].global_coord))
        _direction.append(vector(pair[1].global_coord)-vector(pair[0].global_coord))
    x0 = np.array([_average_vectors(_origin).global_coord, _average_vectors(_direction).global_coord]).flatten()

    res = minimize(_minimize, x0, args=(list_of_points))

    return line(vector([res.x[0],res.x[1],res.x[2]]),
                vector([res.x[3],res.x[4],res.x[5]])*(1/vector([res.x[3],res.x[4],res.x[5]]).abs))


def average_plane(list_of_points):
    """
      returns average plane through a given set of points in the standard basis. The average plane is determined by a
      least square minimization of the distance of all points to the line.

      Parameters
      ----------
      list_of_points: list of vectors
          list of :class:`vector`, not necessarly the same basis

      Returns
      -------
      avg_plane: :class:`plane`
          average plane through points

      """

    def _minimize(x, pnts):
        """
        minimization function: giving the sum of all distances squared to the line defined by x
        Parameters
        ----------
        x: nd.array
            parameters defining the line's origin and direction [o1, o2, o3, d1, d2, d3]
        pnts: list of vectors
              list of :class:`vector`, not necessarly the same basis

        Returns
        -------
            float
                sum of squared distances to line
        """
        _d = []
        _plane = plane(vector([x[0], x[1], x[2]]), vector([x[3], x[4], x[5]]))
        for pnt in pnts:
            _d.append(_plane.distance(pnt))

        return np.sum(np.square(np.array(_d)))

    # estimating starting values from average
    _avg_vec = _average_vectors(list_of_points)
    _origin0 = _avg_vec.global_coord
    _norm0 = np.cross((list_of_points[0]-_avg_vec).global_coord,(list_of_points[-1]-_avg_vec).global_coord)
    x0 = np.array([_origin0, _norm0]).flatten()
    res = minimize(_minimize, x0, args=(list_of_points))

    return plane(vector([res.x[0], res.x[1], res.x[2]]),
                 vector([res.x[3], res.x[4], res.x[5]])*(1/vector([res.x[3], res.x[4], res.x[5]]).abs_global))