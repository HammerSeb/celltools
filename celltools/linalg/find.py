from itertools import combinations
from typing import List

import numpy as np
from scipy.optimize import minimize

from . import Vector, Line, Plane


def _average_vectors(list_of_vectors: List[Vector]) -> Vector:
    """
    auxillary function returning an arthimatic mean vector from a list of vectors
    Parameters
    ----------
    list_of_vectors: list of :class:`Vector`

    Returns
    -------
        :class:`Vector`
    """
    avg_vec = list_of_vectors[0]
    for vec in list_of_vectors[1:]:
        avg_vec += vec
    return avg_vec * (1 / len(list_of_vectors))


def average_line(list_of_points: List[Vector]) -> Line:
    """
    returns average line through a given set of points. The average line is determined by a
    least square minimization of the distance of all points to the line.

    Parameters
    ----------
    list_of_points: list of :class:`vectors`
        list of vectors, not necessarily the same basis, defining points in the global reference frame

    Returns
    -------
    avg_line: :class:`Line`
        average line through points

    """

    def _minimize(x: np.ndarray, pnts: List[Vector]) -> float:
        """
        minimization function: giving the sum of all distances squared to the line defined by x
        Parameters
        ----------
        x: nd.array
            parameters defining the line's origin and direction [o1, o2, o3, d1, d2, d3]
        pnts: list of vectors
              list of :class:`vector`, not necessarily the same basis

        Returns
        -------
            float
                sum of squared distances to line
        """
        _d = []
        _line = Line(Vector([x[0], x[1], x[2]]), Vector([x[3], x[4], x[5]]))
        for pnt in pnts:
            _d.append(_line.distance(pnt))

        return np.sum(np.square(np.array(_d)))

    # estimating starting values as averages
    point_pairs = list(combinations(list_of_points, 2))
    _origin, _direction = [], []
    for pair in point_pairs:
        _origin.append(Vector(pair[0].global_coord))
        _direction.append(Vector(pair[1].global_coord) - Vector(pair[0].global_coord))
    x0 = np.array([_average_vectors(_origin).global_coord, _average_vectors(_direction).global_coord]).flatten()

    res = minimize(_minimize, x0, args=(list_of_points))

    return Line(Vector([res.x[0], res.x[1], res.x[2]]),
                Vector([res.x[3], res.x[4], res.x[5]]) * (1 / Vector([res.x[3], res.x[4], res.x[5]]).abs))


def average_plane(list_of_points: List[Vector]) -> Plane:
    """
      returns average plane through a given set of points. The average plane is determined by a
      least square minimization of the distance of all points to the line.

      Parameters
      ----------
     list_of_points: list of :class:`vectors`
        list of vectors, not necessarily the same basis, defining points in the global reference frame

      Returns
      -------
      avg_plane: :class:`plane`
          average plane through points

      """

    def _minimize(x: np.ndarray, pnts: List[Vector]) -> float:
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
        _plane = Plane(Vector([x[0], x[1], x[2]]), Vector([x[3], x[4], x[5]]))
        for pnt in pnts:
            _d.append(_plane.distance(pnt))

        return np.sum(np.square(np.array(_d)))

    # estimating starting values from average
    _avg_vec = _average_vectors(list_of_points)
    _origin0 = _avg_vec.global_coord
    _norm0 = np.cross((list_of_points[0] - _avg_vec).global_coord, (list_of_points[-1] - _avg_vec).global_coord)
    x0 = np.array([_origin0, _norm0]).flatten()
    res = minimize(_minimize, x0, args=(list_of_points))

    return Plane(Vector([res.x[0], res.x[1], res.x[2]]),
                 Vector([res.x[3], res.x[4], res.x[5]]) * (1 / Vector([res.x[3], res.x[4], res.x[5]]).abs_global))
