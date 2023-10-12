import typing
from typing import List, Union, Optional
import numpy as np
from itertools import chain

import pyqtgraph as pg
import pyqtgraph.opengl as gl

from celltools.linalg import Basis, Vector, standard_basis


class GLPoints(gl.GLScatterPlotItem):
    def __init__(self, pos: Optional[Union[List[Vector], np.ndarray]]=[], size: Optional[np.ndarray]=[],
                 color: Optional[np.ndarray]=[]):
        """
        container class to draw points, inherits from pyqtgraph.opengl.GLScatterPlotItem. See pyqtgraph manual for
        details.
        Parameters
        ----------
        pos: list of N :class:`Vector` or np.ndarray size (N,3)
            list of points. If point positions are given as array, coordinates need to be expressed in global reference
            frame
        size: np.ndarray size (N,)
            list of point sizes in pixels
        color: np.ndarray size (N, 4)
            list of point colors in rgba
        """
        if pos:
            if type(pos[0]) == type(pos[-1]):
                _pos = []
                if isinstance(pos[0], Vector):
                    for p in pos:
                        _pos.append(p.global_coord)
                else:
                    for p in pos:
                        _pos.append(p)
            else:
                raise TypeError("position list must contain elements of same type")
            super().__init__(pos=pos, size=size, color=color, pxMode=False)
        else:
            super().__init__(pxMode=False)

    def add_point(self, pos: Union[ Vector, np.ndarray], size: float, color: np.ndarray):
        """
        adds point to GLPoints object.
        Parameters
        ----------
        pos: :class:`Vector` or np.ndarray size (3,)
            position of point. If given as array, coordinates need to be expressed in global reference frame
        size: float
            size in pixels
        color: np.ndarray size (4,)
            point color as rgba

        Returns
        -------

        """
        _pos, _size, _color = self.pos, self.size, self.color
        if not np.any(_pos):
            _pos, _size, _color = [], [], []
        else:
            _pos, _size, _color = _pos.tolist(), _size.tolist(), _color.tolist()
        if isinstance(pos, Vector):
            _pos.append(pos.global_coord)
        else:
            _pos.append(pos)

        _size.append(size)
        _color.append(color)
        self.setData(pos=np.array(_pos), size=np.array(_size), color=np.array(_color))


class GLLines(gl.GLLinePlotItem):
    def __init__(self, pos: Optional[Union[List[Vector], np.ndarray]] = [], **kwargs):
        """
        container class to draw lines. Inherits from pyqtgraph.opengl.GLLinePlotItem. See pyqtgraph manual for
        details.
        Parameters
        ----------
        pos: List of N :class:`Vector` or np.ndarray size (N,3); N % 2 = 0
            List of pairs of point between which a line is drawn. List length must be multiple of 2. If given as array,
            coordinates need to be expressed in global reference frame
        kwargs: keyword arguments of GLLinePlotItem
            keep size in mind of arrays in mind. Class methods to set color and size might be easier.
            mode is always set to "lines".
        """
        if pos:
            if len(pos) % 2 == 0:
                if type(pos[0]) == type(pos[-1]):
                    _pos = []
                    if isinstance(pos[0], Vector):
                        for p in pos:
                            _pos.append(p.global_coord)
                    else:
                        for p in pos:
                            _pos.append(p)
                else:
                    raise TypeError("position list must contain elements of same type")
            else:
                raise ValueError("only pairs of coordinates accepted")

            super().__init__(pos=_pos, mode="lines", **kwargs)
        else:
            super().__init__( mode="lines", **kwargs)

    def add_line(self, coord1: Union[Vector, np.ndarray], coord2: Union[Vector, np.ndarray],
                 c: Optional[np.ndarray]=None):
        """
        adds a line between the points at coord1 and coord2 to GLLine object.
        Parameters
        ----------
        coord1: Vector or np.ndarray size (3,)
            position of point 1.  If given as array, coordinates need to be expressed in global reference frame
        coord2: Vector or np.ndarray size (3,)
            position of point 2.  If given as array, coordinates need to be expressed in global reference frame
        c: Optional, np.ndarray size (4,)
            line color as rgba
        """
        if isinstance(coord1, Vector) and isinstance(coord2, Vector):
            _c1, _c2 = coord1.global_coord, coord2.global_coord
        elif type(coord1) != type(coord2):
            raise ValueError("coordinates must be of same type")
        else:
            _c1, _c2 = coord1, coord2

        _pos, _color = self.pos, self.color
        if not np.any(_pos):
            _pos, _color = [], []
        else:
            _pos, _color = _pos.tolist(), _color.tolist()

        _pos.append(_c1)
        _pos.append(_c2)
        self.setData(pos=np.array(_pos))
        if c:
            _color.append(c)
            _color.append(c)
            self.setData(color=np.array(_color))

    def set_linewidth(self, lw: float):
        """
        set line width of lines
        Parameters
        ----------
        lw: float
            line width in pixels
        """
        self.setData(width=lw)


def make_figure(grid: bool=True, distance: float =20, title: str="celltools draw", **kwargs) -> gl.GLViewWidget:
    """
    creates pyqtgraph.opengl.GLViewWidget which is used as a "canvas" for 3D drawing.
    Parameters
    ----------
    grid: bool
        weather or not to show grid on the (x,y)-plane of the global reference system
    distance: float
        sets camera position, i.e. how much space can be seen
    title: str
        window title
    kwargs: GLViewWidget keyword arguments
        distance is set by distance argument of function.

    Returns
    -------
        pyqtgraph.opengl.GLViewWidget
    """
    app = pg.mkQApp("CellTools")
    w = gl.GLViewWidget(**kwargs)
    w.show()
    w.setWindowTitle(title)
    w.setCameraPosition(distance=distance)

    if grid:
        g = gl.GLGridItem()
        w.addItem(g)

    return w

def draw_basis(w, basis: Basis, lw: float=3) -> GLLines:
    """
    draws the three basis vectors with their full length in colors red, green and blue, respectively.
    Parameters
    ----------
    w: pyqtgraph.opengl.GLViewWidget
        view widget, use make_figure function to generate
    basis: :class:`Basis`
        basis system to draw
    lw: float
        line width in pixels

    Returns
    -------
        GLLines
    """
    _origin = Vector(basis.offset)
    _pairs = list(chain(*[[_origin, Vector(v)] for v in basis]))
    lns = GLLines(pos=_pairs)
    lns.set_linewidth(lw)
    lns.setData(color=np.array([[1,0,0,1],[1,0,0,1],[0,1,0,1],[0,1,0,1],[0,0,1,1],[0,0,1,1]]))
    w.addItem(lns)
    return lns

def draw_frame(w, basis, lw=1.5):
    _origin = Vector(basis.offset)
    _a, _b, _c = Vector(basis.basis[0]), Vector(basis.basis[1]), Vector(basis.basis[2])
    _pairs = list(chain(
        *[[_origin,_origin+_a], [_origin,_origin+_b], [_origin+_a, _origin+_a+_b], [_origin+_b, _origin+_a+_b], #bottom
        [_origin, _origin+_c], [_origin+_a, _origin+_a+_c], [_origin+_b, _origin+_c+_b], [_origin+_a+_b, _origin+_a+_b+_c], #sides
        [_origin+_c, _origin+_a+_c], [_origin+_c,_origin+_b+_c], [_origin+_a+_c,_origin+_a+_b+_c], [_origin+_b+_c, _origin+_a+_b+_c] #top
        ]
    ))
    lns = GLLines(pos=_pairs)
    lns.set_linewidth(lw)
    w.addItem(lns)
    return lns


def draw_line(w, ln, range = (-1,1), lw=1.5, c=[1,1,1,1]):
    _pos = np.array([ln.point_on_line(range[0]).global_coord, ln.point_on_line(range[1]).global_coord])
    _line = gl.GLLinePlotItem(pos=_pos, width=lw, color=c)
    w.addItem(_line)

def draw_plane(w, pln, range=[(-5,5),(-5,5)], c=[1,0,0,0.7]):
    _corners = np.array([range[0][0]*pln.basis[0], range[0][1]*pln.basis[0], range[1][0]*pln.basis[1], range[1][1]*pln.basis[1]])
    _face = np.array([[0,1,2],[3,1,0]])
    _colors = np.array([c,c])
    _plane = gl.GLMeshItem(vertexes=_corners, faces=_face, faceColors=_colors, shader="balloon", smooth=True, glOptions="additive")
    w.addItem(_plane)
