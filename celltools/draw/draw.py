import numpy as np
from itertools import chain

import pyqtgraph as pg
import pyqtgraph.opengl as gl
from pyqtgraph import functions as fn
from pyqtgraph.Qt import QtCore

from celltools.linalg.basis import basis, vector, standard_basis


class GLPoints(gl.GLScatterPlotItem):
    def __init__(self, pos=[], size=[], color=[]):
        if pos:
            if type(pos[0]) == type(pos[-1]):
                _pos = []
                if isinstance(pos[0], vector):
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

    def add_point(self, pos, size, color):
        _pos, _size, _color = self.pos, self.size, self.color
        if not np.any(_pos):
            _pos, _size, _color = [], [], []
        else:
            _pos, _size, _color = _pos.tolist(), _size.tolist(), _color.tolist()
        if isinstance(pos, vector):
            _pos.append(pos.global_coord)
        else:
            _pos.append(pos)

        _size.append(size)
        _color.append(color)
        self.setData(pos=np.array(_pos), size=np.array(_size), color=np.array(_color))

class GLLines(gl.GLLinePlotItem):
    def __init__(self, pos=[]):
        if pos:
            if len(pos) % 2 == 0:
                if type(pos[0]) == type(pos[-1]):
                    _pos = []
                    if isinstance(pos[0], vector):
                        for p in pos:
                            _pos.append(p.global_coord)
                    else:
                        for p in pos:
                            _pos.append(p)
                else:
                    raise TypeError("position list must contain elements of same type")
            else:
                raise ValueError("only pairs of coordinates accepted")

            super().__init__(pos=_pos, mode="lines")
        else:
            super().__init__( mode="lines")

    def add_line(self, coord1, coord2, c=None):
        if isinstance(coord1, vector) and isinstance(coord2, vector):
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

    def set_linewidth(self,lw):
        self.setData(width=lw)

def make_figure(grid=True, distance=20, title="celltools draw"):
    """
    creates a pyqtgraph opengl widget which can be used for 3D drawing
    Parameters
    ----------
    grid
    distance
    title

    Returns
    -------
    widget
    """
    app = pg.mkQApp("CellTools")
    w = gl.GLViewWidget()
    w.show()
    w.setWindowTitle(title)
    w.setCameraPosition(distance=distance)

    if grid:
        g = gl.GLGridItem()
        w.addItem(g)

    return w

def draw_basis(w, basis, lw=3):
    _origin = vector(basis.offset)
    _pairs = list(chain(*[[_origin, vector(v)] for v in basis]))
    lns = GLLines(pos=_pairs)
    lns.set_linewidth(lw)
    lns.setData(color=np.array([[1,0,0,1],[1,0,0,1],[0,1,0,1],[0,1,0,1],[0,0,1,1],[0,0,1,1]]))
    w.addItem(lns)
    return lns

def draw_frame(w, basis, lw=1.5):
    _origin = vector(basis.offset)
    _a, _b, _c = vector(basis.basis[0]), vector(basis.basis[1]), vector(basis.basis[2])
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
