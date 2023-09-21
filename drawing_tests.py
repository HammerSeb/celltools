import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from pyqtgraph import functions as fn
from pyqtgraph.Qt import QtCore

from draw.draw import make_figure, GLPoints



w = make_figure()


pos = np.array([[0,0,0], [1,0,0], [0,1,0], [1,1,1]])
color = np.array([ [1,1,1,1] for i in range(4) ])
size = np.array([0.4 for i in range(4)])

pts = GLPoints(pos=pos,color=color,size=size)
w.addItem(pts)

lns = gl.GLLinePlotItem(pos=pos, color=color, mode='lines')
w.addItem(lns)

pg.exec()