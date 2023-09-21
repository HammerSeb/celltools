import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from pyqtgraph import functions as fn
from pyqtgraph.Qt import QtCore
from linalg.basis import basis, vector
from draw.draw import make_figure, GLPoints, draw_basis, draw_frame
from draw.cells import draw_cell, draw_supercell
from cell.contents import super_cell
from cell.generate import cell_from_cif
from cell.tools import move


cll = cell_from_cif("testdata/erk.cif")
cll.atoms_to_molecule()
spcll = super_cell(cll, (5,3,1))
move(spcll.molecules[7],vector([0.5,0,0], spcll.lattice))

w = make_figure()

# _ = draw_frame(w, cll.lattice)
# _ = draw_basis(w, cll.lattice)
test = draw_supercell(w, spcll, lw=1)

test[7].setData(color=[1,0,0,1])

pg.exec()