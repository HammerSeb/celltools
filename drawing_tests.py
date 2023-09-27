import numpy as np
import pyqtgraph as pg
from draw.draw import make_figure, draw_line, draw_plane
from draw.cells import draw_cell, draw_supercell
from cell.generate import cell_from_cif
from cell.contents import super_cell
from linalg.find import average_line
from linalg.basis import plane, vector


cll = cell_from_cif("testdata/erk.cif")
cll.atoms_to_molecule()
cll.molecules[0].auto_bonds(rmin = 0.5, rmax=2.05)
# scll = super_cell(cll, (10,3,5))
# molc_coords = [atm.coords for atm in cll.molecules[0].atoms]
# ln = average_line(molc_coords)

w = make_figure()

_cll = draw_cell(w, cll)
# _scll = draw_supercell(w, scll, lw=.2)
# draw_line(w, ln, range=[-10, 10], lw=5, c=[1, 0, 0, 1])
draw_plane(w, pln, range=[(-10,10), (-10,10)])
pg.exec()