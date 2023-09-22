import numpy as np
import pyqtgraph as pg
from draw.draw import make_figure, draw_line, draw_plane
from draw.cells import draw_cell
from cell.generate import cell_from_cif
from linalg.find import average_line
from linalg.basis import plane, vector


cll = cell_from_cif("testdata/erk.cif")
cll.atoms_to_molecule()

molc_coords = [atm.coords for atm in cll.molecules[0].atoms]
ln = average_line(molc_coords)
pln = plane(vector(cll.lattice[1])*0.5, vector(cll.lattice[0]))

w = make_figure()

_cll = draw_cell(w, cll)
draw_line(w, ln, range=[-10, 10], lw=5, c=[1, 0, 0, 1])
draw_plane(w, pln)
pg.exec()