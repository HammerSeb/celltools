from copy import copy
import numpy as np
import pyqtgraph as pg
from celltools.draw.draw import make_figure, draw_line, draw_plane
from celltools.draw.cells import draw_cell, draw_supercell
from celltools.cell.generate import cell_from_cif
from celltools.cell.tools import super_cell, rotate
from celltools.linalg.find import average_line, average_plane
from celltools.linalg.basis import plane, vector, line


cll = cell_from_cif("testdata/erk.cif")
cll.atoms_to_molecule()
cll.molecules[0].auto_bonds(rmin = 0.5, rmax=2.05)
# scll = super_cell(cll, (10,3,5))
molc_coords = [atm.coords for atm in cll.molecules[0].atoms]
pln = average_plane(molc_coords)
lne = line(pln.origin, pln.normal)

cll2 = copy(cll)

print(pln)

# w1 = make_figure(title="original cell", grid=False)
w2 = make_figure(title="rotated molecule", grid=False)

# _cll1 = draw_cell(w1, cll)
# # _scll = draw_supercell(w, scll, lw=.2)
# draw_plane(w1, pln, range=[(-10,10), (-10,10)])
# draw_line(w1, lne, (-10, 10), 3, [1,0,0,1])


rotate(cll2.molecules[0], lne, 45)
_cll2 = draw_cell(w2, cll2)
draw_plane(w2, pln, range=[(-10,10), (-10,10)])
draw_line(w2, lne, (-10, 10), 3, [1,0,0,1])
pg.exec()