import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from crystals import Crystal


from linalg.basis import vector, basis, standard_basis
from draw.draw import make_figure, draw_basis
from draw.cells import draw_atom, draw_molecule, draw_cell, draw_supercell

from cell.contents import atom, molecule, cell, super_cell
from cell.generate import lattice_from_crystal

cryst = Crystal.from_cif("testdata/erk.cif")
uc = lattice_from_crystal(cryst)


f, ax = make_figure(axis="off",xlim=(-2,45),ylim=(-2,45),zlim=(-2,45))
draw_cell(ax, uc)
f.show()