from linalg.basis import basis, standard_basis, vector
from cell.contents import atom, molecule, auto_label_atoms, chemical_formula, lattice, cell, super_cell
from cell.generate import lattice_from_cell_parameters, cell_from_crystal, cell_from_cif
from cell.tools import move

from crystals import Crystal
from draw.draw import make_figure
from draw.cells import draw_supercell

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt


cll = cell_from_cif("testdata/erk.cif")
cll.atoms_to_molecule()

scll = super_cell(cll, (3,3,1))
f, ax = make_figure("off", (-2,30), (-2,30), (-2,10))
draw_supercell(ax, scll)
f.tight_layout()
f.show()
