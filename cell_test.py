from linalg.basis import basis, standard_basis, vector
from cell.contents import atom, molecule, auto_label_atoms, chemical_formula, lattice, cell, super_cell
from cell.generate import lattice_from_cell_parameters
from cell.tools import move

from draw.draw import make_figure, draw_basis

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt

# a = atom('C', vector([0,0,0]))
# atms = [atom(el, pos) for el, pos in zip(['C', 'C', 'N', 'N', 'N', 'H', 'C'], [vector([x,0,0]) for x in range(7)])]
# mol = molecule(atms)
# print(mol.atoms)
# move(mol, vector([1,0,1]))
# print(mol.atoms)

latt = lattice_from_cell_parameters(3.8052, 12.9590, 12.0430, 90.6400, 95.2600, 90.7200)

atms = [atom('C', vector([0,0,0], latt)), atom('O', vector([0.5,.5,.5], latt)), atom('H', vector([0.1,0,0], latt))]
molc = molecule(atms)
uc = cell(latt, [molc])
suc = super_cell(uc, (3,3,3))
f, ax = make_figure(axis="off",xlim=(-2,15),ylim=(-2,15),zlim=(-2,15))
draw_basis(ax, latt)
f.show()