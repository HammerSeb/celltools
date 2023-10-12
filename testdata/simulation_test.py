from copy import deepcopy
from itertools import product
import numpy as np
from simulation.ediff import powder_pattern, diffraction_from_cell, diffraction_from_supercell
from celltools.cell.generate import cell_from_cif, cell_to_crystal
from celltools.cell.tools import rotate, move, basis_transformation, super_cell
from celltools.linalg.basis import vector, line, standard_basis
from celltools.linalg.find import average_plane

from celltools.draw.draw import make_figure
from celltools.draw.cells import draw_cell

from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib

import pyqtgraph as pg


matplotlib.use("Qt5Agg")

cell = cell_from_cif("erk.cif")
cell.atoms_to_molecule()
cell.molecules[0].auto_bonds()

scell = super_cell(cell, (3,3,3))
scell_norm = super_cell(cell, (1,1,1))

# #move atom along a-axis 0.4 Angstrom
# cell_trans = deepcopy(cell)
# trans_vec = vector([1,0,0], cell_trans.lattice)
# trans_vec *= 0.4*(1/trans_vec.abs_global)
# move(cell_trans.molecules[0], trans_vec)
#
#
# #rotate 45 deg around normal through ZnPc atom
# cell_rot = deepcopy(cell)
# molc_coords = [atm.coords for atm in cell_rot.molecules[0].atoms]
# ZnPc_plane = average_plane(molc_coords)
# std_to_lattice = basis_transformation(standard_basis, cell_rot.lattice)
# rotation_axis = line(vector([0,0,0], cell_rot.lattice), std_to_lattice.transform(ZnPc_plane.normal))
# rotate(cell_rot.molecules[0], rotation_axis, 15)
#
# # w = make_figure(grid=False)
# # draw_cell(w, cell_rot)
# # pg.exec()
#
#
# # hkl = [(0,k+1,0) for k in range(4)]
#
hkl = list(product(
    [i-5 for i in range(11)],
    [0],
    [i-5 for i in range(11)]
    ))
hkl = hkl[int((len(hkl)-1)/2):]

_vec_sc, _amp_sc = diffraction_from_supercell(hkl, scell)
_vec, _S = diffraction_from_supercell(hkl, scell_norm)
# _vec, _S = diffraction_from_cell(hkl, cell)
_vec_sc = [np.sqrt(v[0]**2+v[2]**2) for v in _vec_sc]
_vec = [np.sqrt(v[0]**2+v[2]**2) for v in _vec]
q = np.linspace(0.3, 6, 500)
Ic = powder_pattern(q, _vec, _S, w=0.05)
Isc = powder_pattern(q, _vec_sc, _amp_sc, w=0.05)


f, ax = plt.subplots(2,1, sharex=True, figsize=(11, 8.5))
ax[0].plot(q, Ic, c="C0", label="unit cell")
ax[1].plot(q, Isc, c="C2", label="super cell")
ax[1].set_xlim(0.3,6)
ax[1].set_xlabel("q [A-1]")
for a in ax:
    # a.set_yticks([], [])
    # a.set_ylabel("delta I [norm]")
    a.legend()
f.suptitle("diffraction patterns")
f.tight_layout()
f.show()


#
# q = np.linspace(0.3, 6, 500)
# I = []
# for cll in [cell, cell_trans, cell_rot]:
#     _vec, _S = diffraction_from_cell(hkl, cll)
#     _vec = [np.sqrt(v[0]**2+v[2]**2) for v in _vec]
#     I.append(powder_pattern(q, _vec, _S, w=0.05))
#
#
# with PdfPages("testdata/simulation_test.pdf") as pdf:
#     f, ax = plt.subplots(3,1, sharex=True, figsize=(11, 8.5))
#     ax[0].plot(q, I[0], c="C0", label="ZnPc")
#     ax[1].plot(q, I[1], c="C1", label="translated")
#     ax[2].plot(q, I[2], c="C2", label="rotated")
#     ax[2].set_xlim(0.3,6)
#     ax[2].set_xlabel("q [A-1]")
#     for a in ax:
#         a.set_yticks([], [])
#         a.legend()
#     f.suptitle("diffraction patterns")
#     f.tight_layout()
#     f.show()
#     pdf.savefig(f)
#
#     f, ax = plt.subplots(2,1, sharex=True, figsize=(11, 8.5))
#     ax[0].plot(q, I[1]-I[0], c="C1", label="translated")
#     ax[1].plot(q, I[2]-I[0], c="C2", label="rotated")
#     ax[1].set_xlim(0.3,6)
#     ax[1].set_xlabel("q [A-1]")
#     for a in ax:
#         # a.set_yticks([], [])
#         a.set_ylabel("delta I [norm]")
#         a.legend()
#     f.suptitle("difference diffraction patterns")
#     f.tight_layout()
#     f.show()
#     pdf.savefig(f)


