from copy import deepcopy
from itertools import product
import numpy as np
from simulation.ediff import powder_pattern, diffraction_from_cell
from celltools.cell.generate import cell_from_cif, cell_to_crystal
from celltools.cell.tools import rotate, move, basis_transformation
from celltools.linalg.basis import vector, line, standard_basis
from celltools.linalg.find import average_plane

from celltools.draw.draw import make_figure
from celltools.draw.cells import draw_cell

from matplotlib import pyplot as plt
import matplotlib

import pyqtgraph as pg


matplotlib.use("Qt5Agg")

cell = cell_from_cif("testdata/erk.cif")
cell.atoms_to_molecule()
cell.molecules[0].auto_bonds()

#move atom along a-axis 0.4 Angstrom
cell_trans = deepcopy(cell)
trans_vec = vector([1,0,0], cell_trans.lattice)
trans_vec *= 0.4*(1/trans_vec.abs_global)
move(cell_trans.molecules[0], trans_vec)


#rotate 45 deg around normal through ZnPc atom
cell_rot = deepcopy(cell)
molc_coords = [atm.coords for atm in cell_rot.molecules[0].atoms]
ZnPc_plane = average_plane(molc_coords)
std_to_lattice = basis_transformation(standard_basis, cell_rot.lattice)
rotation_axis = line(vector([0,0,0], cell_rot.lattice), std_to_lattice.transform(ZnPc_plane.normal))
rotate(cell_rot.molecules[0], rotation_axis, 45)

# w = make_figure(grid=False)
# draw_cell(w, cell_rot)
# pg.exec()


# hkl = [(0,k+1,0) for k in range(4)]

hkl = list(product(
    [i-5 for i in range(11)],
    [0],
    [i-5 for i in range(11)]
    ))
hkl = hkl[int((len(hkl)-1)/2):]

q = np.linspace(0.3, 6, 500)
I = []
for cll in [cell, cell_trans, cell_rot]:
    _vec, _S = diffraction_from_cell(hkl, cll)
    _vec = [np.sqrt(v[0]**2+v[2]**2) for v in _vec]
    I.append(powder_pattern(q, _vec, _S, w=0.05))


f, ax = plt.subplots(3,1, sharex=True)
ax[0].plot(q, I[0])
ax[1].plot(q, I[1])
ax[2].plot(q, I[2])
f.tight_layout()
f.show()
