from celltools.linalg.basis import basis, standard_basis, vector
from celltools.cell.contents import atom, molecule, auto_label_atoms, chemical_formula, lattice, cell
from celltools.cell.generate import lattice_from_cell_parameters, cell_from_crystal, cell_from_cif, cell_to_crystal
from celltools.cell.tools import move,  super_cell

from crystals import Crystal
from celltools.draw.draw import make_figure
from celltools.draw.cells import draw_supercell


from matplotlib import pyplot as plt


cll = cell_from_cif("erk.cif")
cll.atoms_to_molecule()

crystal = cell_to_crystal(cll)