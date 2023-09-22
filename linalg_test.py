from linalg.basis import vector, basis, standard_basis, line, plane
from linalg.transformations import basis_transformation
from linalg.find import average_line

from cell.generate import cell_from_cif


cll = cell_from_cif("testdata/erk.cif")
cll.atoms_to_molecule()

molc_coords = [atm.coords for atm in cll.molecules[0].atoms]
ln = average_line(molc_coords)
