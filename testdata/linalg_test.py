from numpy import deg2rad

from linalg.basis import vector, basis, standard_basis, line, plane
from linalg.transformations import basis_transformation, rotation
from linalg.find import average_line, average_plane


base = basis([0.5,0,0], [0,0.5,0], [0,0,0.5])
axis = vector([1,1,1], base)

rot = rotation(deg2rad(120), axis)

pnt = vector([1,0,0], base)
new_pnt = rot.rotate(pnt)