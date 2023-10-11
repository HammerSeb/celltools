from numpy import deg2rad

from linalg.basis import vector, basis
from linalg.transformations import rotation


base = basis([0.5,0,0], [0,0.5,0], [0,0,0.5])
axis = vector([1,1,1], base)

rot = rotation(deg2rad(120), axis)

pnt = vector([1,0,0], base)
new_pnt = rot.rotate(pnt)
