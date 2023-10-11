from numpy import deg2rad

from celltools.linalg.basis import Vector, Basis
from celltools.linalg.transformations import Rotation


base = Basis([0.5,0,0], [0,0.5,0], [0,0,0.5])
axis = Vector([1,1,1], base)

rot = Rotation(deg2rad(120), axis)

pnt = Vector([1,0,0], base)
new_pnt = rot.rotate(pnt)
