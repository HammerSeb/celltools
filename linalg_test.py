from linalg.basis import vector, basis, standard_basis
from linalg.transformations import basis_transformation

v1 = [1, 0, 0]
v2 = [0, 1, -1]
v3 = [1, 0, 1]

b = basis(v1,v2,v3)

trans = basis_transformation(b, standard_basis)

a = vector([1,1,0], b)
c = vector([0,1,1], b)

a2 = trans.transform(a)
x = vector([1,1,-1], standard_basis)
a3 = trans.inv_transform(x)

print(a.global_coord)
print(a2.vector)
print(a3.vector, a3.global_coord)