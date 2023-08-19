from basis import basis, vector

v1 = [1, 0, 0]
v2 = [0, 1, -1]
v3 = [1, 0, 1]

b = basis(v1,v2,v3)

a = vector([1,1,0], b)
c = vector([0,1,1], b)

