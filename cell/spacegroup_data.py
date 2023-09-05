from tools import symmetry_operator


identity = symmetry_operator(1,1,1,0,0,0)
inversion = symmetry_operator(-1,-1,-1,0,0,0)

SPACE_GROUP= {
    "1" : [identity],
    "2" : [identity, inversion]
}