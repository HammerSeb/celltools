from cell.tools import symmetry_operator


identity = symmetry_operator(1,1,1,0,0,0, label="id")
inversion = symmetry_operator(-1,-1,-1,0,0,0, label="inversion")

SPACE_GROUP= {
    "1": [identity],
    "2": [identity, inversion]
}