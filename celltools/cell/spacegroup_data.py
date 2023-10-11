from celltools.cell.tools import SymmetryOperator


identity = SymmetryOperator(1,1,1,0,0,0, label="id")
inversion = SymmetryOperator(-1,-1,-1,0,0,0, label="inversion", id=[[0,0,0]])

SPACE_GROUP= {
    "1": [identity],
    "2": [identity, inversion]
}