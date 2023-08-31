import cell.contents as cc
from linalg.basis import basis, vector
from linalg.transformations import basis_transformation
#TODO: function to generate avergage plane through a list of atoms

def move(obj, vec):
    """
    ONLY FUNCTION SKETCH
    takes a cell object, e.g. atom or molecule and moves it in the direction of vector
    Parameters
    ----------
    obj: :class:`atom` or :class:`molecule`
    vec: :class:`vector`
    """
    if isinstance(obj, cc.atom):
        to_atm_basis = basis_transformation(obj.coords.basis, vec.basis)
        obj.coords = obj.coords + to_atm_basis.inv_transform(vec)
    elif isinstance(obj, cc.molecule):
        to_mol_basis = basis_transformation(obj.atoms[0].coords.basis, vec.basis)
        for atm in obj.atoms:
            atm.coords = atm.coords + to_mol_basis.inv_transform(vec)


