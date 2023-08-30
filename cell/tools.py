from contents import atom, molecule
#TODO: function to generate avergage plane through a list of atoms

def move(obj, vector):
    """
    ONLY FUNCTION SKETCH
    takes a cell object, e.g. atom or molecule and moves it in the direction of vector
    Parameters
    ----------
    obj
    vector

    Returns
    -------
    new_obj
    """
    if isinstance(obj, atom):
        pass
    elif isinstance(obj, molecule):
        for atom in molecule:
            _new_atom = move(atom, vector)
    new_obj = None
    return new_obj
