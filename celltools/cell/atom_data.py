import typing
import warnings


class NotImplementedWarning(Warning):
    """ Warning to signalize that a particular element has not been given a color and size yet for drawing """
    pass


property_dict = {
    # 1
    "H": [[0.7, 0.7, 0.7, 1], 0.3],
    "He": [[1, 1, 1, 1], 0.3],
    # 2
    "Li": [[1, 1, 1, 1], 0.3],
    "Be": [[1, 1, 1, 1], 0.3],
    "B": [[1, 1, 1, 1], 0.3],
    "C": [[0.4, 0.4, 0.4, 1], 0.7],
    "N": [[0.43, 0.145, 0.89, 1], 0.7],
    "O": [[0.278, 0.518, 0.851, 1], 0.7],
    "F": [[1, 1, 1, 1], 0.3],
    "Ne": [[1, 1, 1, 1], 0.3],
    # 3
    "Na": [[1, 1, 1, 1], 0.3],
    "Mg": [[1, 1, 1, 1], 0.3],
    "Al": [[1, 1, 1, 1], 0.3],
    "Si": [[1, 1, 1, 1], 0.3],
    "P": [[1, 1, 1, 1], 0.3],
    "S": [[1, 1, 1, 1], 0.3],
    "Cl": [[1, 1, 1, 1], 0.3],
    "Ar": [[1, 1, 1, 1], 0.3],
    # 4
    "K": [[1, 1, 1, 1], 0.3],
    "Ca": [[1, 1, 1, 1], 0.3],
    "Sc": [[1, 1, 1, 1], 0.3],
    "Ti": [[1, 1, 1, 1], 0.3],
    "V": [[1, 1, 1, 1], 0.3],
    "Cr": [[1, 1, 1, 1], 0.3],
    "Mn": [[1, 1, 1, 1], 0.3],
    "Fe": [[1, 1, 1, 1], 0.3],
    "Co": [[1, 1, 1, 1], 0.3],
    "Ni": [[1, 1, 1, 1], 0.3],
    "Cu": [[1, 1, 1, 1], 0.3],
    "Zn": [[0.7, 0, 0.65, 1], 1],
    "Ga": [[1, 1, 1, 1], 0.3],
    "Ge": [[1, 1, 1, 1], 0.3],
    "As": [[1, 1, 1, 1], 0.3],
    "Se": [[1, 1, 1, 1], 0.3],
    "Br": [[1, 1, 1, 1], 0.3],
    "Kr": [[1, 1, 1, 1], 0.3],
    # 5
    "Rb": [[1, 1, 1, 1], 0.3],
    "Sr": [[1, 1, 1, 1], 0.3],
    "Y": [[1, 1, 1, 1], 0.3],
    "Zr": [[1, 1, 1, 1], 0.3],
    "Nb": [[1, 1, 1, 1], 0.3],
    "Mo": [[1, 1, 1, 1], 0.3],
    "Tc": [[1, 1, 1, 1], 0.3],
    "Ru": [[1, 1, 1, 1], 0.3],
    "Rh": [[1, 1, 1, 1], 0.3],
    "Pd": [[1, 1, 1, 1], 0.3],
    "Ag": [[1, 1, 1, 1], 0.3],
    "Cd": [[1, 1, 1, 1], 0.3],
    "In": [[1, 1, 1, 1], 0.3],
    "Sn": [[1, 1, 1, 1], 0.3],
    "Sb": [[1, 1, 1, 1], 0.3],
    "Te": [[1, 1, 1, 1], 0.3],
    "I": [[1, 1, 1, 1], 0.3],
    "Xe": [[1, 1, 1, 1], 0.3],
    # 6
    "Cs": [[1, 1, 1, 1], 0.3],
    "Ba": [[1, 1, 1, 1], 0.3],
    "La": [[1, 1, 1, 1], 0.3],
    "Ce": [[1, 1, 1, 1], 0.3],
    "Pr": [[1, 1, 1, 1], 0.3],
    "Nd": [[1, 1, 1, 1], 0.3],
    "Pm": [[1, 1, 1, 1], 0.3],
    "Sm": [[1, 1, 1, 1], 0.3],
    "Eu": [[1, 1, 1, 1], 0.3],
    "Gd": [[1, 1, 1, 1], 0.3],
    "Tb": [[1, 1, 1, 1], 0.3],
    "Dy": [[1, 1, 1, 1], 0.3],
    "Ho": [[1, 1, 1, 1], 0.3],
    "Er": [[1, 1, 1, 1], 0.3],
    "Tm": [[1, 1, 1, 1], 0.3],
    "Yb": [[1, 1, 1, 1], 0.3],
    "Lu": [[1, 1, 1, 1], 0.3],
    "Hf": [[1, 1, 1, 1], 0.3],
    "Ta": [[1, 1, 1, 1], 0.3],
    "W": [[1, 1, 1, 1], 0.3],
    "Re": [[1, 1, 1, 1], 0.3],
    "Os": [[1, 1, 1, 1], 0.3],
    "Ir": [[1, 1, 1, 1], 0.3],
    "Pt": [[1, 1, 1, 1], 0.3],
    "Au": [[1, 1, 1, 1], 0.3],
    "Hg": [[1, 1, 1, 1], 0.3],
    "Tl": [[1, 1, 1, 1], 0.3],
    "Pb": [[1, 1, 1, 1], 0.3],
    "Bi": [[1, 1, 1, 1], 0.3],
    "Po": [[1, 1, 1, 1], 0.3],
    "At": [[1, 1, 1, 1], 0.3],
    "Rn": [[1, 1, 1, 1], 0.3],
    # 7
    "Fr": [[1, 1, 1, 1], 0.3],
    "Ra": [[1, 1, 1, 1], 0.3],
    "Ac": [[1, 1, 1, 1], 0.3],
    "Th": [[1, 1, 1, 1], 0.3],
    "Pa": [[1, 1, 1, 1], 0.3],
    "U": [[1, 1, 1, 1], 0.3],
    "Np": [[1, 1, 1, 1], 0.3],
    "Pu": [[1, 1, 1, 1], 0.3],
    "Am": [[1, 1, 1, 1], 0.3],
    "Cm": [[1, 1, 1, 1], 0.3],
    "Bk": [[1, 1, 1, 1], 0.3],
    "Cf": [[1, 1, 1, 1], 0.3],
    "Es": [[1, 1, 1, 1], 0.3],
    "Fm": [[1, 1, 1, 1], 0.3],
    "Md": [[1, 1, 1, 1], 0.3],
    "No": [[1, 1, 1, 1], 0.3],
    "Lr": [[1, 1, 1, 1], 0.3],
    "Rf": [[1, 1, 1, 1], 0.3],
    "Db": [[1, 1, 1, 1], 0.3],
    "Sg": [[1, 1, 1, 1], 0.3],
    "Bh": [[1, 1, 1, 1], 0.3],
    "Hs": [[1, 1, 1, 1], 0.3],
    "Mt": [[1, 1, 1, 1], 0.3],
    "Ds": [[1, 1, 1, 1], 0.3],
    "Rg": [[1, 1, 1, 1], 0.3],
    "Cn": [[1, 1, 1, 1], 0.3],
    "Nh": [[1, 1, 1, 1], 0.3],
    "Fl": [[1, 1, 1, 1], 0.3],
    "Mc": [[1, 1, 1, 1], 0.3],
    "Lv": [[1, 1, 1, 1], 0.3],
    "Ts": [[1, 1, 1, 1], 0.3],
    "Og": [[1, 1, 1, 1], 0.3]
}


def ELEM_TO_COLOR(element: str) -> list:
    """
    looks up element in dictionary and returns color as rgba-list. Raises a warning if element is not implemented yet.
    Parameters
    ----------
    element: str
        string specifing the element

    Returns
    -------
        list
    """
    _implemented = ['H', 'C', 'N', 'O', 'Zn']
    if element in _implemented:
        return property_dict[element][0]
    else:
        warnings.warn("element is not implemented yet, return will be white", NotImplementedWarning)
        return property_dict[element][0]


def ELEM_TO_SIZE(element: str) -> float:
    """
    looks up element in dictionary and returns size as float. Raises a warning if element is not implemented yet.
    Parameters
    ----------
    element: str
        string specifing the element

    Returns
    -------
        float
    """
    _implemented = ['H', 'C', 'N', 'O', 'Zn']
    if element in _implemented:
        return property_dict[element][1]
    else:
        warnings.warn("element is not implemented yet, return will be 0.3", NotImplementedWarning)
        return property_dict[element][1]
