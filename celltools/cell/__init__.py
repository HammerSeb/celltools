import typing
from typing import Tuple

from .contents import Atom, Molecule, Lattice, Cell
from .tools import move, rotate, SuperCell


def sort2lists(lst1: list, lst2: list) -> Tuple[list, list]:
    """
    sorts two unsorted lists according to first given list
    Parameters
    ----------
    lst1: list of sortables
    lst2: list

    Returns
    -------
    l1, l2: sorted lists

    """
    return [l1 for l1, _ in sorted(zip(lst1, lst2))], [l2 for _, l2 in sorted(zip(lst1, lst2))]
