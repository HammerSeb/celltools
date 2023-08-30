import numpy as np
from crystals.atom import Element
from crystals.lattice import Lattice

from linalg.basis import vector, basis, standard_basis
from linalg.basis import LinearAlgebraError
# TODO: atom (inherit from crystals), bond, molecule, unit_cell
# TODO: function to generate unit_cell from Crystals object

class atom(Element):
    """
    atom class
    """
    def __init__(self, atm, position=None):
        """

        :param atm:
        :param position:
        """
        super().__init__(self, atm)
        if not position:
            self._v = None
        else:
            self.add_coords(position)

    @property
    def coords(self):
        return self._v

    @coords.setter
    def coords(self, position):
        """
        add coordinates to the atom as a vector
        Parameters
        ----------
        position: either instance of vector or list like [position as ndarray (3,), basis instance]

        """

        if isinstance(position, vector):
            self._v = position
        else:
            self._v = vector(position[0], position[1])

class bond:
    def __init__(self, atm1, atm2):
        """
        bond between atm1 and atm2
        Parameters
        ----------
        atm1: atom instance with position in same basis
        atm2: atom instance with position in same basis
        """
        if isinstance(atm1, atom) and isinstance(atm2, atom):
            if atm1.coords.basis.basis == atm2.coords.basis.basis:
                self._atm1 = atm1
                self._atm2 = atm2
            else:
                raise LinearAlgebraError("positions of atm1 and atm2 must be in the same basis.")
        else:
            raise ValueError("input needs to be instance of atom.")

    @property
    def length(self):
        return (self._atm1.coords + self._atm2.coords).abs

class molecule:
    """
    TODO: needs iterator
    TODO:
    """
    def __init__(self):
        pass