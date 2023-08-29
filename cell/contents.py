import numpy as np
from crystals.atom import Element
from crystals.lattice import Lattice

from linalg.basis import vector, basis, standard_basis
# TODO: atom (inherit from crystals), bond, molecule, unit_cell
# TODO: function to generate unit_cell from Crystals object

class atom(Element):
    def __init__(self, atm, position=None):
        """
        :atm: str or int, define type of atom
        """
        super().__init__(self, atm)

        if not position:
            self.v = None
        else:
            self.add_coords(position)


    def add_coords(self, position):
        """
        :position: either instance of vector or list like [position as ndarray (3,), basis instance]
        """
        if isinstance(position, vector):
            self.v = position

        else:
            self.v = vector(position[0], position[1])

