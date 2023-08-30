import numpy as np
from crystals.atom import Element
from crystals.lattice import Lattice

from linalg.basis import vector, basis, standard_basis
from linalg.basis import LinearAlgebraError

from cell.tools import move
# TODO: unit_cell, super_cell

def auto_label_atoms(atms):
    """
    adds labels to a given list of atoms
    Parameters
    ----------
    atms: iterable of instances of atom

    Returns
    -------
    list of labeled atoms

    """
    labeled_atoms = []
    list_of_elements = []
    for atom in atms:
        list_of_elements.append(atom.element)
    elements = []
    for element in list_of_elements:
        if element not in elements:
            elements.append(element)
    no_of_elements = []
    for element in elements:
        no_of_elements.append(list_of_elements.count(element))

    for element, atom in zip(list_of_elements, atms):
        _idx = elements.index(element)
        atom.label(element + str(no_of_elements[_idx]))
        no_of_elements[_idx] -= 1
    return labeled_atoms

class atom(Element):
    """
    atom class
    """
    def __init__(self, atm, position=None, label=None):
        """

        :param atm:
        :param position:
        """
        super().__init__(self, atm)
        if not position:
            self._v = None
        else:
            self.add_coords(position)

        if label:
            self._label=label
        else:
            self._label=None
    @property
    def coords(self):
        if self._v:
            return self._v

    @property
    def label(self):
        if self._label:
            return self._label

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

    @label.setter
    def label(self, label):
        """

        Parameters
        ----------
        label

        Returns
        -------

        """
        self._label = label

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
    def __init__(self, atms, bonds=None, label=None):
        """

        Parameters
        ----------
        atms: iterable of atoms in same basis
        bonds (optional): iterable of bonds
        """
        self._atoms = auto_label_atoms(atms)
        self._label = self.label(label)

    @property
    def atoms(self):
        return self._atoms

    @property
    def atom_labels(self):
        return [atm.label for atm in self.atoms]

    @property
    def label(self):
        return self._label

    @atom_labels.setter
    def atom_labels(self, labels):
        if len(labels) != len(self.atoms):
            raise ValueError("lengt of list of labels needs to be the same as number of atoms")
        else:
            for label, atom in zip(labels, self._atoms):
                atom.label(label)

    @label.setter
    def label(self, label):
        self._label = label

    def __get_item__(self, item):
        if isinstance(item, str):
            if item in self.atom_labels:
                return self.atoms[self.atom_labels.index(item)]
            else:
                raise IndexError("atom label not in atom labels")
        if isinstance(item, int):
            return self.atoms[item]

class lattice(basis):
    def __init__(self, Basis):
        if isinstance(Basis, basis):
            self._basis = Basis.basis
        else:
            super().__init__(self, Basis[0], Basis[1], Basis[2])
        self._offset = None

class cell:
    def __init__(self, latt, objs):
        """

        Parameters
        ----------
        latt: instance of lattice or linalg.basis
        objs: iterable with atoms or molecules in coordinates of lattice
        """

        self._lattice = self.lattice(latt)
        self.base(objs)

    @property
    def lattice(self):
        return self._lattice

    @property
    def base(self):
        return self._atoms, self._molecules

    @property
    def atoms(self):
        return self._atoms

    @property
    def molecules(self):
        return self._molecules

    @lattice.setter
    def lattice(self, latt):
        if isinstance(latt, basis) or isinstance(latt, lattice):
            self._lattice = latt
        else:
            raise ValueError("latt needs to be instance of basis or lattice")

    @base.setter
    def base(self, objs):
        self._atoms, self._molecules = [], []
        for obj in objs:
            if isinstance(obj, atom):
                self._atoms.append(obj)
            if isinstance(obj, molecule):
                self._molecules.append(obj)
                for atm in obj:
                    self._atoms.append(atm)

        self._atoms = auto_label_atoms(self._atoms)

    def add_atom(self, atm):
        self._atoms.append(atm)

    def add_molecule(self, molc):
        self._molecules.append(molc)

class super_cell(cell):
    def __init__(self, unit_cell, size):
        self.lattice(unit_cell.lattice)

            for atm in unit_cell.atoms:
                pass
