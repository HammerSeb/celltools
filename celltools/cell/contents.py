from itertools import product, combinations
from copy import copy, deepcopy
import numpy as np
from crystals.atom import Element
from crystals.lattice import Lattice

from . import sort2lists

from celltools.linalg.basis import vector, basis, standard_basis
from celltools.linalg.basis import LinearAlgebraError

# from cell.tools import move

def auto_label_atoms(atms):
    """
    adds labels to a given list of atoms. The label will be the element letter and a number from 1 to the number of
    atoms of that elements. E.g., if there is for carbon atoms in the list the labels will be "C1", "C2", "C3", "C4".

    Parameters
    ----------
    atms: iterable of :class:`atom`
        atoms to be labeled
    """
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
        atom.label = element + str(no_of_elements[_idx])
        no_of_elements[_idx] -= 1

def auto_label_molecules(molcs):
    """
    adds labels to a given list of Molecules. Every molecule label will be numbered up to the number of occurrences of
    that label in the list of molecules.
    ----------
    molcs: list of :class:`molecule'
    """
    list_of_labels = []
    for molc in molcs:
        list_of_labels.append(molc.label)
    labels = []
    for label in list_of_labels:
        if label not in labels:
            labels.append(label)
    no_of_labels = []
    for label in labels:
        no_of_labels.append(list_of_labels.count(label))

    for label, molc in zip(list_of_labels, molcs):
        _idx = labels.index(label)
        molc.label = label + str(no_of_labels[_idx])
        no_of_labels[_idx] -= 1

def chemical_formula(atms):
    """
    Returns the chemical formula of a molecule,
    Parameters
    ----------
    atms: iterable of :class:`atom` or class:`molecule`

    Returns
    -------
    chem_form: str
        chemical formula of the given molecule
    """
    if isinstance(atms, molecule):
        atms = atms.atoms
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

    no_of_elements, elements = sort2lists(no_of_elements, elements)

    chem_form = ''
    for e, n in zip(elements, no_of_elements):
        chem_form += e + str(n)

    return chem_form

def auto_bonds(atm_list, rmin=1, rmax=2):
    """
    generates list of bonds from given atom list. A bond between two atoms is generated if the distances of the two
    atoms lies between rmin and rmax
    Parameters
    ----------
    atm_list: list of :class:`atom`
        atom list from which bonds are created
    rmin: float
        minimum distance between two atoms to create a bond
    rmax: float
        maximum distance between two atoms to create a bond

    Returns
    -------
        list of :class:`bond`
    """
    bonds = []
    atm_pairs = combinations(atm_list, 2)
    for pair in atm_pairs:
        _dist = (pair[0].coords - pair[1].coords).abs_global
        if rmin < _dist <rmax:
            bonds.append(bond(pair[0], pair[1]))
    return bonds

class atom(Element):
    """
    atom class
    """
    def __init__(self, atm, position=None, label=None):
        """

        :param atm:
        :param position:
        """
        super().__init__(atm)
        if not position:
            self._v = None
        else:
            self.coords = position

        if label:
            self._label=label
        else:
            self._label=None

    def __repr__(self):
        if self.label:
            return f"< {self.label} @ {self.coords}>"
        else:
            return f"< {self.element} @ {self.coords}>"

    def __eq__(self, other):
        if isinstance(other, atom):
            if self.element == other.element and self.coords == other.coords:
                return True
            else:
                return False
        else:
            raise TypeError("must be compared to atom instance")

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
            if atm1.coords.basis == atm2.coords.basis:
                self._atm1 = atm1
                self._atm2 = atm2
            else:
                raise LinearAlgebraError("positions of atm1 and atm2 must be in the same basis.")
        else:
            raise ValueError("input needs to be instance of atom.")

    def __repr__(self):
        return f"< Bond: [{self._atm1.label}; {self._atm2.label}] >"

    @property
    def length(self):
        return (self._atm1.coords - self._atm2.coords).abs_global

    @property
    def bond(self):
        return (self._atm1, self._atm2)

class molecule:
    def __init__(self, atms, bonds=None, label=None):
        """
        saves a group of atoms as a molecule. molecular bonds can be added if needed
        Parameters
        ----------
        atms: iterable of atoms in same basis
        bonds (optional): iterable of bonds
        """
        self._atoms = atms
        auto_label_atoms(self._atoms)

        if label:
            self.label = label
        else:
           self.label = chemical_formula(self.atoms)
        self._bonds = []

    def __repr__(self):
        return f"< Molecule {self.label} >"

    def __getitem__(self, item):
        if isinstance(item, str):
            if item in self.atom_labels:
                return self.atoms[self.atom_labels.index(item)]
            else:
                raise IndexError("atom label not in atom labels")
        if isinstance(item, int):
            return self.atoms[item]

    @property
    def atoms(self):
        return self._atoms

    @property
    def bonds(self):
        return self._bonds

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

    def add_bond(self,bnd):
        """

        Parameters
        ----------
        bnd

        Returns
        -------

        """
        if isinstance(bnd, bond):
            self._bonds.append(bnd)
        else:
            self._bonds.append(bond(bnd[1], bnd[1]))

    def auto_bonds(self, rmin=1, rmax=2):
        self._bonds = auto_bonds(self.atoms, rmin=rmin, rmax=rmax)

class lattice(basis):
    def __init__(self, Basis):
        if isinstance(Basis, basis):
            self._basis = Basis.basis
        else:
            super().__init__(self, Basis[0], Basis[1], Basis[2])
        self._offset = np.zeros((3,))

class cell:
    """
    Class representing a unit cell with lattice vectors and a basis comprised of atoms or molecules
    Parameters
    ----------
    latt: :class:`lattice` or :class:`linalg.basis`
        defining the lattice vectors in which the basis coordinates are expressed
    objs: iterable of :class:`atom` or :class:`molecule`
        basis of the unit cell, i.e. the atoms in the unit cell. Coordinates must be expressed in the lattice basis
    """
    def __init__(self, latt, objs):
        self._atoms, self._molecules = [], []
        self.set_lattice(latt)
        self.set_base(objs)
        auto_label_atoms(self._atoms)
        auto_label_molecules(self._molecules)

    @property
    def lattice(self):
        """lattice systems"""
        return self._lattice

    @property
    def base(self):
        """
        cell content, i.e. atoms or molecules
        Returns
        -------
        _atoms: list of :class:`atom`
            list of "uncorrelated" atoms in the unit cell
        _molecules: list of :class:`molecule`
            list of molecules in the unit cell
        """
        return self._atoms, self._molecules

    @property
    def atoms(self):
        """ atoms in the unit cell """
        return self._atoms

    @property
    def molecules(self):
        """ molecules in the unit cell"""
        return self._molecules

    def set_lattice(self, latt):
        if isinstance(latt, basis) or isinstance(latt, lattice):
            self._lattice = latt
        else:
            raise ValueError("latt needs to be instance of basis or lattice")

    def set_base(self, objs):
        """
        sets the crystal base of the cell
        Parameters
        ----------
        objs: iterable of :class:`atom` or :class:`molecule`
            atoms or molecules in the unit cell. Atom positions must be given in lattice coordinates.
        """
        for obj in objs:
            if isinstance(obj, atom):
                self._atoms.append(obj)
            if isinstance(obj, molecule):
                self._molecules.append(obj)

    def add_atom(self, atm):
        """
        adds uncorrelated atom to the unit cell
        Parameters
        ----------
        atm: :class:`atom`
        """
        self._atoms.append(atm)

    def add_molecule(self, molc):
        """
        adds molecule to the unit cell
        Parameters
        ----------
        molc: :class:`molecule`
        """
        self._molecules.append(molc)

    def atoms_to_molecule(self):
        """
        converts atom list into one molecule in list
        """
        self.add_molecule(molecule(self.atoms))
        self._atoms = []


