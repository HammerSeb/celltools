import typing
from typing import List, Tuple, Union, Optional
from itertools import combinations

import numpy as np
from crystals.atom import Element, ElementLike

from celltools.linalg import Vector, Basis
from celltools.linalg.basis import LinearAlgebraError
from . import sort2lists


def auto_label_atoms(atms: List['Atom']) -> None:
    """
    adds labels to a given list of atoms. The label will be the element letter and a number from 1 to the number of
    atoms of that elements. E.g., if there is for carbon atoms in the list the labels will be "C1", "C2", "C3", "C4".

    Parameters
    ----------
    atms: iterable of :class:`Atom`
        atoms to be labeled
    """
    list_of_elements = []
    for atm in atms:
        list_of_elements.append(atm.element)
    elements = []
    for element in list_of_elements:
        if element not in elements:
            elements.append(element)
    no_of_elements = []
    for element in elements:
        no_of_elements.append(list_of_elements.count(element))

    for element, atm in zip(list_of_elements, atms):
        _idx = elements.index(element)
        atm.label = element + str(no_of_elements[_idx])
        no_of_elements[_idx] -= 1


def auto_label_molecules(molcs: List['Molecule']) -> None:
    """
    adds labels to a given list of Molecules. Every molecule label will be numbered up to the number of occurrences of
    that label in the list of molecules.
    ----------
    molcs: list of :class:`Molecule`
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


def chemical_formula(atms: Union[List['Atom'], 'Molecule']) -> str:
    """
    Returns the chemical formula of a molecule,
    Parameters
    ----------
    atms: iterable of :class:`Atom` or class:`Molecule`

    Returns
    -------
    chem_form: str
        chemical formula of the given molecule
    """
    if isinstance(atms, Molecule):
        atms = atms.atoms
    list_of_elements = []
    for atm in atms:
        list_of_elements.append(atm.element)
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


def auto_bonds(atm_list: List['Atom'], rmin: float = 1, rmax: float = 2) -> List['Bond']:
    """
    generates list of bonds from given atom list. A bond between two atoms is generated if the distances of the two
    atoms lies between rmin and rmax
    Parameters
    ----------
    atm_list: list of :class:`Atom`
        atom list from which bonds are created
    rmin: float
        minimum distance between two atoms to create a bond
    rmax: float
        maximum distance between two atoms to create a bond

    Returns
    -------
        list of :class:`Bond`
    """
    bonds = []
    atm_pairs = combinations(atm_list, 2)
    for pair in atm_pairs:
        _dist = (pair[0].coords - pair[1].coords).abs_global
        if rmin < _dist < rmax:
            bonds.append(Bond(pair[0], pair[1]))
    return bonds


class Atom(Element):
    """
    Class representing a particular atom. A position and a label can be specified.
    Additionally, it inherits all methods from class:`crystals.atom.Element`
    Parameters
    ----------
    atm: ElementLike
        str or int specifing the element from the periodic table

    position: :class:`Vector` (optional)
        position of the Atom

    label: str (optional)
        label of the atom
    """

    def __init__(self, atm: ElementLike, position: Optional[Vector] = None, label: Optional[str] = None):
        super().__init__(atm)
        if not position:
            self._v = None
        else:
            self.coords = position

        if label:
            self._label = label
        else:
            self._label = None

    def __repr__(self) -> str:
        if self.label:
            return f"< {self.label} @ {self.coords}>"
        else:
            return f"< {self.element} @ {self.coords}>"

    def __eq__(self, other: 'Atom') -> bool:
        if isinstance(other, Atom):
            if self.element == other.element and self.coords == other.coords:
                return True
            else:
                return False
        else:
            raise TypeError("must be compared to atom instance")

    @property
    def coords(self) -> Vector:
        """ returns atom position """
        if self._v:
            return self._v

    @property
    def label(self) -> str:
        """ returns atom label """
        if self._label:
            return self._label

    @coords.setter
    def coords(self, position: Union[Tuple[np.ndarray, Basis], Vector]):
        """
        adds coordinates to the atom as a vector

        Parameters
        ----------
        position: :class:`Vector` or list like [position as ndarray (3,), basis instance]
        """
        if isinstance(position, Vector):
            self._v = position
        else:
            self._v = Vector(position[0], position[1])

    @label.setter
    def label(self, label: str):
        """
        adds label to atom
        Parameters
        ----------
        label: str
            atom label
        """
        self._label = label


class Bond:
    """
    class representing a bond between two atoms.Position of atoms must be given in same basis.

    Parameters
    ----------
    atm1: :class:`Atom`
    atm2: :class:`Atom`
    """

    def __init__(self, atm1: Atom, atm2: Atom):
        if isinstance(atm1, Atom) and isinstance(atm2, Atom):
            if atm1.coords.basis == atm2.coords.basis:
                self._atm1 = atm1
                self._atm2 = atm2
            else:
                raise LinearAlgebraError("positions of atm1 and atm2 must be in the same basis.")
        else:
            raise ValueError("input needs to be instance of atom.")

    def __repr__(self) -> str:
        return f"< Bond: [{self._atm1.label}; {self._atm2.label}] >"

    @property
    def length(self) -> float:
        """
        returns bond length

        Returns
        -------
            float
        """
        return (self._atm1.coords - self._atm2.coords).abs_global

    @property
    def bond(self) -> Tuple[Atom, Atom]:
        """
        returns tuple of atoms of the bond
        Returns
        -------
            tuple
                2-tuple of atoms
        """
        return (self._atm1, self._atm2)


class Molecule:
    """
    class representing a molecule as a container for a list of atoms. All atom positions must be given in same basis
    Parameters
    ----------
    atms: iterable of :class:`atom`
        atoms comprising the molecue
    bonds: iterbale of :class:`bonds` (optional)
        bonds between the molecule's atoms
    label: str (optional)
        label of molecule

    """

    def __init__(self, atms: List[Atom], bonds: Optional[List[Bond]] = None, label: Optional[str] = None):
        self._atoms = atms
        auto_label_atoms(self._atoms)

        if label:
            self.label = label
        else:
            self.label = chemical_formula(self.atoms)

        if bonds:
            self._bonds = bonds
        else:
            self._bonds = []

    def __repr__(self) -> str:
        return f"< Molecule {self.label} >"

    def __getitem__(self, item: Union[str, int]) -> Atom:
        """
        returns atom from atom label or index
        Parameters
        ----------
        item: int or str
            either index of atom list or atom label
        Returns
        -------
            :class:`atom`
        """
        if isinstance(item, str):
            if item in self.atom_labels:
                return self.atoms[self.atom_labels.index(item)]
            else:
                raise IndexError("atom label not in atom labels")
        if isinstance(item, int):
            return self.atoms[item]

    @property
    def atoms(self) -> List[Atom]:
        """ returns atom list """
        return self._atoms

    @property
    def bonds(self) -> List[Bond]:
        """ returns bond list """
        return self._bonds

    @property
    def atom_labels(self) -> List[str]:
        """ returns list of atom labels """
        return [atm.label for atm in self.atoms]

    @property
    def label(self) -> str:
        """ returns label of molecule """
        return self._label

    @atom_labels.setter
    def atom_labels(self, labels: List[str]):
        """
        sets atom labels for all atoms according to list of strings
        Parameters
        ----------
        labels: list of str
            labels of each atom in same order as atom list

        Raises
        -------
            ValueError
                if list of labels is not the same length as list of atoms
        """
        if len(labels) != len(self.atoms):
            raise ValueError("lengt of list of labels needs to be the same as number of atoms")
        else:
            for label, atom in zip(labels, self._atoms):
                atom.label(label)

    @label.setter
    def label(self, label: str):
        """ set label of molecule """
        self._label = label

    def add_bond(self, bnd: Union[Bond, Tuple[Atom, Atom]]):
        """
        adds a bond to the molecule
        Parameters
        ----------
        bnd: :class:`bond` or list [:class:`atom`, :class:`atom`]
            bond to be added. (make sure atoms of bond are atoms of the molecule, no additional check performed!)
        """
        if isinstance(bnd, Bond):
            self._bonds.append(bnd)
        else:
            self._bonds.append(Bond(bnd[1], bnd[1]))

    def auto_bonds(self, rmin: float = 1, rmax: float = 2):
        """
        uses function auto_bonds on molecule's atoms. rmin, rmax same as in auto_bonds.
        """
        self._bonds = auto_bonds(self.atoms, rmin=rmin, rmax=rmax)


class Lattice(Basis):
    """
    crystal lattice, inherits from basis. After initalization the offset is always set to lattice.offset = [0,0,0].
    Parameters
    ----------
    Basis: :class:basis or list of 3 np.ndarray size (3,)
        define the lattice Vectors either as a basis or as three basis Vectors
    """

    def __init__(self, basis: Union[Basis, Tuple[np.ndarray, np.ndarray, np.ndarray]]):

        if isinstance(basis, Basis):
            self._basis = basis.basis
        else:
            super().__init__(self, basis[0], basis[1], basis[2])

        self._offset = np.zeros((3,))


class Cell:
    """
    Class representing a unit cell with lattice Vectors and a basis comprised of atoms or molecules
    Parameters
    ----------
    latt: :class:`Lattice` or :class:`linalg.Basis`
        defining the lattice Vectors in which the basis coordinates are expressed
    objs: iterable of :class:`Atom` or :class:`Molecule`
        basis of the unit cell, i.e. the atoms in the unit cell. Coordinates must be expressed in the lattice basis
    """

    def __init__(self, latt: Lattice, objs: Union[List[Atom], List[Molecule]]):
        self._atoms, self._molecules = [], []
        self.set_lattice(latt)
        self.set_base(objs)
        auto_label_atoms(self._atoms)
        auto_label_molecules(self._molecules)

    @property
    def lattice(self) -> Lattice:
        """ lattice systems """
        return self._lattice

    @property
    def base(self) -> Tuple[List[Atom], List[Molecule]]:
        """
        cell content, i.e. atoms or molecules
        Returns
        -------
        _atoms: list of :class:`Atom`
            list of "uncorrelated" atoms in the unit cell
        _molecules: list of :class:`Molecule`
            list of molecules in the unit cell
        """
        return self._atoms, self._molecules

    @property
    def atoms(self) -> List[Atom]:
        """ atoms in the unit cell """
        return self._atoms

    @property
    def molecules(self) -> List[Molecule]:
        """ molecules in the unit cell"""
        return self._molecules

    def set_lattice(self, latt: Union[Lattice, Basis]):
        """
        set lattice of unit cell
        Parameters
        ----------
        latt: :class:`Lattice` or :class:`Basis`

        Raises
        -------
            Value Error
        """
        if isinstance(latt, Basis) or isinstance(latt, Lattice):
            self._lattice = latt
        else:
            raise ValueError("latt needs to be instance of basis or lattice")

    def set_base(self, objs: Union[List[Atom], List[Molecule]]):
        """
        sets the crystal base of the cell
        Parameters
        ----------
        objs: iterable of :class:`Atom` or :class:`Molecule`
            atoms or molecules in the unit cell. Atom positions must be given in lattice coordinates.
        """
        for obj in objs:
            if isinstance(obj, Atom):
                self._atoms.append(obj)
            if isinstance(obj, Molecule):
                self._molecules.append(obj)

    def add_atom(self, atm: Atom):
        """
        adds uncorrelated atom to the unit cell
        Parameters
        ----------
        atm: :class:`atom`
        """
        self._atoms.append(atm)

    def add_molecule(self, molc: Molecule):
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
        self.add_molecule(Molecule(self.atoms))
        self._atoms = []
