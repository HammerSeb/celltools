from typing import Optional, List

from pyqtgraph import opengl as gl

from celltools.cell import Atom, Molecule, Cell, SuperCell
from celltools.cell.atom_data import ELEM_TO_COLOR, ELEM_TO_SIZE
from celltools.linalg.basis import Basis
from .draw import GLPoints, GLLines, draw_basis, draw_frame


def _add_atom(GLPts: GLPoints, atom: Atom):
    """
    adds atom to GLPoints object
    Parameters
    ----------
    GLPts: :class:`GlPoints`
    atom: :class:`Atom`
    """
    GLPts.add_point(atom.coords, ELEM_TO_SIZE(atom.element), ELEM_TO_COLOR(atom.element))


def _add_molecule(GLPts: GLPoints, molc: Molecule, GLLns: Optional[GLLines]):
    """
    adds the atoms of a molecule to GLPoints objects. If a GLLines object is given, the bonds of the molecule are added
    as lines.
    Parameters
    ----------
    GLPts: :class:`GLPoints`
    molc: :class:`Molecule`
    GLLns: :class:`GLLines` (Optional)
    """
    for atm in molc.atoms:
        GLPts.add_point(atm.coords, ELEM_TO_SIZE(atm.element), ELEM_TO_COLOR(atm.element))
    if GLLns:
        for bond in molc.bonds:
            if bond.bond[0].atomic_number > bond.bond[1].atomic_number:
                color = ELEM_TO_COLOR(bond.bond[0].element)
            else:
                color = ELEM_TO_COLOR(bond.bond[1].element)
            GLLns.add_line(bond.bond[0].coords, bond.bond[1].coords, c=color)


def draw_cell(w: gl.GLViewWidget, cell: Cell, lw: float = 3) -> List:
    """
    draws a unit cell into view widget
    Parameters
    ----------
    w: pyqtgraph.opengl.GLViewWidget
        view widget, use make_figure function to generate
    cell: :class:`Cell`
        unit cell to draw
    lw: float
        line width for frame indicating unit cell extent

    Returns
    -------
        list of :class:`GLPoints` and :class:`GLLines`
            containing all atoms and molecular atoms and bonds
    """
    _frame = draw_frame(w, cell.lattice, lw=lw)
    _lattice = draw_basis(w, cell.lattice, lw=lw + 4)
    _content = []
    if cell.atoms:
        _content.append(GLPoints())
        for atm in cell.atoms:
            _add_atom(_content[0], atm)
    if cell.molecules:
        for molc in cell.molecules:
            _content.append(GLPoints())
            if molc.bonds:
                _content.append(GLLines())
                _add_molecule(_content[-2], molc, GLLns=_content[-1])
            else:
                _add_molecule(_content[-1], molc)
    for cnt in _content:
        w.addItem(cnt)
    return _content


def draw_supercell(w: gl.GLViewWidget, supercell: SuperCell, lw: float = 3) -> List:
    """
    draws a super cell into view widget
    Parameters
    ----------
    w: pyqtgraph.opengl.GLViewWidget
        view widget, use make_figure function to generate
    supercell: :class:`SuperCell`
        super cell to draw
    lw: float
        line width for frame indicating unit cell extent

    Returns
    -------
        list of :class:`GLPoints` and :class:`GLLines`
            containing all atoms and molecular atoms and bonds
    """
    for trans_vec in supercell._translation_vector:
        _base = Basis(*supercell.lattice.basis, offset=trans_vec.global_coord)
        __ = draw_frame(w, _base, lw=lw)
    __ = draw_basis(w, supercell.lattice, lw=lw + 4)
    _content = []
    if supercell.atoms:
        _content.append(GLPoints())
        for atm in supercell.atoms:
            _add_atom(_content[0], atm)
    if supercell.molecules:
        for molc in supercell.molecules:
            _content.append(GLPoints())
            if molc.bonds:
                _content.append(GLLines())
                _add_molecule(_content[-2], molc, GLLns=_content[-1])
            else:
                _add_molecule(_content[-1], molc)
    for cnt in _content:
        w.addItem(cnt)
    return _content
