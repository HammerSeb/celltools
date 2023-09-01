from draw.draw import draw_point, draw_line, draw_vector, draw_basis, draw_frame
from cell.atom_data import ELEM_TO_COLOR
from linalg.basis import basis
# TODO: draw_atom, draw_molecule, draw_lattice, draw_cell, draw_supercell

from crystals import atom_data


def draw_atom(ax, atm, s=100):
    """
    Draws instance of :class:`atom` into given axis
    Parameters
    ----------
    ax: :class:
    atm: :class:`atom`
    s: float (default 100)
        point size for atom depiction
    """
    draw_point(ax, atm.coords, s=s, c=ELEM_TO_COLOR(atm.element))


def draw_bond(ax, bnd):
    pass

def draw_molecule(ax, molc, s=100):
    """
    Draws molecule into given axis
    Parameters
    ----------
    ax
    molc
    s
    """
    for atm in molc.atoms:
        draw_atom(ax, atm, s=s)
    for bnd in molc.bonds:
        draw_bond(ax, bnd)

def draw_cell(ax, cell, c="black", lw=1.5, s=100):
    draw_frame(ax, cell.lattice, c=c, lw=lw)
    draw_basis(ax, cell.lattice)
    for atm in cell.atoms:
        draw_atom(ax, atm, s=s)
    for molc in cell.molecules:
        draw_molecule(ax, molc, s=s)

def draw_supercell(ax, supercell, c="black", lw=1.5, s=100):
    for trans_vec in supercell._translation_vector:
        _base = basis(*supercell.lattice.basis, offset=trans_vec.global_coord)
        draw_frame(ax, _base, c=c, lw=lw)
    draw_basis(ax, supercell.lattice)
    for atm in supercell.atoms:
        draw_atom(ax, atm, s=s)
    for molc in supercell.molecules:
        draw_molecule(ax, molc, s=s)