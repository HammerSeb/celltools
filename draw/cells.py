from draw.draw import GLPoints, GLLines, draw_basis, draw_frame
from cell.atom_data import ELEM_TO_COLOR, ELEM_TO_SIZE
from linalg.basis import basis
# TODO: draw_atom, draw_molecule, draw_lattice, draw_cell, draw_supercell


def _add_atom(GLPts, atm):
    GLPts.add_point(atm.coords, ELEM_TO_SIZE(atm.element), ELEM_TO_COLOR(atm.element))


def draw_cell(w, cell, lw=3, s=100):
    _frame = draw_frame(w, cell.lattice, lw=lw)
    _lattice = draw_basis(w, cell.lattice, lw=lw+1)
    _content = GLPoints()
    for atm in cell.atoms:
        _add_atom(_content, atm)
    w.addItem(_content)
    # for molc in cell.molecules:
    #     draw_molecule(w, molc, s=s)
#
# def draw_supercell(ax, supercell, c="black", lw=1.5, s=100):
#     for trans_vec in supercell._translation_vector:
#         _base = basis(*supercell.lattice.basis, offset=trans_vec.global_coord)
#         draw_frame(ax, _base, c=c, lw=lw)
#     draw_basis(ax, supercell.lattice)
#     for atm in supercell.atoms:
#         draw_atom(ax, atm, s=s)
#     for molc in supercell.molecules:
#         draw_molecule(ax, molc, s=s)