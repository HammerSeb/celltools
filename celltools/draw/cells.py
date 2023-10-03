from .draw import GLPoints, GLLines, draw_basis, draw_frame
from celltools.cell.atom_data import ELEM_TO_COLOR, ELEM_TO_SIZE
from celltools.linalg.basis import basis
# TODO: draw_atom, draw_molecule, draw_lattice, draw_cell, draw_supercell


def _add_atom(GLPts,atm):
    GLPts.add_point(atm.coords, ELEM_TO_SIZE(atm.element), ELEM_TO_COLOR(atm.element))

def _add_molecule(GLPts, molc, GLLns=None):
    for atm in molc.atoms:
        GLPts.add_point(atm.coords, ELEM_TO_SIZE(atm.element), ELEM_TO_COLOR(atm.element))
    if GLLns:
        for bond in molc.bonds:
            if bond.bond[0].atomic_number > bond.bond[1].atomic_number:
                color = ELEM_TO_COLOR(bond.bond[0].element)
            else:
                color = ELEM_TO_COLOR(bond.bond[1].element)
            GLLns.add_line(bond.bond[0].coords, bond.bond[1].coords, c=color)



def draw_cell(w, cell, lw=3):
    _frame = draw_frame(w, cell.lattice, lw=lw)
    _lattice = draw_basis(w, cell.lattice, lw=lw+4)
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

def draw_supercell(w, supercell, lw=3):
    for trans_vec in supercell._translation_vector:
            _base = basis(*supercell.lattice.basis, offset=trans_vec.global_coord)
            __ = draw_frame(w, _base, lw=lw)
    __ = draw_basis(w, supercell.lattice, lw=lw+4)
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