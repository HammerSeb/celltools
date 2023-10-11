from celltools.cell.generate import cell_from_cif, cell_to_crystal


cll = cell_from_cif("testdata/erk.cif")
cll.atoms_to_molecule()

crystal = cell_to_crystal(cll)
