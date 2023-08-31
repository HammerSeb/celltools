from linalg.basis import basis, standard_basis, vector
from cell.contents import atom, molecule, auto_label_atoms, chemical_formula
from cell.tools import move

a = atom('C', vector([0,0,0]))
atms = [atom(el, pos) for el, pos in zip(['C', 'C', 'N', 'N', 'N', 'H', 'C'], [vector([x,0,0]) for x in range(7)])]
mol = molecule(atms)
print(mol.atoms)
move(mol, vector([1,0,1]))
print(mol.atoms)
