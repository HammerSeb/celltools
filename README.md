***LEGACY BRANCH - THIS BRANCH PRESERVER CELLTOOLS AS USED FOR THE ZNPC EXCIMER STUDY.***

# CellTools
This python package provides the tools to build and manipulate unit cells and supercells and to calculate diffraction 
patterns from these cells. Using ```pyqtgraph``` the cells and supercells can be drawn in 3D. 

The package has been developed as part of the project "Analyzing the Intermolecular Dynamics of Excited States in Molecular Semiconductors" funded by the German Research Foundation (DFG) within in the Walter Benjamin Program (project 490894053) which is greatfully acknowledged. 

***This package is under active development and collaborates are welcome. I will try to ensure consistency so that old code
can run on newer version, but cannot guarantee it at this point.*** 

*At some point a readthedocs will also come up, so far you have to work with the docstrings of the functions, 
unfortunately. All classes and functions are well documented.*



## Installation
Either download the whole package from here and install into any environment by executing 
```shell
python3 -m pip install -e celltools
```
inside the download folder, or use
```shell
python3 -m pip install git+https://github.com/HammerSeb/celltools.git
```
to install directly from Github. 

### Known issues
For drawing the unit cells the package uses the opengl module of ```pyqtgraph```, which can be hard to get to work 
sometimes. Before you start, make sure that you can run all the 3D examples of pyqtgraph. Run
```shell
python3 -m pyqtgraph.examples
```
to start the example module and try all 3D examples. Unfortunately I cannot help you if you have problems there. What 
has been known to work is to install ```pyqt5```, ```pyopengl``` and ```pyqtgraph``` into a clean environment and start 
from there. Good Luck!

*For anaconda, loading errors of iris or swrast can sometimes be fixed by running*
```shell
conda install -c conda-forge libstdcxx-ng
```

## How do I use it
The package comes in two parts, ```celltools``` and ```simulation```. The former contains all tools to build a unit cell
and manipulate it, while the latter contains diffraction and pair distribution function simulations. 

### celltools
```celltools``` cotains a linear algebra part ```linalg```, a unit and supercell part ```cell``` and a drawing part ```
draw```. A quick explanation follows of the important features follows

#### linalg
```Basis```: class describing a basis system spanned by three linear independent vectors
```python
from celltools.linalg import Basis
standard_basis = Basis([1, 0, 0], [0, 1, 0], [0, 0, 1])
```

```Vector```: class representing a vector defined by coordinates in a given ```Basis```
```python
from celltools.linalg import Basis, Vector
basis = Basis([2, 0, 0], [1,1,0], [1, 0, -3])
v = Vector([1,0,0], basis)
w = Vector([1,1,1], basis)

z = 3*v + w
print(f"coordinates: {z.vector}")
print(f"global coordinates: {z.global_coord}, global length: {z.abs_global}")
```
vector addition is possible for vectors in the same basis. Scaler multiplication is possible as well. Coordinates in 
basis system, coordinates  and length in global reference frame can be returned by ```Vector.vector```, 
```Vector.global_coord``` and ```Vector.abs_global```, respectively. 
*If no basis is given, vectors are defined in the standard basis.*

```Line``` and ```Plane```: Classes representing lines and planes in 3D space. They are defined by a point of origin and
a directional/normal vector. Can be used to define rotational axis for example. 
```python
from celltools.linalg import Vector, Line, Plane

origin = Vector([1, 1, 1])
vec = Vector([1, 2, 3])

# Line through origin in direction to vec
line = Line(origin, vec)

# Plane through origin, normal to vec
plane = Plane(origin, vec)
```

#### cell
```Latice```, ```Atom``` and ```Cell```: Classes representing lattice vectors, atoms and unit cells. 

```python
from celltools.linalg import Basis, Vector
from celltools.cell import Lattice, Atom, Cell

#Aluminum lattice
lattice = Lattice(Basis([4.59, 0, 0], [0, 4.59, 0], [0, 0, 4.59]))
origin = Vector([0, 0, 0], lattice)
a, b, c = Vector([1, 0, 0], lattice), Vector([0, 1, 0], lattice), Vector([0, 0, 1], lattice)
#Atom positions in an fcc lattice
atoms = [Atom('Al', origin), Atom('Al', a), Atom('Al', b), Atom('Al', c), #corners
         Atom('Al', 0.5*(a+b)), Atom('Al', 0.5*(a+c)), Atom('Al', 0.5*(b+c)), # faces
         ]

al_unit_cell = Cell(lattice, atoms)
```
```SuperCell```: Class contructing a super cell of given size from a ```Cell``` object.
```python
from celltools import SuperCell

supercell = SuperCell(al_unit_cell, (3,3,3))
```

```Molecule```: Container class for atoms, can also be content of a cell
```python
from celltools.linalg import Vector
from celltools.cell import Atom, Molecule
#Atoms of CO2 molecule along x-axis of standard basis
atoms = [Atom('C', Vector([0,0,0])), Atom('O', Vector([-1.16,0,0])), Atom('O', Vector([1.16,0,0]))]
#Define molecule
carbondioxide = Molecule(atoms)
carbondioxide.auto_bonds() #add bonds between atoms automatically
```

#### Moving and Rotation
There are explicit functions to move atoms and molecules in ```celltools.cell.tools```. Using the CO2 molecule from 
above as an example shows how it's done:

```python
from celltools.linalg import Line
from celltools.cell import move, rotate
#Move atom 2 Angstrom along x-axis
move(carbondioxide, Vector([2,0,0]))

#Rotate molecule, now centered at (2,0,0) 45 degrees around its z-axis
zaxis = Line(Vector([2,0,0]), Vector([0,0,1]))
rotate(carbondioxide, zaxis, 45, mode="deg")
```

#### Draw
Drawing is done using ```pyqtgraph```. There are a number of function provided in ```draw``` to facilitate the 3D 
depiction of the cells. We will give an example with the aluminum supercell from above
```python
import pyqtgraph as pg
from celltools import make_figure, draw_cell

w = make_figure()
draw_cell(w, al_unit_cell)

pg.exec()
```
Using ```draw_supercell``` instead of ```draw_cell``` draws a given super cell. With ```draw_line``` and ```draw_plane```
Line and Plane objections can be drawn, to e.g. show molecular planes or rotation axis.

#### import and export
Instead of adding each atom to the unit cell yourself, cif-files or ```crystals.Crystal``` objects can be used to 
generate ```Cell``` objects. Usually not all atoms in the unit cell are listed in the cif file, but only the irreducible
basis is listed and the rest is generated from the space group symmetry. For the ```cell_from_cif``` function, there
is the keyword argument ```mode``` which accepts ```"file"```, for which only the listed atoms of the cif-file are added
, and ```"sym"```, for which the rest of the atoms is generated from the crystal symmetry. If the symmetry operations 
provided in the cif-file, the function will use these, otherwise, if the space group has been implemented, the 
necessary operations are taken from ```SPACE_GROUP``` dictionary in ```celltools\cell\spacegroup_data.py```. 

*So far, there's no export option for the ```Cell``` or ```SuperCell``` class. Needs to be implemented*

```python
from celltools import cell_from_cif
from celltools.cell.generate import cell_from_crystal
from crystals import Crystal

cell_cif = cell_from_cif("./your_file.cif")

crystal = Crystal() #dummy Crystal object
cell_crystal = cell_from_crystal(crystal)
```


### simulation
So far only the **kinematic** electron diffraction simulation has been implemented. *I plan to implement electron pair distribution 
function and x-ray diffraction simulations as well.* 

Diffraction can either be simulated from a ```Cell``` or a ```Supercell```. The scattering vectors as a ```Vector``` and
the complex scattering amplitude can be simulated from a list of miller (hkl) indices. The package also provides a 
powder diffraction simulation. 

```python
import numpy as np
from celltools import cell_from_cif
from simulation.ediff import diffraction_from_cell, powder_pattern

cell = cell_from_cif("./your_file.cif")

# simulate diffraction for given miller indices
hkl = [(1,0,0), (1,-1,0), (0,2,0), (0,0,3)] #exmplary miller indices
scatter_vec, amlpitude = diffraction_from_cell(hkl, cell)

# simulate poweder instensity of given q-range
q = np.linspace(0,3,200)
q_scatt = [vec.abs_global for vec in scatter_vec]
intensity = powder_pattern(q, q_scatt, amlpitude, w=0.03)
```
