import pathlib
from os.path import join
from re import findall
from celltools.cell.tools import SymmetryOperator, create_SymmetryOperator

def read_sym_file(file):
    """

    Parameters
    ----------
    file

    Returns
    -------

    """

    path_to_sg_symmetries = join(pathlib.Path(__file__).parent.resolve(),"sg_symmetries")

    sg_operator = []
    with open(join(path_to_sg_symmetries,file), "r") as symfile:
        for line in symfile.readlines():
            sg_operator.append(
                create_SymmetryOperator(findall("\S{1,20},\S{1,20},\S{1,20}",line)[0])
            )
    return sg_operator

SPACE_GROUP= {
    "1": read_sym_file("1.sym"),
    "2": read_sym_file("2.sym"),
}