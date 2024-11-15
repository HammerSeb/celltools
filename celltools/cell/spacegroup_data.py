import typing
import pathlib
from os.path import join
from re import findall
from celltools.cell.tools import SymmetryOperator, create_SymmetryOperator


def read_sym_file(file: str) -> SymmetryOperator:
    """
    This function reads a symmetry file in celltools/cell/sg_symmetries and returns a SymmetryOperator. It's an
    auxiliary function that is used in cell_from_cif if no symmetry is provided in the cif-file.
    Parameters
    ----------
    file: str
        symmetry file

    Returns
    -------
        SymmetryOperator
    """

    path_to_sg_symmetries = join(
        pathlib.Path(__file__).parent.resolve(), "sg_symmetries"
    )

    sg_operator = []
    with open(join(path_to_sg_symmetries, file), "r") as symfile:
        for line in symfile.readlines():
            sg_operator.append(
                create_SymmetryOperator(findall(r"\S{1,20},\S{1,20},\S{1,20}", line)[0])
            )
    return sg_operator


SPACE_GROUP = {
    """ This dictionary points to all the symmetry files for the respective space group that have been implemented"""
    "1": "1.sym",
    "2": "2.sym",
    "14": "14.sym",
}
