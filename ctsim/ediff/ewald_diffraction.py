import typing
from typing import List, Tuple, Union, Literal
import numpy as np

from celltools.linalg import Vector
from celltools import Cell


def ewald_diffraction(k_vector: Vector, cell: Cell, *,
                      ds: float = 1e-6, 
                      distance_correction: Union[None, Literal["exponential"], Literal["gaussian"]] = None) -> [List[Vector], List[float]]:
    """_summary_

    Parameters
    ----------
    k_vector : Vector
        _description_
    cell : Cell
        _description_
    ds : float, optional
        _description_, by default 1e-6, keyword-only
    distance_correction : Union[None, Literal[&quot;exponential&quot;], Literal[&quot;gaussian&quot;]], optional
        _description_, by default None, , keyword-only

    Returns
    -------
    [List[Vector], List[float]]
        _description_
    """
    #raise NotImplementedError("Nothing has been done yet")
    
    reflections = None
    structure_factor = None
    return reflections, structure_factor
