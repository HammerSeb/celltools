import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from ctsim.ediff.ewald_diffraction import DiffractionExperiment


def gaussian_2d(x: np.ndarray, y: np.ndarray, A: float, xc: float, yc: float, wx: float, wy: float, offset: float = 0) -> np.ndarray:
    """
    _summary_

    Parameters
    ----------
    x : np.ndarray
        _description_
    y : np.ndarray
        _description_
    A : float
        _description_
    xc : float
        _description_
    yc : float
        _description_
    wx : float
        _description_
    wy : float
        _description_
    offset : float, optional
        _description_, by default 0

    Returns
    -------
    np.ndarray
        _description_
    """
    return A * np.exp( -1 * ( (x - xc)**2 / (2 * wx**2) + (y - yc)**2 / (2 * wy**2)) )



def show_diffraction_experiment(diffexp: DiffractionExperiment, qrange: float, *, kwargs=None) -> [Figure, Axes]:
    """
    show diffraction Experiment using matplotlib's imshow
    
    Parameters
    ----------
    diffexp : DiffractionExperiment
        diffraction Experiment to plot
    qrange : float
        reciprocal space range to plot
    kwargs : _type_
        kwargs of imshow

    Returns
    -------
    Figure, Axes
        figure and axes of plot
    """

    qx, qy = np.meshgrid(
        np.linspace(-1*qrange, qrange, 500),
        np.linspace(-1*qrange, qrange, 500)
    )

    intensity = np.zeros(qx.shape)

    for q, S in zip(diffexp.reflections, diffexp.structure_factor):
        intensity += gaussian_2d(qx, qy, np.abs(S)**2, q[0], q[1], 0.1, 0.1)


    f, ax = plt.subplots()

    ax.imshow(intensity, cmap = "inferno", vmin = 0, vmax = np.percentile(intensity, 98)
              )

    return f, ax
