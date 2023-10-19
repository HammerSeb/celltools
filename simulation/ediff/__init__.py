
from typing import List, Tuple
import numpy as np
from .diffraction_from_cell import diffraction_from_cell
from .diffraction_from_supercell import diffraction_from_supercell



def _gaussian(x: np.ndarray, A: float, w: float, xc: float) -> np.ndarray:
    return A / (w * np.sqrt(2 * np.pi)) * np.exp(- np.square(x - xc) / (2 * w **2))

def powder_pattern(q: np.ndarray, scatt_vec: List[float], amp: List[complex], w: float = 0.05) -> np.ndarray:
    """
    simulate a powder pattern from a list of scattering vectors and complex scattering amplitudes. A gaussian peak is
    used to imitate finite size broadening
    Parameters
    ----------
    q: np.ndarray
        q-range over which to simulate the powder pattern
    scatt_vec: list of floats
        list of scattering vectors
    amp: list of complex
        list of complex amplitudes at scattering vectors
    w: float (default 0.05)
        standard deviation of gaussian peaks

    Returns
    -------
        np.ndarray size of q
            powder pattern over q
    """
    intensity = np.zeros_like(q)
    for _scatt_vec, _amp in zip(scatt_vec, amp):
        intensity += _gaussian(q, np.abs(_amp)**2, w, _scatt_vec)
    return intensity
