from py_fft_project2.py_fft_iterative import fft_basic
from py_fft_project2.py_fft_iterative import complex_dft
import numpy as np
import math

def test_complex_dft():
    # Two sine waves, 40 Hz + 90 Hz
    n = 2 ** 10
    t = np.linspace(0, 0.5, n)
    s = np.sin(40 * 2 * np.pi * t) + 0.5 * np.sin(90 * 2 * np.pi * t)

    xr = s[:]
    xi = s[:]
    rex, imx = complex_dft(xr, xi, n)


    np_fft = np.fft.fft(s)

    assert all(np.isclose([abs(i) for i in np.fft.fft(s)], [abs(i) for i in fft_basic(s)]))