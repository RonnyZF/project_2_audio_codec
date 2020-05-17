from py_fft_project2.py_fft_list import fft1
import numpy as np

def test_fft1():
    # Two sine waves, 40 Hz + 90 Hz
    t = np.linspace(0, 0.5, 2 ** 10)
    s = np.sin(40 * 2 * np.pi * t) + 0.5 * np.sin(90 * 2 * np.pi * t)

    assert all(np.isclose([abs(i) for i in np.fft.fft(s)], [abs(i) for i in fft1(s)]) )