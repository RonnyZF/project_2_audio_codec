from py_fft_project2.py_fft_np import FFT_vectorized

import numpy as np

def test_FFT_vectorized():
    x = np.random.random(1024)
    np.allclose(FFT_vectorized(x), np.fft.fft(x))
