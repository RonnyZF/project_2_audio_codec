import numpy as np
from py_fft_project2.wave_lib import leer_wave,guardar_wave
from py_fft_project2.plot import plot_fft, plot_signal

factor_escala = 22

*_,data = leer_wave('../audio_samples_8kHz/sample_1-8kHz.wav')
fft_vector = np.fft.fft(data)

def normalize_complex_arr(a):
    a_oo = a - a.real.min() - 1j*a.imag.min() # origin offsetted
    return a_oo/np.abs(a_oo).max()
# aplicando normalizacion
# aplicando escalamiento
fft_vector= fft_vector * (2 ** -factor_escala)

import matplotlib.pyplot as plt
plt.plot(fft_vector.real)
plt.grid()
plt.show()

plt.plot(fft_vector.imag)
plt.grid()
plt.show()

ifft_vector = np.real(np.fft.ifft(fft_vector))

ifft_vector =ifft_vector * (2 ** factor_escala)

guardar_wave(file_name='from_ifft.wav', arreglo=ifft_vector, sampleRate=8000,width=2)
plot_signal(ifft_vector)
plot_fft(ifft_vector)
pass