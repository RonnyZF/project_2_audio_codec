# imports
import math
import numpy as np
from first_deliverable.wave_lib import leer_wave,guardar_wave
from first_deliverable.plot import plot_fft, plot_signal
from first_deliverable.array_store_read_binary import store_coeffs, read_coeffs

# constantes
coder_level = 'float16'

settings = {
    'float16': {'factor_escala': 22,'fixed_e' : '','width_coded':'<ee'}
           }
factor_escala = settings[coder_level]['factor_escala']
fixed_e = settings[coder_level]['fixed_e']
width_coded = settings[coder_level]['width_coded']

*_,data = leer_wave('../audio_samples_8kHz/sample_1-8kHz.wav')
data=np.array(data)
size = math.log(data.shape[0],2)
# recorte para muestras de tamaño de potencias de 2
if size % 2 != 0:
    to_2_pow = math.floor(size)
    data=data[:2 ** to_2_pow]
fft_vector = np.fft.fft(data)

# aplicando escalamiento
fft_vector= fft_vector * (2 ** -factor_escala)

import matplotlib.pyplot as plt
plt.plot(fft_vector.real)
plt.grid()
plt.show()

plt.plot(fft_vector.imag)
plt.grid()
plt.show()

h, _ = np.split(fft_vector,2)
store_coeffs(filename='file.mpx', dato=h,width_coded=width_coded)
recovered = read_coeffs(filename='file.mpx',width_coded=width_coded)
recovered = np.concatenate((recovered,recovered[::-1].conj()))

ifft_vector = np.real(np.fft.ifft(recovered))
ifft_vector =ifft_vector  * (2 ** (factor_escala))

guardar_wave(file_name='from_ifft.wav', arreglo=ifft_vector, sampleRate=8000,width=2)
