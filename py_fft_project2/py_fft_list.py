import cmath
import numpy

# Handmade implementation of the Cooley-Tukey algorithm, from here: https://jeremykun.com/2012/07/18/the-fast-fourier-transform/
# Make sure that input signal is a power of two!
def omega(p, q):
  return cmath.exp((2.0 * cmath.pi * 1j * q) / p)

def fft1(signal):
  n = len(signal)
  if n % 1 > 0:
      raise ValueError("size of x must be a power of 2")
  if n == 1:
     return signal
  else:
     Feven = fft1([signal[i] for i in range(0, n, 2)])
     Fodd = fft1([signal[i] for i in range(1, n, 2)])
     combined = [0] * n
     for m in range(n // 2):
        combined[m] = Feven[m] + omega(n, -m) * Fodd[m]
        combined[m + n // 2] = Feven[m] - omega(n, -m) * Fodd[m]
     return combined

# Simple samples
s = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]

fft_a = fft1(s)
fft_b = numpy.fft.fft(s)

# Same
print(fft_a[3].real)
print(fft_b[3].real)

# Two sine waves, 40 Hz + 90 Hz
t = numpy.linspace(0, 0.5, 2**10)
s = numpy.sin(40 * 2 * numpy.pi * t) + 0.5 * numpy.sin(90 * 2 * numpy.pi * t)

fft_a = fft1(s)
fft_b = numpy.fft.fft(s)

# Not the same
print(fft_a[3].real)
print(fft_b[3].real)