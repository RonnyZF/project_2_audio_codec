import cmath

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
