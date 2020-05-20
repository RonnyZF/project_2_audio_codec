from scipy.fftpack import fft
import numpy as np
import matplotlib.pyplot as plt

def plot_signal(data, N = 80,T = 1.0 / 8000.0):
    """

    :param data:
    :param N: Number of samplepoints
    :param T: sample spacing
    :return:
    """
    yf = np.array(data)
    xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
    plt.plot(xf, yf[0:N//2])
    plt.grid()
    plt.show()

def plot_fft(data, N = 80000,T = 1.0 / 8000.0):
    """

    :param data:
    :param N: Number of samplepoints
    :param T: sample spacing
    :return:
    """
    yf = fft(data)
    xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
    plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))
    plt.grid()
    plt.show()