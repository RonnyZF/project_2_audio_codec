"""
x is the input fixed number which is of integer datatype
e is the number of fractional bits for example in Q1.15 e = 15
"""
from bitstring import BitStream, BitArray
import math
import os
import numpy as np
import os.path
from os import path
import wave, struct, math, random
from collections import namedtuple

coeficiente = namedtuple('coeficiente', 'real imaginario')

def escribir_binario(filename, data_h, data_l, width_coded='<ii'):
    with open(filename, 'ab') as f:
        data=struct.pack(width_coded, data_h, data_l)
        f.write(data)

def leer_binario(filename, width_coded='<ii'):
    with open(filename, 'rb') as f:
        data = f.read()
    return struct.iter_unpack(width_coded,data)

def store_coeffs(filename, dato,width_coded='<ii'):
    if path.exists(filename):
        os.remove(filename)
    n = [coeficiente._make((i.real, i.imag)) for i in dato]
    for nn in n:
        coef_towrite= coeficiente._make([nn.real, nn.imaginario])
        escribir_binario(filename=filename, data_h=coef_towrite.real,data_l=coef_towrite.imaginario,width_coded=width_coded)

def read_coeffs(filename, width_coded='<ii'):
    dato = leer_binario(filename=filename,width_coded=width_coded)
    coef = [coeficiente._make(i) for i in dato]
    coef_float = [coeficiente._make([i.real, i.imaginario]) for i in coef]
    return np.array([i.real + i.imaginario*1j for i in coef_float])

if __name__ == "__main__":
    """
    Para verificar el funcionamiento de del conversor de matrices de punto flotante
    a punto fijo de e=10
    La tolerancia obtenida es alrededor de 1e-3
    """
    from numpy.random import rand
    A = rand(100, 2).view(dtype=np.complex128)
    n = np.array([i[0] for i in A])
    store_coeffs(filename='file.mpx', dato=n)
    recovered = read_coeffs(filename='file.mpx')

    assert all(np.isclose(n,recovered, atol=1e-05))

    pass

