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

def to_float(x,e):
    c = abs(x)
    sign = 1
    if x < 0:
        # convert back from two's complement
        c = x - 1
        c = ~c
        sign = -1
    f = (1.0 * c) / (2 ** e)
    f = f * sign
    return f

"""
f is the input floating point number 
e is the number of fractional bits in the Q format. 
    Example in Q1.15 format e = 15
"""
def to_fixed(f,e):
    a = f* (2**e)
    b = int(round(a))
    if a < 0:
        # next three lines turns b into it's 2's complement.
        b = abs(b)
        b = ~b
        b = b + 1
    return b

def escribir_binario(filename, data_h, data_l, width_coded='<ii'):
    with open(filename, 'ab') as f:
        # (-32768) <= number <= 32767
        g = lambda x : x if ((-32768) <= x <= 32767) else (np.sign(x)*32767)
        data = struct.pack(width_coded, g(data_h), g(data_l))
        f.write(data)

def leer_binario(filename, width_coded='<ii'):
    with open(filename, 'rb') as f:
        data = f.read()
    return struct.iter_unpack(width_coded,data)

def store_coeffs(filename, dato,fixed_e = 16,width_coded='<ii'):
    if path.exists(filename):
        os.remove(filename)
    n = [coeficiente._make((i.real, i.imag)) for i in dato]
    for nn in n:
        coef_towrite= coeficiente._make([to_fixed(f=nn.real,e=fixed_e), to_fixed(f=nn.imaginario,e=fixed_e)])
        escribir_binario(filename=filename, data_h=coef_towrite.real,data_l=coef_towrite.imaginario,width_coded=width_coded)

def read_coeffs(filename, fixed_e = 16,width_coded='<ii'):
    dato = leer_binario(filename=filename,width_coded=width_coded)
    coef = [coeficiente._make(i) for i in dato]
    coef_float = [coeficiente._make([to_float(i.real, fixed_e), to_float(i.imaginario, fixed_e)]) for i in coef]
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
