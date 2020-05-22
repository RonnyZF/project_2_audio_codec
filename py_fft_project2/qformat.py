"""
x is the input fixed number which is of integer datatype
e is the number of fractional bits for example in Q1.15 e = 15
"""
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

if __name__ == "__main__":
    from bitstring import BitStream, BitArray
    import math
    import wave, struct, math, random
    from collections import namedtuple
    coeficiente = namedtuple('coeficiente', 'real imaginario')


    def escribir_binario(filename, data_h, data_l):
        with open(filename, 'ab') as f:
            data = struct.pack('<hh', data_h, data_l)
            f.write(data)

    def leer_binario(filename):
        with open(filename, 'rb') as f:
            data = f.read()
        return struct.iter_unpack('<hh',data)
    n = [1,2,3]
    i = 12
    for nn in n:
        coef_towrite= coeficiente._make([to_fixed(f=nn,e=i), to_fixed(f=-nn,e=i)])
        escribir_binario(filename='file.mpx', data_h=coef_towrite.real,data_l=coef_towrite.imaginario)

    dato = leer_binario(filename='file.mpx')
    coef = [coeficiente._make(i) for i in dato]
    pass

