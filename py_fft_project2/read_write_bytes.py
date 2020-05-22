import numpy as np

def complex_array_binary_store(file_name, data):
    array = []
    with open(file_name, 'w+b') as f:
        for number in data:
            array.extend([int(number.real), int(number.imag)])
        binary_format = bytearray(array)
        f.write(binary_format)

def complex_array_binary_read(file_name):
    data = []
    convert = lambda x: int.from_bytes(x, "little")
    with open(file_name, 'r+b') as f:
        for ch in iter(lambda: f.read(1), ""):
            if ch == b'':
                break
            data.append(convert(ch))
    return np.array(data).reshape((len(data)//2),2)

Z = np.random.randint(low=0,high=(2**8)-1,size=10) + np.random.randint(low=0,high=(2**8)-1,size=10) * 1j

complex_array_binary_store(file_name='my_file.bin',data=Z)
complex_array_binary_read(file_name='my_file.bin')