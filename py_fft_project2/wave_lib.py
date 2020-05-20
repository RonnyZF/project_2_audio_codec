import wave, struct

def guardar_wave(file_name, arreglo=None, sampleRate=44100):
    obj = wave.open(file_name,'w')
    obj.setnchannels(1) # mono
    obj.setsampwidth(2) # 2 Bytes = 16 bits
    obj.setframerate(sampleRate)
    for i in arreglo:
       data = struct.pack('<h', i)
       obj.writeframesraw( data )
    obj.close()

def leer_wave(file_name):
    waveFile = wave.open(file_name, 'r')
    length = waveFile.getnframes()
    channels = waveFile.getnchannels()
    width = waveFile.getsampwidth()
    frame_rate = waveFile.getframerate()
    number_of_frames = waveFile.getnframes()
    parameters = waveFile.getparams()
    wave_list = []
    for i in range(0,length):
        waveData = waveFile.readframes(1)
        data = struct.unpack("<h", waveData)
        wave_list.append(data[0])
    waveFile.close()
    return channels,width,frame_rate,number_of_frames,parameters,wave_list

wav = leer_wave(file_name='leer.wav')

guardar_wave('salvar.wav',arreglo=wav,sampleRate=8000)