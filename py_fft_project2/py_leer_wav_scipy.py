import wave
import os

dir = os.path.realpath(os.path.normpath(os.path.join(os.path.dirname(__file__), r'../', 'audio_samples_8kHz')))
filename = 'sample_1-8kHz.wav'
with wave.open (os.path.join(dir,filename)) as f:
    pass
pass
