import numpy as np

##########################################
# Cooley–Tukey FFT algorithm
##########################################
def FFT(x):
    N = len(x)
    if N <= 1: return x
    X_even = FFT(x[::2])
    X_odd = FFT(x[1::2])
    T= [np.exp(-2j*np.pi*k/N)*X_odd[k] for k in range(N//2)]
    return [X_even[k] + T[k] for k in range(N//2)] + [X_even[k] - T[k] for k in range(N//2)]

##########################################
# Test data
##########################################
x = np.random.random(1024)

##########################################
# Output
##########################################
print(np.allclose(FFT(x), np.fft.fft(x)))