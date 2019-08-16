import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import get_window
import time, math
from scipy import signal

def power_two(n):
    return math.ceil((math.log(n, 2)))

def periodogram(ax, x):
    "Argument x is a NumPy array containing the time series data"
    n = len(x)
    I_w = np.abs(np.fft.fft(x))**2 / n
    #w = n / (2 * np.pi * np.arange(n))      # Fourier frequencies
    w = 1.0/np.fft.fftfreq(n,1)
    w, I_w = w[:int(n/2)], I_w[:int(n/2)]  # Truncate to interval [0, pi]
    ax.plot(w, I_w, 'b-', lw=2, alpha=0.5, label='periodogram')
    #return w, I_w

"""
def getNumpylibPeriod(ax, seq, isevents=False):
    if (isevents == True):
        s = convertToTimeSeries(seq)
    else:
        s = seq
    N = len(s)
    if (N < 4):
        return 1
    W = abs(np.fft.fft(s))
    freq = np.fft.fftfreq(N,1)
    periods = 1.0/freq[1:N/2]
    if (N == 2):
        return 1
    else:
        idx = np.argmax(W[1:N/2])
    period = periods[idx]
    if (ax != None):
        ax.plot(periods, abs(W[1:N/2]))
        #plt.scatter([1/max_f,], [np.abs(W[idx]),], s=100,color='r')
        ax.set_xlabel(r"$1/f$")
        ax.set_title("T0:"+str(period))
    return period
"""

def getSignallibPeriod(ax, s):
    fs = 1
    f, Pxx_spec = signal.periodogram(s, fs, 'hanning', scaling='spectrum')
    Pxx_den = np.sqrt(Pxx_spec)
    freq = f[np.argmax(Pxx_den[1:])+1]
    period = 1.0/freq
    #print "FFT period:", period
    if (ax != None):
        #ax.semilogy(f, np.sqrt(Pxx_spec))
        ax.plot(1/f[1:], Pxx_den[1:], color='black')
        ax.set_xlabel('Potential Period')
        ax.set_ylabel('Energy')
        ax.set_title('FFT Periodogram T0:'+str(period))
    return period