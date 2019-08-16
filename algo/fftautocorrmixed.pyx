import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from random import randint
from utils import util
import math

def getFFTEnergy(seq):
    fs = 1
    f, Pxx_spec = signal.periodogram(seq, fs, 'hanning', scaling='spectrum')
    Pxx_den = np.sqrt(Pxx_spec)
    return Pxx_den, f

def getRandomSeriesEnergy(length, nOnes):
    nOnes = int(nOnes)
    randS = np.zeros(length+1)
    for i in [randint(0,length) for p in range(nOnes)]:
        randS[i] = 1
    Pxx_den, f = getFFTEnergy(randS)
    return max(Pxx_den)

def crossCorrel(x, y):
    mx = x - x.mean()
    my = y - y.mean()
    denom = math.sqrt(np.sum(np.inner(mx, mx)) * np.sum(np.inner(my,my)))
    sxy = np.sum(mx * my)
    result = sxy/denom
    return result

def autocorr(x):
    """
    result = np.correlate(x, x, mode='full')
    return result[result.size/2:]
    """
    n = len(x)
    #print "n:", n
    result = np.array([crossCorrel(x, x)])
    for delay in range(1, len(x)/2):
        dx =x[delay:]
        y = x[:-delay]
        a = crossCorrel(y, dx)
        #print delay, a
        result = np.append(result, a)
    return result

def isPeak(corr, p):
    print "None"

def getPeriod(sequence):
    baseFreq = np.array([])
    for i in range(10):
        baseFreq = np.append(baseFreq, getRandomSeriesEnergy(len(sequence), sum(sequence)))
    sortedBaseFreq = np.sort(baseFreq)
    threshold = sortedBaseFreq[1]
    #print threshold
    Pxx_den, f = getFFTEnergy(sequence)
    #print "Pxx_den:", Pxx_den
    il = np.where(Pxx_den > threshold)[0]
    ir = il - 1
    pranges = zip(1/f[il], 1/f[ir])
    #print "potential_freqs:", pranges
    corr = autocorr(sequence,)
    #print corr
    maxcorr, period = -1e3, 0
    for p in pranges:
        x = corr[int(math.floor(p[0])):int(math.ceil(p[1]))+1]
        if (maxcorr < max(x)):
            maxcorr = max(x)
            period = np.argmax(x)+p[0]
    #print "period:", period
    return period

"""
def main():
    seq = util.convertToTimeSeries(events)
    getPeriod(seq)


events = np.array([3.31392011e-01,9.41780361e-01,1.00921065e+01,1.54730736e+01,2.04843281e+01,
                   2.07257649e+01,2.38216496e+01,2.66483390e+01,2.95364541e+01,3.60736219e+01,
                   3.91941372e+01,4.82894484e+01,4.90027096e+01,5.84145721e+01,6.36107130e+01,
                   6.85559051e+01,7.71114636e+01,7.89729782e+01,8.34009359e+01,8.35741043e+01,
                   8.47771601e+01,8.74012924e+01,9.39992263e+01,9.80998642e+01,9.87489452e+01,
                   1.05490022e+02,1.08140146e+02,1.11923933e+02,1.18598411e+02,1.29653776e+02,
                   1.34340199e+02,1.40098828e+02,1.49584102e+02,1.51991942e+02,1.59654449e+02,
                   1.69458179e+02,1.72388943e+02,1.76125180e+02,1.76857526e+02,1.80227347e+02,
                   1.87724558e+02,1.88118842e+02,1.89947448e+02,1.91214140e+02,1.91291394e+02,
                   2.01347125e+02,2.06017469e+02,2.11787742e+02,2.20271829e+02,2.20392503e+02,
                   2.23108102e+02,2.30099470e+02,2.30688063e+02,2.30802134e+02,2.32275744e+02,
                   2.32589931e+02,2.34235764e+02,2.41372407e+02,2.42073255e+02,2.50796489e+02,
                   2.53138509e+02,2.53642952e+02,2.59445292e+02,2.62138026e+02,2.71381241e+02,
                   2.81517791e+02,2.90217306e+02,2.92427209e+02,2.95657573e+02,3.02129752e+02,
                   3.07267655e+02,3.07415309e+02,3.10450247e+02,3.12350952e+02,3.15219727e+02,
                   3.22706662e+02,3.30940653e+02,3.38775102e+02,3.49772129e+02,3.51065501e+02,
                   3.52243971e+02,3.59623985e+02,3.68931519e+02,3.77098024e+02,3.79296282e+02,
                   3.79648306e+02,3.86875085e+02,3.89991512e+02,3.99181065e+02,4.09392491e+02,
                   4.09489538e+02,4.11522296e+02,4.11624542e+02,4.12381765e+02,4.18031040e+02,
                   4.26499185e+02,4.37074959e+02,4.47646625e+02,4.47959523e+02,4.56074356e+02,
                   4.58102981e+02,4.61266517e+02,4.67100521e+02,4.76322354e+02,4.86061876e+02,
                   4.87217574e+02,4.98078940e+02])
if __name__=="__main__":
    main()
"""