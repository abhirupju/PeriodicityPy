import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

def getPeriods(ax, events):
    d = np.diff(events)
    hist, bins = np.histogram(d, bins=100)
    if (ax != None):
        ax.bar(bins[:-1], hist)
    return bins[np.argmax(hist)]
    #return stats.mode(d)[0][0]

def getOnlinePeriods(events):
    d = 0
    periods = np.array([])
    for i in range(1, len(events)):
        periods = np.append(periods, getPeriods(None, events[:i]))
    return periods