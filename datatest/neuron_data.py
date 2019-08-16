import numpy as np
import matplotlib.pyplot as plt
import math

def getax():
    fig = plt.figure()
    return fig.add_subplot(111)

def convertToTimeSeries(events, resolution=1):
    n = int(math.ceil(max(events)*resolution)+1)
    s = np.zeros(n)
    for i in range(0, len(events)):
        e = events[i]
        x = int(math.floor(e*resolution))
        s[x] = 1
    return s

def analyze(data):
    ax = getax()
    d = data[:,0][:100]
    print d*1000
    s = convertToTimeSeries(d, 1000)
    ax.bar(range(len(s)), s)
    plt.show()

if __name__ == "__main__":
    analyze(np.load("/home/abhirup/Documents/Research/dataset/neuron/data1.npy"))