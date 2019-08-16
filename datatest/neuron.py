import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from utils import util
from algo import segment_periodicity_gcd
from algo import ParticleFilterBuildRelease as pfbrProcess

def getax():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    return ax

def plotSingleNeuron(events):
    seq = util.convertToTimeSeq(events, resolution=1)
    ax = getax()
    ax.bar(range(len(seq)), seq, align='center', color='black')

def loadData(filename):
    data = np.load(filename)
    ns = np.unique(data[:,1])
    arr = np.array([])
    for d in data:
        if (d[1] == 15):
            arr = np.append(arr, d[0])
    return arr
        
if __name__ == "__main__":
    events = loadData("/home/abhirup/Documents/Research/dataset/neuron/data1.npy")*100
    #plotSingleNeuron(events[:100])
    print "len(events):", len(events)
    pf = pfbrProcess.ParticleFilter(100, [0.0, 20.0],showgraph=True)
    pf.run(events)
    #segment_periodicity_gcd.getPeriods(getax(), events)
    #seq = util.convertToTimeSeq(events, 1)
    #fftperiod.getSignallibPeriod(getax(), seq)
    plt.show()