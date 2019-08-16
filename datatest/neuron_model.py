from brian2 import *
import matplotlib.pyplot as plt
import math

def convertToTimeSeries(events, resolution=1):
    n = int(math.ceil(max(events)*resolution)+1)
    s = np.zeros(n)
    for i in range(0, len(events)):
        e = events[i]
        x = int(math.floor(e*resolution))
        s[x] = 1
    return s

def convertToEvents(allevents):
    d = np.diff(allevents)
    threshold = np.percentile(d, 98)
    #for dx in range(len(d)):
    #    print dx, allevents[dx], d[dx]
    #print "threshold:", threshold
    indexes = np.where(d>threshold)
    #print allevents[indexes]
    return allevents[indexes]

def getEvents():
    return simulateNeuron()

def simulateNeuron():
    N = 1000
    Vr = 10*mV
    theta = 20*mV
    tau = 20*ms
    delta = 2*ms
    taurefr = 2*ms
    duration = 1*second
    C = 1000
    sparseness = float(C)/N
    J = .1*mV
    muext = 25*mV
    sigmaext = 1*mV
    eqs = """
    dV/dt = (-V+muext + sigmaext * sqrt(tau) * xi)/tau : volt
    """
    
    group = NeuronGroup(N, eqs, threshold='V>theta',
                        reset='V=Vr', refractory=taurefr, method='euler')
    group.V = Vr
    conn = Synapses(group, group, on_pre='V += -J', delay=delta)
    conn.connect(p=sparseness)
    M = SpikeMonitor(group)
    LFP = PopulationRateMonitor(group)
    
    run(duration)
    
    events = M.t/ms
    es = convertToEvents(events)
    #seq = convertToTimeSeries(es)
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.bar(range(len(seq)), seq)
    #plt.show()
    return es
    """
    seq = convertToTimeSeries(events)
    plt.subplot(211)
    plt.bar(range(len(seq)), seq)
    #plt.plot(M.t/ms, M.i, '.')
    #plt.xlim(0, duration/ms)
    
    plt.subplot(212)
    plt.plot(LFP.t/ms, LFP.smooth_rate(window='flat', width=0.5*ms)/Hz)
    plt.xlim(0, duration/ms)
    
    plt.show()
    """
if __name__ == "__main__":
    events = simulateNeuron()
    np.save("neuron_data/events", events)