from gen import buildrelease
from Queue import Queue
from threading import Thread
from threading import Semaphore
from algo import ParticleFilterBuildRelease as pfbrProcess
from multiprocessing import Pool
import numpy as np
import math
from time import sleep
import matplotlib.pyplot as plt

from algo import fftperiod
from algo import histogram
from algo import autocorrelation

import scipy.stats as stats

"""
class PhaseDriftTest:
    def __init__(self):
        self.eventQ = Queue(maxsize=1)
    
    def generateStreamingEvents(self, series, semaphore, nLevels):
        for i in range(1000):
            event = series.getNext()
            semaphore.acquire(nLevels)
            self.eventQ.put(event)
    def run(self, ratio):
        T = 15.0
        series = buildrelease.StreamingSeries(period=T, phase=0.0, std=T/ratio, 
                                  mode='phaseDrift', rate=0)
        semaphore = Semaphore()
        ts = Thread(target=self.generateStreamingEvents, args=(series, semaphore, 1,)) 
        ts.start()
        pf1 = pfbrProcess.ParticleFilter(100, [0.0, 50.0], self.eventQ, semaphore)
        t1 = Thread(target=pf1.run)
        t1.start()
    def experiment(self):
        Parallel(n_jobs=2)(delayed(self.run)(i) for i in np.linspace(7.5,15,3))
"""
eventQ = Queue(maxsize=1)
def generateStreamingEvents(series, semaphore, nLevels):
    for i in range(1000):
        event = series.getNext()
        semaphore.acquire(nLevels)
        eventQ.put(event)

def plotData(xr, data, filename):
    ############################
    fig = plt.figure()
    ax = fig.add_subplot(111)
    width = 0.35
    ax.bar(xr,         data[:,0], width=width, label="Particle Filter", color='b', edgecolor='black', ecolor='r')
    #TODO: Why this 1.5 times? :(
    ax.bar(xr+1.5*width,   data[:,1], width=width, label="FFT", color='g', edgecolor='black', ecolor='r', align='center')
    ax.bar(xr+2*width, data[:,2], width=width, label="Histogram", color='m', edgecolor='black', ecolor='r')
    ax.bar(xr+3*width, data[:,3], width=width, label="Autocorrelation", color='c', edgecolor='black', ecolor='r')
    ax.set_xticks(xr+1.5*width)
    ax.set_xticklabels(xr)
    ax.set_yscale("log")
    ax.set_xlabel("Drift standard deviation")
    ax.set_ylabel("RMS Error")
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    # Put a legend below current axis
    lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2)
    savefile(filename, lgd)
    #plt.show()

def savefile(filename, lgd=None):
    if (lgd!=None):
        plt.savefig(filename+'.eps', format='eps', dpi=1600, bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.savefig(filename+'.eps', format='eps', dpi=1600, bbox_inches='tight')

def getError(pred_period, true_period):
    x = math.pow(((pred_period - true_period)/true_period),2)
    return x

def runPhaseDrift(std):
    T = 10.0
    print "std:", std
    series = buildrelease.StreamingSeries(period=T, phase=0.0, std=std, 
                              mode='phaseDrift', rate=0, length = T*500)
    events = series.getEvents()
    seq = fftperiod.convertToTimeSeries(events)
    pf1 = pfbrProcess.ParticleFilter(100, [0.0, 20.0], eventQ=None, semaphore=None)
    pfperiod, pf1periods = pf1.run(events)
    pfError = getError(pfperiod, T)
    #print "\tpfError:", pfError 
    fftError = getError(fftperiod.getSignallibPeriod(None, seq), T)
    #print "\tfftError:", fftError
    histError = getError(histogram.getPeriods(None, events), T)
    #print "\thistError:", histError
    autocorrError = getError(autocorrelation.getPeriods(None, seq), T)
    return np.array([pfError, fftError, histError, autocorrError])

def runPeriodDrift(std):
    T = 10.0
    print "std:", std
    series = buildrelease.StreamingSeries(period=T, phase=0.0, std=std, 
                              mode='periodDrift', rate=0, length = std*500)
    events = series.getEvents()
    seq = fftperiod.convertToTimeSeries(events)
    pf1 = pfbrProcess.ParticleFilter(100, [0.0, 20.0])
    pfperiod, pf1periods = pf1.run(events)
    pfError = getError(pfperiod, std)
    return
    #print "\tpfError:", pfError 
    fftError = getError(fftperiod.getSignallibPeriod(None, seq), std)
    #print "\tfftError:", fftError
    histError = getError(histogram.getPeriods(None, events), std)
    #print "\thistError:", histError
    autocorrError = getError(autocorrelation.getPeriods(None, seq), std)
    return np.array([pfError, fftError, histError, autocorrError])

def runPeriodJump():
    T1 = 10
    T2 = 100
    series1 = buildrelease.StreamingSeries(period=T1, phase=0.0, std=T1/3, 
                              mode='phaseDrift', rate=0, length = T1*500)
    series2 = buildrelease.StreamingSeries(period=T2, phase=0.0, std=T2/3, 
                              mode='phaseDrift', rate=0, length = T2*500)
    events = series1.getEvents()
    seq2 = series2.getEvents() + max(events)
    events = np.append(events, seq2)
    pf1 = pfbrProcess.ParticleFilter(100, [0.0, 20.0])
    pf1period, estimatedPeriods1 = pf1.run(events)
    pf2 = pfbrProcess.ParticleFilter(100, [21.0, 200.0])
    pf2period, estimatedPeriods2 = pf2.run(events)
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(estimatedPeriods1, 'g', label="Period bound 0-20")
    ax.plot(estimatedPeriods2, 'r', label="Period bound 21-200")
    ymin, ymax = ax.get_ylim()
    ax.plot(np.ones(10)*(500),np.linspace(ymin, ymax, 10), 'b--')
    ax.set_xlabel("# of samples")
    ax.set_ylabel("# Estimated period")
    """
    errors1 = np.power(np.power((estimatedPeriods1[:501]-T1)/T1, 2), 0.5)
    errors1 = np.append(errors1, np.power(np.power((estimatedPeriods1[500:]-T2)/T2, 2), 0.5))
    errors2 = np.power(np.power((estimatedPeriods2[:501]-T1)/T1, 2), 0.5)
    errors2 = np.append(errors2, np.power(np.power((estimatedPeriods2[500:]-T2)/T2, 2), 0.5))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(errors1, 'g', label="Period bound 0-20")
    ax.plot(errors2, 'r', label="Period bound 21-200")
    ymin, ymax = ax.get_ylim()
    ax.plot(np.ones(10)*(500),np.linspace(ymin, ymax, 10), 'b--')
    ax.set_xlabel("# of samples")
    ax.set_ylabel("# Error")
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    # Put a legend below current axis
    lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2)
    savefile("plot/PeriodJump")
    plt.show()

def experiment():
    """
    print ("##################PERIOD JUMP#######################")
    #runPeriodJump()
    print ("##################PHASE DRIFT#######################")
    xr = np.linspace(0, 12, 7)
    errors = np.zeros((len(xr),4))
    for i in range(10):
        print "Iteration:",i
        errors += np.array(Parallel(n_jobs=8)(delayed(runPhaseDrift)(std) for std in xr))
    errors = np.power(errors/10, 0.5)
    print errors
    plotData(xr, errors, "plot/phaseDrift_drift_vary")
    print ("##################PERIOD DRIFT#######################")
    xr = np.linspace(10, 20, 5)
    errors = np.zeros((len(xr),4))
    for i in range(10):
        print "Iteration:",i
        errors += np.array(Parallel(n_jobs=8)(delayed(runPeriodDrift)(std) for std in xr))
    errors = np.power(errors/10, 0.5)
    print errors
    plotData(xr, errors, "plot/periodDrift_drift_vary")
    """