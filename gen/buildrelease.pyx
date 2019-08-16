import scipy.stats as stats
import numpy as np
import math, scipy, os
from time import time
import utils.util as util
import matplotlib.pyplot as plt

class StreamingSeries:
    lasteventIndex = 0
    def generateevnets(self, length, lastPoint=0):
        if (self.mode == 'periodDrift'):
            periodicevents = self.periodDrift(length) + lastPoint
        else:
            periodicevents = self.phaseDrift(length) + lastPoint
        #print "periodicevents:", periodicevents, len(periodicevents)
        noiseevents = self.backgroundNoise(int(math.floor(max(periodicevents)))) + lastPoint
        #print "noiseevents:", noiseevents, len(noiseevents)
        self.events = np.unique(np.sort(np.append(periodicevents, noiseevents)))
        
    def __init__(self, period, phase, std, mode='periodDrift', rate=0, length=1000, distribution='normal'):
        self.period = period
        self.phase = phase
        self.std = std
        self.rate = rate
        self.mode = mode
        self.distribution = distribution
        
        self.generateevnets(length)

    def updatePeriod(self, period):
        self.period = period

    def getEvents(self):
        return self.events
    
    def periodDrift(self, n):
        fs = np.array([])
        fs = np.append(fs, stats.uniform.rvs())
        while(True):
            norm = stats.norm(loc = self.phase + fs[-1], scale=self.std)
            nextx = util.boundedSample(norm, fs[-1], 1e9)
            if (nextx > n):
                break
            if (fs[-1] != n):
                fs = np.append(fs, nextx)
        return fs

    def phaseDrift(self, n):
        fs = np.array([])
        fs = np.append(fs, stats.uniform.rvs())
        while(True):
            if (self.distribution == "normal"):
                rrv = stats.norm(loc = self.phase + self.period + fs[-1], scale=self.std)
            elif (self.distribution == "laplace"):
                rrv = stats.laplace(loc = self.phase + self.period + fs[-1], scale=self.std)
            elif (self.distribution == "gamma"):
                rrv = stats.gamma(2, loc = self.phase + self.period + fs[-1], scale=self.std)
            nextx = util.boundedSample(rrv, fs[-1], 1e9)
            if (nextx > n):
                break
            if (nextx >= nextx):
                fs = np.append(fs, nextx)
        return fs

    def backgroundNoise(self, n):
        if (self.rate == 0):
            return np.array([])
        fs = np.array([np.random.uniform()])
        expon = stats.expon(loc=0, scale=1.0/self.rate)
        while True:
            x = fs[-1] + expon.rvs()
            if (x > n):
                break
            fs = np.append(fs, x)
        return fs

    def getNext(self):
        if (self.lasteventIndex < len(self.events)):
            event = self.events[self.lasteventIndex]
            self.lasteventIndex += 1
        else:
            self.generateevnets(lastPoint=max(self.events))
            event = self.events[0]
            self.lasteventIndex = 1
        return event

class test:
    def plotHistogram(self, events):
        d = np.diff(events)
        n, bins, patches = plt.hist(d, 50, normed=1, 
                                    facecolor='green', alpha=0.75)
        plt.show()

    def periodDriftTest(self):
        brProcess = StreamingSeries(period=2.0, phase=0.0, std=1.0, 
                                 mode='periodDrifft', rate=0)
        events = brProcess.getEvents()
        self.plotHistogram(events)
    def phaseDriftTest(self):
        brProcess = StreamingSeries(period=2.0, phase=0.0, std=1.0, 
                                 mode='phaseDrifft', rate=0)
        events = brProcess.getEvents()
        self.plotHistogram(events)