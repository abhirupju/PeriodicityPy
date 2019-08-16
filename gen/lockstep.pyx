import scipy.stats as stats
import numpy as np
import math

def boundedSample(rv, base, mmin, mmax):
    xx = rv.rvs() + base
    tryCount = 0
    while(xx < mmin or xx > mmax):
        tryCount -= 1
        if (tryCount < 0):
            return -1
        xx = xx = rv.rvs() + base
    return xx

class StreamingSeries:
    lastPeriodicEvent = 0
    segNum = 0
    lasteventIndex = 0
    
    def generateevnets(self, lastPoint=0):
        periodicevents = self.obsNoise(500) + lastPoint
        noiseevents = self.backgroundNoise(int(math.floor(max(periodicevents)))) + lastPoint
        self.events = np.unique(np.sort(np.append(periodicevents, noiseevents)))
        
    def __init__(self, period, phase, std, rate=0):
        self.period = period
        self.phase = phase
        self.std = std
        self.rate = rate
        self.generateevnets()
    def updatePeriod(self, period):
        self.period = period
    
    #def createCombinedSeries(self):
    def getNextPeriodic(self):
        if (self.segNum == 0):
            event = stats.uniform().rvs()
        else:
            rdn = stats.norm(loc=0, scale=self.std)
            event = boundedSample(rdn, self.phase+
                                       (self.segNum * self.period), 
                                       self.lastPeriodicEvent, 1e10)
        self.segNum += 1
        self.lastPeriodicEvent = event
        return event
    def obsNoise(self, nt):
        fs = np.array([])
        s = np.arange(self.period, (nt+1)*self.period, self.period) + self.phase
        norm = stats.norm(loc=0, scale=self.std)
        for i in range(nt):
            x = s[i] + norm.rvs()
            while (x < 0):
                x = s[i] + norm.rvs()
            if (i>0):
                while (x < s[i-1]):
                    x = s[i] + norm.rvs()
            fs = np.append(fs, x)
        return fs

    def backgroundNoise(self, n):
        fs = np.array([np.random.uniform()])
        expon = stats.expon(loc=0, scale=1.0/self.rate)
        while True:
            x = fs[-1] + expon.rvs()
            trycount = 10
            while (x > n):
                trycount -= 1
                if (trycount <= 0):
                    return fs
                x = fs[-1] + expon.rvs()
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