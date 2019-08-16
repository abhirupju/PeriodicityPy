import numpy as np
import matplotlib.pyplot as plt
import random, math
import scipy.stats as stats
from threading import Lock
import utils.util as util
from multiprocessing import Pool
import platform, random
from time import time
import pdist

MAX = 1e9 #float('-inf')
DEBUG = True
lock = Lock()

def logprint(*argv):
    lock.acquire()
    if (DEBUG == False):
        lock.release()
        return
    print "<",
    for arg in argv:
        print arg,
    print ">"
    lock.release()

###### plots ########
def getax():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #plt.ion()
    return ax

prevx, prevy=0, 0
def drawPeriodPosterior2(ax1, ax2, ax3, event, particles, weights):
    global prevx, prevy
    periods = [p.period for p in particles]
    stds = [p.std for p in particles]
    rates = [p.rate for p in particles]
    rgba_colors = np.zeros((len(particles),4))
    rgba_colors[:, 3] = weights
    #hist, bins = np.histogram(periods, bins=50, normed=True)
    #ax.bar(bins[:-1], hist, color='green', alpha=0.2)
    if (ax1 != None and len(particles)>1):
        ax1.set_xlim([0,400])
        ax1.set_ylim([0,30])
        ax1.scatter(np.ones(len(periods))*event, periods, alpha=0.2)#, color=rgba_colors)
        #estimate = np.mean(np.inner(periods, weights))
        estimate = np.median(periods)
        ax1.plot([event], estimate, 'o', color='r')
        #ax1.plot([event], np.median(periods), 'o', color='y')
        ax1.plot([prevx, event], [prevy, estimate], color='g')
        prevx, prevy = event, estimate
        ax1.set_title("periods")
    if (ax2 != None):
        ax2.set_xlim([0,400])
        ax2.scatter(np.ones(len(stds))*event, stds, alpha=0.2)#, color=rgba_colors)
        Stdestimate = np.median(stds)
        ax2.plot([event], Stdestimate, 'o', color='r')
        ax2.set_title("stds")
    if (ax3 != None):
        ax3.set_xlim([0,400])
        #ax3.set_ylim([0,1])
        ax3.scatter(np.ones(len(rates))*event, rates, alpha=0.2)#, color=rgba_colors)
        Noiseestimate = np.median(rates)
        ax3.plot([event], Noiseestimate, 'o', color='r')
        ax3.set_title("noise rates")
        ax3.set_xlabel("Event times")
        
    plt.draw()
    plt.pause(0.00001)

def boundedSample(rv, base, mmin, mmax, tryCount=10):
    xx = base + rv.rvs()
    #tryCount = 100
    if(xx > mmin and xx < mmax):
        return xx
    if (tryCount < 0):
        return base
    tryCount -= 1
    return boundedSample(rv, base, mmin, mmax, tryCount)

######## Particle ###########
def unwrap_self_f(arg, **kwarg):
    return Particle.getWeight(*arg, **kwarg)

class Particle:
    weight = 0
    pid = 0
    def __init__(self, period=0.0, std=0.0, rate=0.0,
                 lastPeriodic=0.0, lastPoint=0.0, lastNoiseEvent=0.0, 
                 periodbound=[], scheme="Plain"):
        self.scheme=scheme
        self.period = period
        self.std = std
        self.rate = rate
        self.lastPeriodic = lastPeriodic
        if (len(periodbound) == 0):
            periodbound = [0, MAX]
        self.periodbound = periodbound
        self.lastPoint = lastPoint #last event seen by this particle
        self.lastNoiseEvent = lastNoiseEvent
        #self.dynamics = stats.t(df=1, loc=0, scale=0.1)
        self.dynamics = pdist.t(df=1, loc=0, scale=0.1)
    def copyfrom(self, p):
        self.scheme = p.scheme
        self.period = p.period
        self.std = p.std
        self.rate = p.rate
        self.lastPeriodic = p.lastPeriodic
        self.lastPoint = p.lastPoint
        self.lastNoiseEvent = p.lastNoiseEvent
        self.periodbound = p.periodbound
        self.pi1 = p.pi1
        self.pi2 = p.pi2
        self.pi3 = p.pi3
        self.pid = p.pid
    def setPrior(self):
        if ((self.periodbound[1] - self.periodbound[0]) < 1e3):
            self.period = random.uniform(self.periodbound[0], self.periodbound[1])
        else:
            self.period = boundedSample(stats.expon(scale=1), 0, self.periodbound[0], self.periodbound[1])
        #self.std = stats.expon(scale=1).rvs() #boundedSample(rdv, 1e-10, self.periodbound[0], self.period)
        self.std = random.uniform(0, self.period)
        self.rate = stats.expon(scale=1/self.period).rvs()
        #random.seed(13214124)
        self.pi1 = random.uniform(0.8,1.0-1e-14)
        self.pi2 = random.uniform(0,1.0-self.pi1)
        self.pi3 = random.uniform(0,1.0-self.pi1-self.pi2)

    def applyDynamics(self, mw, w, ess, nParticles):
        x = stats.uniform.rvs()
        loc = self.period
        scale = 0.1
        if(util.Scheme_Dynamics_Mixed_PHalf in self.scheme):
            #print util.Scheme_Dynamics_Mixed_PHalf, "in", self.scheme
            if (x < self.pi1):
                pass
            else:
                loc = self.period/2
        if (util.Scheme_Dynamics_Mixed_PTwoSided in self.scheme):
            #print util.Scheme_Dynamics_Mixed_PTwoSided, "in", self.scheme
            if (x < self.pi1):
                pass
            elif (self.pi1 < x < self.pi2+self.pi1):
                loc = self.period/2
            else:
                loc = self.period*2
        if (util.Scheme_Dynamics_ESS in self.scheme):
            #print util.Scheme_Dynamics_ESS, "in", self.scheme
            scale = 1/ess
        if (util.Scheme_Dynamics_Inverse_ESS in self.scheme):
            scale = (ess/(nParticles*10.0))
        #self.period = util.boundedSample(stats.t(df=1, loc=loc, scale=scale), self.periodbound[0], self.periodbound[1])
        self.period = util.boundedSample(pdist.t(df=1, loc=loc, scale=scale), self.periodbound[0], self.periodbound[1])
        
        #self.period = boundedSample(self.dynamics, self.period, self.periodbound[0], self.periodbound[1])
        if (util.Scheme_Reweight_Std in self.scheme):
            self.std = boundedSample(self.dynamics, self.std, 0, MAX)
        else:
            self.std = boundedSample(self.dynamics, self.std, 0, self.period)
        self.rate = boundedSample(self.dynamics, self.rate, 0, MAX)
        
        #self.pi1 = boundedSample(stats.norm(loc=0, scale=0.01), self.pi1, 0, 1)
        #self.pi2 = boundedSample(stats.norm(loc=0, scale=0.01), self.pi2, 0, 1)
        #self.pi3 = boundedSample(stats.norm(loc=0, scale=0.01), self.pi3, 0, 1)
        #self.rate = util.boundedSample(stats.t(df=1,loc=self.rate,scale=0.1), 0, MAX)

    def getupPerupBg2(self, e):
        #r11 = stats.norm(loc=self.lastPeriodic+self.period, scale=self.std)
        r1 = pdist.norm(loc=self.lastPeriodic+self.period, scale=self.std)
        #q11 = stats.expon(loc=self.lastNoiseEvent, scale=1/self.rate)
        q1 = pdist.expon(loc=self.lastNoiseEvent, scale=1/self.rate)
        #q22 = stats.expon(loc=self.lastPoint, scale=1/self.rate)
        q2 = pdist.expon(loc=self.lastPoint, scale=1/self.rate)
        #print "e:", e
        upPer = (1-q2.cdf(e)) * r1.pdf(e)
        upBg = (1-r1.cdf(e)) * q1.pdf(e)
        return upPer, upBg

    def getupPerupBgStats(self, e):
        r1 = stats.norm(loc=self.lastPeriodic+self.period, scale=self.std)
        q1 = stats.expon(loc=self.lastNoiseEvent, scale=1/self.rate)
        q2 = stats.expon(loc=self.lastPoint, scale=1/self.rate)
        upPer = (1-q2.cdf(e)) * r1.pdf(e)
        upBg = (1-r1.cdf(e)) * q1.pdf(e)
        return upPer, upBg

    def getWeight(self, e):
        #prob = stats.norm(loc=self.lastPeriodic + self.period, scale=self.std).pdf(e)
        #self.lastPeriodic = e
        #return prob
        upPer, upBg = self.getupPerupBg2(e)
        #upPer1, upBg1 = self.getupPerupBgStats(e)
        pSum = upBg + upPer
        #pSum1 = upBg1 + upPer1
        #approxError = abs(pSum - pSum1)
        #if (approxError > 1e-7):
        #    print "ERROR! approxError:", approxError
        if (pSum == 0.0):
            #logprint ("*******************", e, self)
            pBg = 0.5
        else:
            pBg = upBg / pSum
        if (stats.uniform().rvs() < pBg):
            #print "noise"
            self.lastNoiseEvent = e
        else:
            #print "periodic"
            self.lastPeriodic = e
        self.lastPoint = e
        if (util.Scheme_Reweight_Noise in self.scheme):
            if (self.rate * self.period > 2):
                #pSum = pSum * (stats.expon(loc=2/self.period, scale=0.1).pdf(self.rate)/10)
                pSum = pSum * (pdist.expon(loc=2/self.period, scale=0.1).pdf(self.rate)/10)
        if (util.Scheme_Reweight_Std in self.scheme):
            if (self.std > self.period*0.95):
                pSum = pSum * (pdist.expon(loc=0.95*self.period, scale=1/self.period).pdf(self.std)/10)
        return pSum
    
    def __str__(self):
        return str(self.period)+","+str(self.std)+","+str(self.rate)+ \
            ","+str(self.lastPeriodic)+","+str(self.lastPoint)+","+str(self.lastNoiseEvent)

######## Weights ########
def getParticleWithMaxWeight(particles, weights):
    return particles[np.argmax(weights)]

def normalizeWeights(weights):
    s = sum(weights)
    if(s != 0):
        return weights/s
    else:
        print "sum of weights is 0"
        return np.repeat(1.0/float(len(weights)), len(weights))
    
######## Resample ############
def displayParticles(particles):
    for p in particles:
        logprint( p)

######## Basic Particle Filter ##########
class ParticleFilter:
    periodPlotter = None
    def __init__(self, nParticles, periodbound, scheme, showgraph=False, verbose=False):
        global DEBUG
        self.nParticles = nParticles
        if (len(periodbound) == 0):
            periodBound = [0, MAX]
        self.periodbound = periodbound
        self.scheme = scheme
        self.showgraph = showgraph
        DEBUG = verbose
        if (util.isDICE() == True):
            self.showgraph = False
            DEBUG = False
        self.particles = self.initParticles()
    
    def initParticles(self):
        p = np.array([])
        #rdp = stats.uniform(periodbound[0], periodbound[1])
        for i in range(self.nParticles):
            x = Particle(periodbound=self.periodbound, scheme=self.scheme)
            x.setPrior()
            x.pid = i
            p = np.append(p, x)
            #logprint( "Created particle:", x)
        """
        for lucky in range(10):
            p[lucky].period = 10
            p[lucky].std = 2
            p[lucky].rate = 0.1
            print "lucky particle pid:", p[lucky].pid, "index:", lucky
        #"""
        return p

    def getEss(self, weights):
        return 1/sum(np.power(weights,2))
    
    def needResample(self, weights):
        if (util.Scheme_Resample_ESS not in self.scheme):
            #print util.Scheme_Resample_ESS, "not in", self.scheme
            return True
        #print util.Scheme_Resample_ESS, "in", self.scheme
        ess = self.getEss(weights)
        #print "ess:", ess, "(self.nParticles/2):", (self.nParticles/2)
        if (ess < (self.nParticles/2)):
            return True
        return False
    def systematic_resample(self, weights):
        newparticles = np.array([])
        ess = self.getEss(weights)
        mw = max(weights)
        c = np.array([0])
        for i in range(self.nParticles):
            c = np.append(c, c[-1] + weights[i])
        u1 = random.uniform(0, 1/float(self.nParticles))
        i = 0
        for j in range(self.nParticles):
            u = u1 + (j/float(self.nParticles))
            while (u > c[i]):
                i += 1
            index = i-1
            p = self.particles[index]
            selected = Particle()
            selected.copyfrom(p)
            if (util.Scheme_Dynamics_MW_W in self.scheme):
                if (weights[index] < mw/2.0):
                    selected.applyDynamics(mw, weights[index], ess, self.nParticles)
            else:
                selected.applyDynamics(mw, weights[index], ess, self.nParticles)
            newparticles = np.append(newparticles, selected)
        return newparticles

    def resample(self, weights):
        newparticles = np.array([])
        N = len(self.particles)
        index = int(random.random() * N)
        beta = 0.0
        ess = self.getEss(weights)
        mw = max(weights)
        for i in range(N):
            beta += random.random() * 2.0 * mw
            while beta > weights[index]:
                beta -= weights[index]
                index = (index + 1) % N
            #logprint( "sampling index:", index)
            p = self.particles[index]
            #if (p.pid in range(10)):
            #    stillLucky = True
            selected = Particle()
            selected.copyfrom(p)
            if (util.Scheme_Dynamics_MW_W in self.scheme):
                #print util.Scheme_Dynamics_MW_W, "in", self.scheme
                if (weights[index] < mw/2.0):
                    selected.applyDynamics(mw, weights[index], ess, self.nParticles)
            else:
                selected.applyDynamics(mw, weights[index], ess, self.nParticles)
            newparticles = np.append(newparticles, selected)
        return newparticles
    def getEstimatedPeriod(self, weights):
        es = 0
        for i in range(self.nParticles):
            es += self.particles[i].period * weights[i]
        return es
    
    def run(self, events=[], isTimeExp=False):
        #displayParticles(particles)
        #if (self.showgraph == True):
        #    fig = plt.figure(1)
        #    ax1 = fig.add_subplot(311)
        #    ax1.plot(np.arange(0,400), np.ones(400)*10, 'r')
        #    ax2 = fig.add_subplot(312)
        #    ax3 = fig.add_subplot(313)
        #    ax1.grid(True)
        #    ax2.grid(True)
        #    ax3.grid(True)
        estimatedPeriod, estimatedPeriods = 0, np.array([])
        i = 0
        weights = np.ones(self.nParticles)
        times = np.array([])
        sumweights = np.array([])
        eventTypes = np.array([])
        for event in events:
            logprint(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            logprint(i, "event:", event)
            i += 1
            starttime = time()
            new_weights = np.array([])
            noiseEventVote = 0
            for p in self.particles:
                savedLastNoise = p.lastNoiseEvent
                new_weights = np.append(new_weights, p.getWeight(event))
                if (savedLastNoise != p.lastNoiseEvent):
                    #print "run():noisevote"
                    noiseEventVote = noiseEventVote + 1
            percentageVote = noiseEventVote/float(len(self.particles))
            #print "noiseEventVote:", noiseEventVote, "percentageVote:", percentageVote
            if (percentageVote > 0.5):
                print "      noise, percentageVote:", percentageVote
                eventTypes = np.append(eventTypes, 0)
            else:
                print "      periodic, percentageVote:", percentageVote
                eventTypes = np.append(eventTypes, 1)
            
            weights = np.array(new_weights) * weights
            ww = sum(weights)
            sumweights = np.append(sumweights, ww)
            normed_weights = normalizeWeights(weights)
            #if (i % 1 == 0 and self.showgraph == True):
            #    drawPeriodPosterior2(ax1, ax2, ax3, event, self.particles, normed_weights)
            estimatedPeriod = np.median(np.array([p.period for p in self.particles]))
            
            #esStd = np.median(np.array([p.std for p in self.particles]))
            #esNoise = np.median(np.array([p.rate for p in self.particles]))
            #print "estimatedPeriod:",estimatedPeriod
            estimatedPeriods = np.append(estimatedPeriods, estimatedPeriod)
            if (self.needResample(normed_weights) == True):
                logprint("Resample")
                if (util.Scheme_Systematic_Resample in self.scheme):
                    self.particles = self.systematic_resample(normed_weights)
                else:
                    self.particles = self.resample(normed_weights)
                weights = np.ones(self.nParticles)
            else:
                logprint("Skipped resample!")
                for j in range(self.nParticles):
                    p = self.particles[j]
                    p.applyDynamics(max(normed_weights), normed_weights[j], self.getEss(weights), self.nParticles)
            #times = np.append(times, time()-starttime)
            print event, "time:", time()-starttime, "estimatedPeriod:",estimatedPeriod, "sum(weights):", ww
        #if (self.showgraph==True):
        #    plt.ioff()
        if (isTimeExp == True):
            return times
        print("lengths", len(estimatedPeriods), len(sumweights))
        #self.drawPlot(estimatedPeriods, sumweights)
        self.drawPlot(events, estimatedPeriods)
        return estimatedPeriod, estimatedPeriods

    def drawPlot(self, estimatedPeriods, sumweights):
        f, axarr = plt.subplots(2, sharex=True)
        #seq = util.convertToTimeSeries(events)
        #axarr[0].bar(range(len(seq)), seq, color='black')
        #axarr[0].set_xlabel("Time (0.8 mSec)")
        #axarr[0].set_ylim([-0.2, 1.2])
        estimatedPeriods = estimatedPeriods[1:]
        axarr[0].plot(range(len(estimatedPeriods)), estimatedPeriods, color='black')
        axarr[0].set_ylabel("Estimated Period")
        
        sumweights = sumweights[1:]
        axarr[1].plot(range(len(sumweights)), sumweights, color='black')
        axarr[1].fill_between(range(len(sumweights)), 0, sumweights, facecolor='green', alpha=0.5)
        axarr[1].set_ylabel("$\sum w_i$")
        axarr[1].set_xlabel("# of samples")
