import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import math, timeit, os.path
import traceback

def buildnormcdftable():
    nd = stats.norm()
    table = np.array([])
    for x in np.arange(0, 4.0+0.01, 0.01):
        t = nd.cdf(x)
        table = np.append(table, t)
    np.save("normtable", table)

#table = None

class norm:
    def __init__(self, loc=0, scale=1):
        self.loc = loc
        self.scale = scale
        #print "self.loc:", self.loc, "loc:", loc
        self.var = pow(self.scale,2)
        self.denom = (1/math.sqrt(2*self.var*np.pi))
        if not os.path.exists("normtable.npy"):
            buildnormcdftable()
        self.table = np.load("normtable.npy")
    def pdf(self, x):
        return self.denom * math.exp(-1*pow((x-self.loc),2)/(2*self.var))

    def cdfUtility(self, x):
        x1 = (math.floor(x/0.01)) * 0.01
        x3 = (math.ceil(x/0.01)) * 0.01
        i1, i3 = int(x1/0.01), int(x3/0.01)
        if (i3 - i1 == 1):
            y1 = self.table[i1]
            y3 = self.table[i3]
            y = ((y3-y1)/(x3-x1))*(x - x1) + y1
        else:
            y = self.table[i1]
        return y

    def cdf(self, mx):
        #print "cdf:(",mx,") = ", type(mx), "self.loc:", self.loc, "self.scale:", self.scale
        x = (mx-self.loc)/self.scale
        if (-0.001 < x < 0.001):
            return 0.5
        elif (0 < x < 4):
            #print "+ cdf:(",x,") = ", type(x) 
            return self.cdfUtility(x)
        elif (-4 < x < 0):
            #print "- cdf:(",x,")", type(x)
            #traceback.print_stack()
            return 1 - self.cdfUtility(-1*x)
        elif ( x > 4.0):
            return 1.0
        else:
            return 0.0

class expon:
    def __init__(self, loc=0, scale=1):
        self.loc, self.scale = loc, scale
        self.l = 1/float(self.scale)
    def pdf(self, x):
        if (x < self.loc):
            return 0.
        return self.l * math.exp(-1*self.l * (x-self.loc))
    def cdf(self, x):
        if (x < self.loc):
            return 0.
        return 1-math.exp(-1*self.l* (x-self.loc))

class t:
    def __init__(self, df=1, loc=0, scale=1):
        self.loc, self.scale = loc, scale
    
    def icdf(self, p):
        return np.tan(np.pi * (p-0.5))

    def rvs(self, size=1):
        if (size > 1):
            return stats.t(df=1, loc=self.loc, scale=self.scale).rvs(size)
            #u = np.random.random(size=size)
            #return np.array([(self.icdf(x)*self.scale)+self.loc for x in u])
        else:
            return stats.t(df=1, loc=self.loc, scale=self.scale).rvs()
            #u = np.random.random()
            #return (self.icdf(u)*self.scale) + self.loc

def mynormpdf():
    norm(loc=0.1, scale=0.2).pdf(1)

def scinormpdf():
    stats.norm(loc=0.1, scale=0.2).pdf(1)

def mynormcdf():
    norm(loc=0.1, scale=0.2).cdf(1)

def scinormcdf():
    stats.norm(loc=0.1, scale=0.2).cdf(1)

def myexponpdf():
    expon(loc=0.1, scale=0.2).pdf(1)

def sciexponpdf():
    stats.expon(loc=0.1, scale=0.2).pdf(1)

def myexponcdf():
    expon(loc=0.1, scale=0.2).cdf(1)

def sciexponcdf():
    stats.expon(loc=0.1, scale=0.2).cdf(1)

def myTLocTime():
    x = t().rvs()

def scipyTlocTime():
    x = stats.t(df=1).rvs()

if __name__ == "__main__":
    print "normpdf", timeit.timeit(mynormpdf, number=256), timeit.timeit(scinormpdf, number=256)
    print "normcdf", timeit.timeit(mynormcdf, number=256), timeit.timeit(scinormcdf, number=256)
    print "exponpdf", timeit.timeit(myexponpdf, number=256), timeit.timeit(sciexponpdf, number=256)
    print "exponcdf", timeit.timeit(myexponcdf, number=256), timeit.timeit(sciexponcdf, number=256)
    print "Tdraw", timeit.timeit(myTLocTime, number=256*3), timeit.timeit(scipyTlocTime, number=256*3)