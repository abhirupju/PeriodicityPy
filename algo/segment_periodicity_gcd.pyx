import numpy as np
import matplotlib.pyplot as plt
import time, random, math, decimal
from sklearn.mixture import GMM
from scipy.signal import argrelextrema
import scipy.fftpack
from algo import fftperiod as fftp
from algo import autocorrelation as acorr

from scipy import stats
from scipy import signal
from time import time
from time import sleep

L = 6
U = 17*5*3
NUM_TRIAL = 10
randseqs = np.array([])

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def getax():
    fig = plt.figure()
    return fig.add_subplot(111)

def fitGMM(s, T, draw=False):
    X = np.array([])
    for i in s:
        X = np.append(X, i%T)
    X = X.reshape(len(X), 1)
    #------------------------------------------------------------
    # Learn the best-fit GMM models
    #  Here we'll use GMM in the standard way: the fit() method
    #  uses an Expectation-Maximization approach to find the best
    #  mixture of Gaussians for the data
    # fit models with 1-10 components
    N = np.arange(1, 21)
    models = [None for i in range(len(N))]
    
    for i in range(len(N)):
        print "N[i]", N[i]
        models[i] = GMM(N[i]).fit(X)
    
    # compute the AIC and the BIC
    AIC = [m.aic(X) for m in models]
    BIC = [m.bic(X) for m in models]
    
    #------------------------------------------------------------
    # Plot the results
    #  We'll use three panels:
    #   1) data + best-fit mixture
    #   2) AIC and BIC vs number of components
    
    #M_best = models[np.argmin(AIC)]
    M_best = models[12]
    x = np.linspace(0, T-1, T)
    logprob, responsibilities = M_best.score_samples(x.reshape(len(x),1))
    pdf = np.exp(logprob)
    pdf_individual = responsibilities * pdf[:, np.newaxis]
    print M_best
    maximas = argrelextrema(pdf, np.greater)[0]
    minimas = argrelextrema(pdf, np.less)[0]
    c = np.empty((maximas.size + minimas.size,), dtype=maximas.dtype)
    c[0::2] = maximas
    c[1::2] = minimas
    c = np.insert(c, 0, 0)
    c = np.append(c, len(pdf)-1)
    heights = np.array([])
    for i in range(len(maximas)):
        j = np.where(c==maximas[i])[0][0]
        heights = np.append(heights, pdf[c[j]] - min(pdf[c[j-1]], pdf[c[j+1]]))
    #Config value: What is a significant height of peaks?
    counter = np.where(heights > max(heights)/5.0)[0].size
    if (draw==True):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(X, T, normed=True, histtype='stepfilled', alpha=0.4)
        ax.plot(x, pdf, '-k')
        ax.plot(x, pdf_individual, '--k', color='red')
        ax.text(0.04, 0.96, "Best-fit Mixture", ha='left', va='top', transform=ax.transAxes)
        ax.set_xlabel('$x$')
        ax.set_ylabel('$p(x)$')
        """
        # plot 2: AIC and BIC
        ax = fig.add_subplot(212)
        ax.plot(N, AIC, '-k', label='AIC')
        ax.plot(N, BIC, '--k', label='BIC')
        ax.set_xlabel('n. components')
        ax.set_ylabel('information criterion')
        ax.legend(loc=2)
        """
        plt.show()
    return counter

def getPeakCountMassCompare(dist):
    #print dist
    for i in range(3, 20):
        dividers = np.linspace(0, len(dist), i).astype(int)
        #print dividers
        ds = np.array([])
        for index in range(1, len(dividers)):
            ds = np.append(ds, sum(dist[dividers[index-1]:dividers[index]]))
        m = np.mean(ds)
        diff = np.fabs(ds - m)
        if (np.where(diff<1e-2)[0].size == diff.size):
            print "possible peakcount:", (i-1), ds
    return 1

def getPeakCountFFTSmoothing(dist):
    x = np.arange(len(dist))
    rft = np.fft.rfft(dist)
    rft[5:] = 0   # Note, rft.shape = 21
    y_smooth = np.fft.irfft(rft)
    if (y_smooth.size != x.size):
        y_smooth = np.append(y_smooth, 0)
    #print x.shape, dist.shape, y_smooth.shape
    plt.bar(x, dist, align='center', label='Original')
    plt.plot(x, y_smooth, label='Smoothed')
    plt.legend(loc=0).draggable()
    #plt.show()

def getPeakCountStd(dist):
    m = np.mean(dist)
    std = np.std(dist)
    x = np.where(dist > (m+(std)))
    diff = np.diff(x)
    peakcount = np.where(diff>3)[0].size + 1
    print diff, peakcount
    return peakcount

def getPeakCountMixed(dist):
    return 0

def alignDist(distribution, orient="mean"):
    if (orient == "mean"):
        a = int(mean(distribution))
    elif (orient == "max"):
        a = np.argmax(distribution)
    else:
        return distribution
    if (a > len(distribution)/2):
        rotateIndex = (len(distribution) - a + len(distribution)/2)
    else:
        rotateIndex = len(distribution)/2 - a
    return np.roll(distribution, rotateIndex)

def getDistributionScore(dist):
        val=0
        m = np.mean(dist)
        for d in dist:
            v = math.pow(d-m, 2)
            val += v
        return math.sqrt(val)

def isUnifomDistribution(dist, n):
    bmean, bstd = getBaseScorePeriod(len(dist))
    v = getDistributionScore(dist)
    #Config value: What is a significant value to be non-uniform distribution?
    if (v > bmean+bstd):
        return True, v
    else:
        return False, v

def mean(distribution):
    mean = 0
    for i in range(len(distribution)):
        mean = mean + (i*distribution[i])
    return mean

def std(distribution):
    std = 0
    m = mean(distribution)
    for i in range(len(distribution)):
        std = std + (math.pow((i-m),2) * distribution[i])
    return math.sqrt(std)

def getDistribution(s, T):
    if (T <= 1):
        return np.array([])
    T = int(T)
    dist = np.zeros(T)
    for item in s:
        offset = int(item % T)
        #changed_offset = convertToThinnerDistribution(offset, T)
        dist[offset] = dist[offset] + 1
    norm = len(s)
    nDist = alignDist(dist/float(norm))
    return nDist

def getGcdFFT(dist):
    #g = fftp.getPeriods(None, dist)
    g = fftp.getSignallibPeriod(None, dist)
    #G = fftp.plotBuiltinPeriodogram(getax(), dist)
    return g

def getGcdAutoCorr(ax, dist):
    g = acorr.getPeriods(ax, dist, len(dist)/2+1)
    return g

def getgcd(dist):
    gcd = getGcdFFT(dist)
    return gcd

def confirmgcd(dist, g):
    return True
    if (g != len(dist)):
        g1 = getGcdAutoCorr(None, dist)
        if (g != g1):
            print "FFT g:", g, " autocorr g:", g1 
            return False
    else:
        print "Cannot confirm g:", g, " len(dist):", len(dist)
    return True

def analyzePeriod(s, T, draw=False):
    nDist = getDistribution(s, T)
    if (len(nDist) == 0):
        return 1
    #g = T/getPeakCountMassCompare(nDist)
    #g = T/fitGMM(s, T, draw)
    #g = T/getPeakCountStd(nDist)
    g = getgcd(nDist)
    return g

def annotate(ax, mx, my, annotationText, dx=10, dy=-0.5):
    ax.plot([mx], [my], 'o', color='red')
    ax.annotate(annotationText, xy=(mx, my), xytext=(mx+dx, my+dy),arrowprops=dict(arrowstyle="->"), fontsize=20)

def plotDistribution(ax, nDist):
    odist = alignDist(nDist, "max")
    ax.bar(range(len(odist)), odist, align='center', color='black')
    ax.set_xlabel("Time")
    ax.set_ylabel("Probability of event")
    ax.set_xlim([0,len(nDist)])
    return ax

primes = np.array([3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997])

def getBaseScorePeriod(T):
    global randseqs
    scores = np.array([])
    for i in range(NUM_TRIAL):
        randDist = getDistribution(randseqs[i], T)
        scores = np.append(scores, getDistributionScore(randDist))
    return np.mean(scores), np.std(scores)

def getBaseScore(x):
    baseMeanScores, baseStdScores = {}, {}
    for T in np.arange(L,U):
        baseMeanScores[T], baseStdScores[T] = getBaseScorePeriod(T)
    return baseMeanScores, baseStdScores

def createRandomSeqs(n):
    global randseqs
    randseqs = getRandomSeries(n)
    for i in range(1, NUM_TRIAL):
        x = getRandomSeries(n)
        randseqs = np.vstack((randseqs, x))

#### Recursive Multiplicative progressive algorithm #############
def factorize(n):
    factors = np.array([1, n])
    for i in range(2, int(np.sqrt(n))+1):
        if ((n%i) == 0):
            factors = np.append(factors, i)
            factors = np.append(factors, n/i)
    factors = np.sort(np.unique(factors))
    return factors

def checkNpeaksValidity(f, npeaks, base, basenpeaks, prefix = ""):
    x = int(f/base)
    factors = factorize(x)*basenpeaks
    #print prefix, f, npeaks, base, basenpeaks, " factors:", factors, "len(np.where(factors == npeaks)[0]):", len(np.where(factors == npeaks)[0])
    if (len(np.where(factors == npeaks)[0]) == 1):
        return True
    return False

def updateSupportnCount(f, npeaks, base, basenpeaks, support, nCount, MAXT, prefix):
    gcd = f/npeaks
    if (checkNpeaksValidity(f, npeaks, base, basenpeaks, prefix) == True and gcd > 1):
        support += 1.0 / int(float(MAXT)/float(base))
        nCount = MAXNBADPERIOD
    else:
        if (nCount > 0):
            nCount -= 1
    return support, nCount

potentialPeriods, seen = np.array([[0,0],[1,0]]), np.array([[0,0],[1,0]])
threshold, MAXNBADPERIOD = 0, 2

def recursivePeriodDetect(seq, base, incr, basenpeaks, MAXT, prefix=""):
    global seen, potentialPeriods, threshold
    support, nCount, isBad = 0, MAXNBADPERIOD, False
    prefix = prefix+"\t"
    for f in range(base, MAXT, incr):
        xx = np.where(seen[:,0]==f)[0]
        if (len(xx) > 0):
            peaks = seen[xx[0]][1]
            support, nCount = updateSupportnCount(f, peaks, base, basenpeaks, support, nCount, MAXT, prefix)
            #print prefix, f, "peaks", peaks, "Continue.", "nCount:", nCount
            if (nCount <= 0 and incr != 1):
                #print prefix, base, incr, basenpeaks, MAXT, support, " Bad period RETURN.."
                return
            continue
        dist = getDistribution(seq, f)
        score = getDistributionScore(dist)
        g = getgcd(dist)
        gcd = int(round(g))
        #print prefix, "f:", f, "g:", g, "gcd:", gcd, "p1:", (len(dist)/g)
        npeaks = len(dist)/gcd
        seen = np.vstack((seen, np.array([f, npeaks])))
        support, nCount = updateSupportnCount(f, npeaks, base, basenpeaks, support, nCount, MAXT, prefix)
        #print prefix, f, npeaks, support, "nCount:", nCount
        if (nCount <= 0 and incr != 1):
            #print prefix, base, incr, basenpeaks, MAXT, support, " Bad period RETURN.."
            potentialPeriods = np.vstack((potentialPeriods, np.array([base, support])))
            return
        if (score > threshold and gcd > 1 and gcd != base):
            threshold = score
            #print prefix, "gcd:", gcd, "Recursive call start >>"
            recursivePeriodDetect(seq, gcd, gcd, npeaks, MAXT, prefix)
    #print prefix, base, incr, basenpeaks, MAXT, support, " period investigation done RETURN.."
    potentialPeriods = np.vstack((potentialPeriods, np.array([base, support])))

def getProbablePeriods(potentialPeriods):
    bestsupport, bestperiod = 0, 0
    for i in potentialPeriods:
        #print i
        if (i[1] > bestsupport):
            bestsupport = i[1]
            bestperiod = i[0]
            if (bestsupport == 1.0):
                break
    return bestperiod

#### Multiresolution algorithm #############
def getLowResolutionDistribution(dist, r):
    lows = np.array([])
    for index in range(0,len(dist),r):
        #print index, index+r
        lows = np.append(lows, sum(dist[index:index+r]))
    return lows

def getRandomSeries(n):
    return np.random.uniform(0, n, size=n)

def createLowerResolutionSeries(s, r, n):
    #Converts event series
    ls = np.array([])
    for index in range(min(n, len(s))):
        i = s[index]
        x = int(i / r)
        ls = np.append(ls, x)
    return np.unique(ls.astype(int))

def multiResolutionAlgorithm(seq):
    invdensity = max(seq)/float(len(seq))
    h = int(np.log2(invdensity))
    x = np.unique(np.linspace(int(max(1,h-4)), h).astype(int))
    ps = np.array([])
    t = time()
    for r in x[::-1]:
        res = math.pow(2,r)
        #n = int((2000*(h-1))/(r))
        ls = createLowerResolutionSeries(seq, res, len(seq))
        n = max(ls)
        print "n:", n
        #fig, ax = plt.subplots()
        #ax.set_color_cycle(['black', 'black'])
        period = getFirstDivisor(None, ls, res, 10, n/4)
        #plt.show()
        ps = np.append(ps, period)
    period = np.median(ps)
    print ">>>>Time(sec):", time()-t, " period estimate:", period
    return period
    #Refine the process with highest resolution
    t = time()
    e = 0.08
    lb, hb = period/float(1+e), period/float(1-e)
    finalperiod = getFirstDivisor(None, seq, 1, lb, hb)
    if (finalperiod < lb or finalperiod > hb):
        print "FFT is bad! It found a bad peak!!!!!"
        return period
    print "Time(sec):", time()-t, "final period:", finalperiod
    return finalperiod

def peakDynamicThreshold(data):
    """
    http://ce-publications.et.tudelft.nl/publications/877_developing_and_implementing_peak_detection_for_realtime_ima.pdf
    """
    m = np.mean(data)
    std = np.std(data)
    deviation = sum(np.fabs(data-m))/len(data)
    k = 0.5
    threshold = (max(data)+m)/2 + k*deviation
    return threshold

def getFirstDivisor(ax, seq, res, lb, hb):
    #print "Investigating <",lb,",",hb,"> resolution:", res
    ps, pcs, runningThresh, potentialTs = np.array([]), np.array([]), np.array([]), np.array([])
    x = np.arange(lb, hb)
    for f in x:
        dist = getDistribution(seq, f)
        score = getDistributionScore(dist)
        ps = np.append(ps, score)
        g = analyzePeriod(seq, f, draw=False)
        pcs = np.append(pcs, f/g)
        threshold = peakDynamicThreshold(ps)
        print "f:", f, "threshold:", threshold
        runningThresh = np.append(runningThresh, threshold)
        if (score > threshold):
            print f, "score:", score, threshold
            potentialTs = np.append(potentialTs, f)
    amax = x[np.argmax(ps)]
    g = analyzePeriod(seq, amax, draw=False)
    npeak = pcs[np.argmax(ps)]
    print "potentialTs:", potentialTs
    print "amax:", amax, "#peak:", npeak, "g:", g,
    period = res*g
    maximas = argrelextrema(ps, np.greater)[0]
    #print "maximas:", maximas
    threshold = peakDynamicThreshold(ps)
    if (ax != None):
        ax.plot(x, ps, label='resolution:'+str(res))
        ax.plot(x, np.ones(len(x))*threshold, color='red')
        ax.plot(x, runningThresh, color='green')
        for mx in maximas:
            my = ps[mx]
            #print ">>>mx:", mx, "my:", my, "threshold:", threshold
            if (my > threshold):
                significance = (my-threshold)/my
                potentialPeriod = (mx/pcs[mx])*res
                #print significance, "mx(amax):", mx, "g:", "peaks:",pcs[mx], "mx/pcs[mx]:",mx/pcs[mx],"potentialPeriod:", potentialPeriod
                ax.annotate(str(int(pcs[mx])), xy=(mx+lb, my))
        ax.legend(loc='best')
        ax.set_xlabel("Potential Period")
        ax.set_ylabel("Periodicity Score")
    #print "\t\t","res:",res,"period:", period #, "amax:", amax, "g:", g
    return period

def experimentMultiResolution(seq):
    invdensity = max(seq)/float(len(seq))
    h = int(np.log2(invdensity))
    res = 8#math.pow(2, h-1)
    print "res:", res,
    ls = createLowerResolutionSeries(seq, res, len(seq))
    period = getFirstDivisor(getax(), ls, res, 2, 200)
    #factorize(period)
    dist = getDistribution(seq, period)
    #plotDistribution(getax(), dist)
    #fftp.plotBuiltinPeriodogram(getax(),dist)
    return period

#### Brute Force algorithm #############
def exhaustiveAlgorithm(ax, seq):
    L = 2
    ps = np.array([])
    U = max(10, int(max(seq)/2))
    print "exhaustiveAlgorithm L:",L,"U:",U
    x = np.arange(L, U)
    prevscore = 0
    for f in x:
        dist = getDistribution(seq, f)
        val = getDistributionScore(dist)
        ps = np.append(ps, val)
        #if (prevscore == 0):
        #    prevscore = val
        #g = getgcd(dist)
        #ax1.plot([f-1, f],[prevscore, val],color='black')
        #ax1.annotate(str(int(f/g)), xy=(f+1,val))
        #plt.draw()
        #prevscore = val
    if (len(ps) <= 0):
        return 0
    mx = x[np.argmax(ps)]
    my = ps[mx-2]
    if (ax != None):
        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        ax.plot(x, ps, color='black')
        #ax.plot([mx], [my], 'o', color='black')
        #ax.annotate("24 hours", xy=(mx, my), xytext=(mx+20, my-0.01),arrowprops=dict(arrowstyle="->"), fontsize=25)
        #ax.set_xticks(np.arange(0,200, 24))
        ax.set_xlabel("Potential Period")
        ax.set_ylabel("Periodicity Score")
        ax.set_title("GCD algorithm T0="+str(mx))
    gcd = getgcd(getDistribution(seq, mx))
    npeaks = mx/gcd
    print "period:", mx
    return mx, npeaks

#### Fractional resolution algorithm #############
def createLowerResolutionSeriesFractional(seq, res):
    #Handles fractional resolution
    ls = np.array([])
    for e in seq:
        x = e/res
        xf = int(math.floor(x))
        ls = np.append(ls, int(math.floor(x)))
        xc = int(math.ceil(x))
        d = (res*xc) - e
        #print "e:", e, "x:", x, "xf:", xf, "xc:", xc, "d:", d,
        if (0< d < 1):
            #print "xc added",
            ls = np.append(ls, xc)
        #print "xf added"
    return np.unique(ls.astype(int))

def recurssiveAlgorithmAdaptation(seq):
    global seen, potentialPeriods, threshold
    MAXT = int(max(seq)/2)
    recursivePeriodDetect(seq, 2, 1, 0, MAXT)
    period = getProbablePeriods(potentialPeriods)
    if (period == 0):
        return 0, 0
    xx = np.where(seen[:,0]==period)[0]
    peaks = seen[xx[0]][1]
    seen, potentialPeriods, threshold = np.array([[0,0],[1,0]]), np.array([[0,0],[1,0]]), 0
    return period, peaks

def experimentFractionalResolution(seq):
    invdensity = max(seq)/float(len(seq))
    h = int(np.log2(invdensity))
    x = np.unique(np.linspace(int(max(1,h-4)), h).astype(int))
    ps = np.array([])
    for r in x[::-1]:
        res = r + np.random.uniform(0,1)
        ls = createLowerResolutionSeriesFractional(seq, res)
        period1, npeaks1 = exhaustiveAlgorithm(None, ls)
        period, npeaks = recurssiveAlgorithmAdaptation(ls)
        if (npeaks != 0):
            trueperiodextimate = (period/npeaks) * res
        else:
            trueperiodextimate = None
        trueperiodextimate1 = (period1/npeaks1) * res
        print "Recursive algorithm: res:", res, "period:", period, "peaks:", npeaks, "period estimate:", trueperiodextimate
        print "Exhaustive algorithm: res:", res, "period:", period1, "peaks:", npeaks1, "period estimate:", trueperiodextimate1

def getPeriods(ax, seq, prefix=""):
    #experimentFractionalResolution(seq)
    """
    dist = getDistribution(seq, 105)
    fig = plt.figure()
    #ax = fig.add_subplot(513)
    plotDistribution(ax, dist)
    ax = fig.add_subplot(515)
    print "getgcd(dist):", getgcd(dist)
    print "New FFT algo perod:", fftp.getPeriods2(ax, dist)
    plt.show()
    exit(0)
    """
    #fig = plt.figure()
    #ax = fig.add_subplot(211)
    period, npeaks = exhaustiveAlgorithm(ax, seq)
    ax = getax()
    plotDistribution(ax, getDistribution(seq, period))
    #period = iteratePeriod(seq)
    #MAXT = int(max(seq)/20)
    #recursivePeriodDetect(seq, 2, 1, 0, MAXT)
    #period = getProbablePeriods(potentialPeriods)
    
    #period = multiResolutionAlgorithm(seq)
    #ax = fig.add_subplot(212)
    #plotDistribution(ax, getDistribution(seq, period))
    #plt.savefig('plot/'+prefix+'.png', bbox_inches='tight')
    #plt.close(fig)
    #getFirstDivisor(seq)
    #period = iteratePeriod(seq)
    #period = multiResolutionAlgorithm(seq)
    #print "gcd algorithm Period:", period
    return period
