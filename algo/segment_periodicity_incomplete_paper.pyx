import numpy as np
import matplotlib.pyplot as plt
import time, random

#Counts 1 in fixed period and notes down how many ones are there in each period.
#Implements the algorithm from Incomplete observation paper.

def plotInputSeries(series):
    plt.axis([0,1500,-1, 2])
    plt.plot(series)
    plt.show()

def getOneOffset(x, period):
    countOneOffsets = np.zeros(period)
    for index in range(len(x)):
        if (x[index] == 1):
            offset = index % period
            countOneOffsets[offset] = countOneOffsets[offset] + 1
            #print "index:", index, "offset:", offset
    norm = sum(x)
    normalize_x = countOneOffsets / norm
    #print "countOneOffsets:", countOneOffsets
    return normalize_x

def getZeroOffset(x, period):
    countOffsets = np.zeros(period)
    for index in range(len(x)):
        if (x[index] == 0):
            offset = index % period
            countOffsets[offset] = countOffsets[offset] + 1
            #print "index:", index, "offset:", offset
    norm = len(x) - sum(x)
    normalize_x = countOffsets / norm
    #print "countOneOffsets:", countOneOffsets
    return normalize_x

def getDiscrepancy(countZeroOffsets, countOneOffsets):
    Delta = countOneOffsets - countZeroOffsets
    return Delta

def getPeriodicityScore(Delta):
    score = 0
    for item in Delta:
        if item > 0:
            score = score + item
    return score

def getOneOffsetFromRawTimestamp(infilename):
    f = open(infilename)
    countOneOffsets = np.zeros(24)
    for line in f:
        values = line.split(":")
        timestamp = int(values[3])
        hour = int(time.strftime("%H", time.localtime(int(timestamp))))
        countOneOffsets[hour] = countOneOffsets[hour] + 1
    return countOneOffsets

def buildRandomPermuteSeq(N, numofOnes):
    randS = np.zeros(N)
    indexes = random.sample(range(0, N-1), numofOnes)
    randS[indexes] = 1
    return randS

def getBaseScore(s, L, U):
    numofOnes = int(sum(s))
    N = len(s)
    NUM_TRIAL = 2
    baseScores = {}
    for T in np.arange(L,U):
        s = buildRandomPermuteSeq(N, numofOnes)
        scores = np.array([])
        for i in range(NUM_TRIAL):
            muPos = getOneOffset(s,T)
            muMinus = getZeroOffset(s,T)
            Delta = getDiscrepancy(muPos, muMinus)
            score = getPeriodicityScore(Delta)
            scores = np.append(scores, score)
        baseScores[T] = np.mean(scores)
    #print "baseScores:", baseScores
    return baseScores

def plotRatioGraph(ax, s, T):
    muPos = getOneOffset(s,T)
    ax.bar(range(len(muPos)), muPos, align='center')
    ax.set_xlabel("Time offset (hours)")
    ax.set_title("Ratio graph")

def getPeriods(ax, s):
    L = 1
    U = 1000
    #plotRatioGraph(ax, s, 96)
    #return
    baseScores = getBaseScore(s, L, U)
    scores = np.array([])
    for T in np.arange(L,U):
        muPos = getOneOffset(s,T)
        muMinus = getZeroOffset(s,T)
        Delta = getDiscrepancy(muPos, muMinus)
        score = getPeriodicityScore(Delta) - baseScores[T]
        scores = np.append(scores, score)
        #print "T:", T, "score:", score
    if (ax != None):
        ax.plot(np.arange(len(scores)), scores)
        m = max(scores)
        am = np.argmax(scores)
        ax.plot([am],[m],'o',color='red')
        ax.set_ylabel("Periodicity Score")
        #ax.set_xlabel("Potential Time Period (Hours)")
        ax.set_title("Incomplete paper algorithm")
    #print "Max score:", max(scores), np.argwhere(abs(scores - np.amax(scores))<1e-10)
    return np.argwhere(abs(scores - np.amax(scores))<1e-10)[0][0]+1#np.argmax(scores)
