from gen import buildrelease
from algo import ParticleFilterBuildRelease as pfbrProcess
import math
from utils import util
import numpy as np
import matplotlib.pyplot as plt
import matplotlib, multiprocessing
from time import time
from multiprocessing import Pool
import math, scipy, os, random

from algo import fftperiod
from algo import histogram
from algo import autocorrelation
from algo import fftautocorrmixed
from algo import segment_periodicity_incomplete_paper as segment_perodicity
from _abcoll import Iterable
from algo.ParticleFilterBuildRelease import getax

matplotlib.rcParams['legend.numpoints'] = 1
matplotlib.rcParams['lines.linewidth'] = 3
matplotlib.rcParams['font.size'] = 20

def saveresults(filename, events, pfperiods, eventTypes):
    z = zip(events, pfperiods, eventTypes)
    np.save("results/"+filename, z)

def plotsignal(s):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(range(len(s)), s, color='black')
    ax.set_ylim([-0.2,1.2])
    ax.set_xlim([0,len(s)])
    ax.set_xlabel("Time")
    ax.set_ylabel("Event")
    #ax.set_title("Time window:"+str(EPOCH_DIV)+" ms")
    #plt.show()

def getNumberOfCores():
    return multiprocessing.cpu_count()*2

def generateSeries(T, std, rate, NT=40, distribution='normal'):
    series = buildrelease.StreamingSeries(period=T, phase=0.0, std=std, 
                              mode='phaseDrift', rate=rate, length = T*NT, distribution=distribution)
    return series.getEvents()[1:]

def setSeed():
    seed = (os.getpid()*int(time()*1e6 %100)) % 4294967295
    scipy.random.seed(seed)
    random.seed(seed)

def runPhaseDrift(T: float, std: float, rate: float,# runType, scheme):
                  NP, PR, runType, scheme):
    #return [rate, 1.0,5.0,10.0,20.0]
    setSeed()
    NPARTICLES = int(NP)
    PERIOD_RANGE = int(PR)
    print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    print T, std, rate, NPARTICLES, PERIOD_RANGE, runType, scheme
    events = generateSeries(T, std, rate)
    print "len(events):", len(events)
    """
    pf1 = pfbrProcess.ParticleFilter(NPARTICLES, [0.0, PERIOD_RANGE], 
                                     scheme, showgraph=False, verbose=True)
    pfperiod, pfperiods, eventTypes = pf1.run(events)
    saveresults(runType, events, pfperiods, eventTypes)
    if (runType == "numParticles"):
        print "NPARTICLES:", NPARTICLES, " pfperiod:", pfperiod 
        return np.array([NPARTICLES, pfperiod, 0, 0, 0])
    """
    seq = util.convertToTimeSeries(events, 10)
    fperiod = fftperiod.getSignallibPeriod(None, seq)/10.0
    #histperiod = histogram.getPeriods(None, events)
    autocorrperiod = autocorrelation.getPeriods(None, seq)/10.0
    seg_period = segment_perodicity.getPeriods(None, seq)/10.0
    
    print "FFT period:", fperiod, "AutoCorrelation period:", autocorrperiod, "seg_period:", seg_period
    #mWindowAvg = sum(pfperiods[-10:])/10
    #mWindowMedian = np.median(pfperiods[-10:])
    """
    print runType, T, std, rate, "pfperiod:", pfperiod, "fperiod:",fperiod, \
        "histperiod:", histperiod, "autocorrperiod:", autocorrperiod, "fftautocorrmixedperiod:", fftautocorrmixedperiod#, "seg_period:", seg_period
    if (runType == "noise"):
        return np.array([rate, pfperiod, fperiod, histperiod, autocorrperiod, fftautocorrmixedperiod])
    elif (runType == "std"):
        return np.array([std, pfperiod, fperiod, histperiod, autocorrperiod, fftautocorrmixedperiod])
    elif (runType == "period"):
        return np.array([T, pfperiod, fperiod, histperiod, autocorrperiod, fftautocorrmixedperiod])
    """
def runPhaseDriftOnline(T: float, std: float, rate: float, NP, PR, runType, scheme):
    setSeed()
    NPARTICLES = int(NP)
    PERIOD_RANGE = int(PR)
    print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    print T, std, rate, NPARTICLES, PERIOD_RANGE, runType, scheme
    events = generateSeries(T, std, rate)
    print "len(events):", len(events)
    pf1 = pfbrProcess.ParticleFilter(NPARTICLES, [0.0, PERIOD_RANGE], scheme, showgraph=False, verbose=False)
    pfperiod, pfperiods = pf1.run(events)
    if (runType == "numParticles"):
        print "NPARTICLES:", NPARTICLES, " pfperiod:", pfperiod 
        return np.array([NPARTICLES, pfperiod, 0, 0, 0])
    
    onlinefft, onlinehist, onlineautocorr = np.array([]),np.array([]),np.array([])
    for i in range(2, len(events)):
        contevents = events[:i]
        seq = util.convertToTimeSeries(contevents, 10)
        onlinefft = np.append(onlinefft, fftperiod.getSignallibPeriod(None, seq)/10.0)
        onlinesegment=np.append(onlinefft, segment_perodicity.getPeriods(None, seq)/10.0)
        onlineautocorr = np.append(onlineautocorr, autocorrelation.getPeriods(None, seq)/10.0)
    
    xyz = random.randint(0,100)
    print "pfperiods:", pfperiods, "onlinefft:", onlinefft
    np.save("results/online/pf"+str(xyz), pfperiods)
    np.save("results/online/fft"+str(xyz), onlinefft)
    np.save("results/online/segment"+str(xyz), onlinesegment)
    np.save("results/online/autocorr"+str(xyz), onlineautocorr)
    
    
def unwrap_f(args):
    return runPhaseDrift(*args)
def unwraponline_f(args):
    runPhaseDriftOnline(*args)

def experiment(runType, nParticles, prange, scheme):
    if (runType=="noise"):
        print "========NOISE EXPERIMENT========="
        #print "type(nParticles):", type(nParticles), "type(prange):", type(prange)
        noiseExperiment(runType, scheme, nParticles, prange)
    elif (runType=="std"):
        print "========STD EXPERIMENT========="
        stdExperiment(runType, scheme, nParticles, prange)
    elif (runType=="period"):
        print "========PERIOD EXPERIMENT========="
        periodExperiment(runType, scheme, nParticles, prange)
    elif (runType=="numParticles"):
        print "========NUMBER OF PARTCILES EXPERIMENT=========="
        numParticlesExperiment(runType, scheme, nParticles, prange)
    elif (runType=="network"):
        if (util.isDICE() == False):
            print "========Network dataset EXPERIMENT========="
            networkDataExperiment(runType, scheme, nParticles, prange)
        else:
            print "DICE is not configured to run"
    elif (runType=="pulserate"):
        print "pulserate experiment"
        pulserateExperiment(runType, scheme, nParticles, prange)
    elif (runType=="steprate"):
        print "steprate experiment"
        steprateExperiment(runType, scheme, nParticles, prange)
    elif(runType == "online"):
        print "========ONLINE EXPERIMENT========="
        onlineExperiment(runType, scheme, nParticles, prange)
    elif(runType=="neuron"):
        print "========NEURON EXPERIMENT========="
        neuronDataExperiment(runType, scheme, nParticles, prange)
    elif(runType == "cputime"):
        print "========CPU TIME EXPERIMENT========="
        cpuTimeExperiment(runType, scheme, nParticles, prange)
    elif(runType == "preiodvary"):
        print "========Period vary EXPERIMENT========="
        periodVaryExperiment(runType, scheme, nParticles, prange)
    elif(runType == "diffdistribution"):
        print "========Different inter-event distributions EXPERIMENT========="
        diffDistributionExperiment(runType, scheme, nParticles, prange)

def saveConfiguration(runType, stdRange, noiseRange, NPARTICLES, PERIOD_RANGE, scheme):
    configuration = {}
    configuration["runType"]=runType
    configuration['T']=T
    configuration['stdRange'] = stdRange
    configuration['noiseRange'] = noiseRange
    configuration['NITERATION'] = NITERATION
    configuration['PERIOD_RANGE'] = PERIOD_RANGE
    configuration['NPARTICLES'] = NPARTICLES
    import csv
    w = csv.writer(open("results/configuration_"+runType+"_"+scheme+"_"+\
                        str(NPARTICLES)+"_"+str(PERIOD_RANGE)+".csv", "w"))
    for key in configuration.keys():
        w.writerow([key, configuration[key]])

T = 10.0
std = 1.0
noise = 0.05
NITERATION = 16

def runDataSetOnline(filepath):
    print filepath
    events = np.loadtxt(filepath)
    #"""
    pf1 = pfbrProcess.ParticleFilter(256, [0.0, 100000], "Reweight_Noise|Systematic_Resample")
    pfperiod, pfPeriods, sumweights = pf1.run(events)
    np.save(filepath+"_PF1", pfPeriods)
    #return
    #"""
    onlinefft, onlineseg, onlineautocorr = np.array([]),np.array([]),np.array([])
    for i in range(2, len(events)):
        contevents = events[:i]
        seq = util.convertToTimeSeries(contevents, 10)
        fft_period = fftperiod.getSignallibPeriod(None, seq)/10.0
        onlinefft = np.append(onlinefft, fft_period)
        #onlineseg = np.append(onlineseg, segment_perodicity.getPeriods(None, seq)/1000.0)
        #onlineautocorr = np.append(onlineautocorr, autocorrelation.getPeriods(None, seq)/1000.0)
        print i, filepath, "online FFT:", fft_period
    #"""
    print "Saving:", (filepath+"_FFT"), (filepath+"_SEG"), (filepath+"_AutoCorr")
    np.save(filepath+"_FFT1", onlinefft)
    np.save(filepath+"_SEG1", onlineseg)
    np.save(filepath+"_AutoCorr1", onlineautocorr)
    #"""

def diffDistributionExperiment(runType, scheme, NPARTICLES, PERIOD_RANGE):
    print ">>>>>>>>>>>>>>>>>>>NORMAL DISTRIBUTION>>>>>>>>>>>>>>>>>>>>"
    print T, std, noise, NPARTICLES, PERIOD_RANGE, runType, scheme
    
    """
    events = generateSeries(T, std, noise, distribution='normal')
    print "len(events):", len(events)
    pf1 = pfbrProcess.ParticleFilter(NPARTICLES, [0.0, PERIOD_RANGE], scheme, showgraph=False, verbose=False)
    pfperiod, pfperiods = pf1.run(events)
    np.save("results/diffdistributions/normal", pfperiods)
    del pfperiods
    """
    print ">>>>>>>>>>>>>>>>>>>LAPLACE DISTRIBUTION>>>>>>>>>>>>>>>>>>>>"
    events = generateSeries(T, std, noise, distribution='laplace')
    print "len(events):", len(events)
    pf1 = pfbrProcess.ParticleFilter(NPARTICLES, [0.0, PERIOD_RANGE], scheme, showgraph=False, verbose=False)
    pfperiod, pfperiods = pf1.run(events)
    np.save("results/diffdistributions/laplace", pfperiods)
    del pfperiods
    """
    print ">>>>>>>>>>>>>>>>>>> Gamma DISTRIBUTION>>>>>>>>>>>>>>>>>>>>"
    events = generateSeries(T, std, noise, distribution='gamma')
    print "len(events):", len(events)
    pf1 = pfbrProcess.ParticleFilter(NPARTICLES, [0.0, PERIOD_RANGE], scheme, showgraph=False, verbose=False)
    pfperiod, pfperiods = pf1.run(events)
    np.save("results/diffdistributions/gamma", pfperiods)
    """
def steprateExperiment(runType, scheme, NPARTICLES, PERIOD_RANGE):
    pool = Pool(processes=6)
    filePaths = np.array([])
    baseFolder = "/home/abhirup/Documents/Research/dataset/steprate/population/signals/"
    for folder in os.listdir(baseFolder):
        #print "folder:"+folder
        for filename in os.listdir(baseFolder+folder):
            if ("std" in filename and "npy" not in filename and "png" not in filename):
                fullfilename = baseFolder+folder+"/"+filename
                #print "Processing:", fullfilename
                filePaths = np.append(filePaths, fullfilename)
    
    results = pool.map(func=runDataSetOnline, iterable=filePaths, chunksize=None)

def pulserateExperiment(runType, scheme, NPARTICLES, PERIOD_RANGE):
    pool = Pool(processes=6)
    filePaths = np.array([])
    baseFolder = "/home/abhirup/Documents/Research/dataset/pulserate/signals/"
    for folder in os.listdir(baseFolder):
        #print "folder:"+folder
        for filename in os.listdir(baseFolder+folder):
            if (filename.endswith("csv") == True):
                fullfilename = baseFolder+folder+"/"+filename
                #print "Processing:", fullfilename
                filePaths = np.append(filePaths, fullfilename)
    
    results = pool.map(func=runDataSetOnline, iterable=filePaths, chunksize=None)

def networkDataExperiment(runType, scheme, NPARTICLES, PERIOD_RANGE):
    from datatest import network_analyze
    events = network_analyze.getEvents()
    np.save("results/network_evnets", events)
    pf = pfbrProcess.ParticleFilter(NPARTICLES, [0.0, PERIOD_RANGE], scheme, showgraph=False, verbose=True)
    pfperiod, pfperiods = pf.run(events)
    np.save("results/network_pfperiods", pfperiods)
    histperiods = histogram.getOnlinePeriods(events)
    np.save("results/network_histperiods", histperiods)

def neuronDataExperiment(runType, scheme, NPARTICLES, PERIOD_RANGE):
    #from datatest import neuron_model
    events = np.load("datatest/neuron_data/events.npy")
    print "len(events):", len(events)
    seq = util.convertToTimeSeries(events)
    ax = (plt.figure()).add_subplot(111)
    print "FFT period:", fftperiod.getSignallibPeriod(ax, seq)
    plt.show()
    pf = pfbrProcess.ParticleFilter(NPARTICLES, [0.0, PERIOD_RANGE], scheme, showgraph=False, verbose=True)
    pfperiod, pfperiods = pf.run(events)
    np.save("results/neuron_pfperiods", pfperiods)
    histperiods = histogram.getOnlinePeriods(events)
    np.save("results/neuron_histperiods", histperiods)

def cpuTimeSingle(dummy_var):
    events = np.load("datatest/neuron_data/events.npy")
    events = np.repeat(events, 14)
    #print "len(events):", len(events)
    fileindex = random.randint(0,100)
    #pf = pfbrProcess.ParticleFilter(256, [0.0, 1e5], "Reweight_Noise_Systematic_Resample", showgraph=False, verbose=False)
    #pftimes = pf.run(events, isTimeExp=True)
    #np.save("results/cputime/pf"+str(fileindex), pftimes)
    ffttimes, autotimes, segtimes = np.array([]), np.array([]), np.array([])
    for i in range(10, len(events)):
        contevents = events[:i]
        seq = util.convertToTimeSeries(contevents, 10)
        #start = time()
        #fftperiod.getSignallibPeriod(None, seq)
        #x = time()-start
        #print "fft time:", x
        #ffttimes = np.append(ffttimes, x)
        #start = time()
        #autocorrelation.getPeriods(None, seq)
        #x = time()-start
        #print "autocorr time:", x
        #autotimes = np.append(autotimes, x)
        start = time()
        segment_perodicity.getPeriods(None, seq)
        x = time()-start
        print "seg time:", x
        segtimes = np.append(segtimes, x)
        #np.save("results/cputime/fft"+str(fileindex), ffttimes)
        #np.save("results/cputime/auto"+str(fileindex), autotimes)
        np.save("results/cputime/seg"+str(fileindex), segtimes)

def cpuTimeExperiment(runType, scheme, NPARTICLES, PERIOD_RANGE):
    pool = Pool(processes=5)
    pool.map(func=cpuTimeSingle, iterable=np.arange(10), chunksize=None)

def periodVaryExperiment(runType, scheme, NPARTICLES, PERIOD_RANGE):
    T1, T2, T3, T4, T5 = 4.0, 7.0, 8.0, 10.0, 2.0
    events1 = generateSeries(T1, 0.2*T1, 0.01/T1,  50)
    events2 = generateSeries(T2, 0.2*T2, 0.01/T2, 20)
    events3 = generateSeries(T3, 0.2*T3, 0.01/T3, 20)
    events4 = generateSeries(T4, 0.2*T4, 0.01/T4, 80)
    events5 = generateSeries(T5, 0.2*T5, 0.01/T5, 80)
    info = [[T1, len(events1)],
            [T2, len(events2)],
            [T3, len(events3)],
            [T4, len(events4)],
            [T5, len(events5)]]
    np.save("results/periodVary_info", info)
    events = np.append(events1, events2)
    events = np.append(events, events3)
    events = np.append(events, events4)
    events = np.append(events, events5)
    print "len(events1):", len(events1)
    print "len(events2):", len(events2)
    print "len(events3):", len(events3)
    print "len(events4):", len(events4)
    print "len(events5):", len(events5)
    print "len(events):", len(events)
    for i in range(3):
        print ">>>>>>>>>>>>> Start of Partcile filter <<<<<<<<<<<<<<<<<"
        pf1 = pfbrProcess.ParticleFilter(NPARTICLES, [0.0, PERIOD_RANGE], 
                                         scheme, showgraph=False, verbose=True)
        pfperiod, pfperiods = pf1.run(events)
        index = random.randint(0,100)
        np.save("results/periodVary_pfperiods"+str(index), pfperiods)

def numParticlesExperiment(runType, scheme, NPARTICLES, PERIOD_RANGE):
    powers = np.arange(3,11)
    xr = np.power(2, powers).astype(int)
    std = 6.0
    noise = 0.15
    saveConfiguration('numParticles', std, xr, NPARTICLES, PERIOD_RANGE, scheme)
    t1 = time()
    pool = Pool(processes=min(getNumberOfCores(), NITERATION*len(xr)))
    Ts = np.ones(NITERATION*len(xr))*float(T)
    Stds = np.ones(NITERATION*len(xr))*float(std)
    Rates = np.ones(NITERATION*len(xr))*float(noise)
    NPs = np.repeat(xr, NITERATION)
    PRs = np.ones(len(Ts))*PERIOD_RANGE
    z = zip(Ts, Stds, Rates, NPs, PRs, ['numParticles'] * len(Ts), [scheme] * len(Ts))
    results = pool.map(func=unwrap_f, iterable=z, chunksize=None)
    pool.close()
    pool.join()
    print "Time lapsed:", (time()-t1)
    np.save("results/numParticles_"+scheme+"_"+\
                        "_"+str(PERIOD_RANGE), results)
    xxx = np.load("results/numParticles_"+scheme+"_"+\
                        "_"+str(PERIOD_RANGE)+".npy")
    print xxx

def onlineExperiment(runType, scheme, NPARTICLES, PERIOD_RANGE):
    #runPhaseDriftOnline(T,std,noise,NPARTICLES,PERIOD_RANGE,None, scheme)
    #return
    noise=0.01
    std=2.0
    pool = Pool(processes=min(getNumberOfCores(), NITERATION))
    Ts = np.ones(NITERATION)*float(T)
    Stds = np.ones(NITERATION)*float(std)
    Rates = np.ones(NITERATION)*float(noise)
    NPs = np.ones(NITERATION)*NPARTICLES
    PRs = np.ones(NITERATION)*PERIOD_RANGE
    z = zip(Ts, Stds, Rates, NPs, PRs, [runType] * len(Ts), [scheme] * len(Ts))
    results = pool.map(func=unwraponline_f, iterable=z, chunksize=None)
    

def periodExperiment(runType, scheme, NPARTICLES, PERIOD_RANGE):
    #runPhaseDrift(10.0, 1.0, 0.05, NPARTICLES, PERIOD_RANGE, "noise", scheme)
    #return
    powerval = np.arange(1, 7)
    xr = np.power(2, powerval)
    saveConfiguration(runType, std, xr, NPARTICLES, PERIOD_RANGE, scheme)
    t1 = time()
    pool = Pool(processes=min(getNumberOfCores(), NITERATION*len(xr)))
    Ts = np.repeat(xr, NITERATION)
    Stds = np.ones(NITERATION*len(xr))*float(std)
    Rates = np.ones(NITERATION*len(xr))*float(noise)
    NPs = np.ones(len(Ts))*NPARTICLES
    PRs = np.ones(len(Ts))*PERIOD_RANGE
    z = zip(Ts, Stds, Rates, NPs, PRs, [runType] * len(Ts), [scheme] * len(Ts))
    results = pool.map(func=unwrap_f, iterable=z, chunksize=None)
    pool.close()
    pool.join()
    print "Time lapsed:", (time()-t1)
    np.save      ("results/"+runType+"_"+scheme+"_"+str(NPARTICLES)+"_"+str(PERIOD_RANGE), results)
    xxx = np.load("results/"+runType+"_"+scheme+"_"+str(NPARTICLES)+"_"+str(PERIOD_RANGE)+".npy")
    print xxx

def noiseExperiment(runType, scheme, NPARTICLES, PERIOD_RANGE):
    runPhaseDrift(10.0, 1.0, 0.05, NPARTICLES, PERIOD_RANGE, "noise", scheme)
    return
    xr = np.arange(0.0, 0.26, 0.05)
    saveConfiguration(runType, std, xr, NPARTICLES, PERIOD_RANGE, scheme)
    t1 = time()
    pool = Pool(processes=min(getNumberOfCores(), NITERATION*len(xr)))
    Ts = np.ones(NITERATION*len(xr))*T
    Vs = np.repeat(xr, NITERATION)
    NPs = np.ones(len(Ts))*NPARTICLES
    PRs = np.ones(len(Ts))*PERIOD_RANGE
    z = zip(Ts, np.ones(len(Vs)) * std, Vs, NPs, PRs, [runType] * len(Ts), [scheme] * len(Ts))
    #z = zip(Ts, np.ones(len(Vs)) * std, Vs, ['noise'] * len(Ts), [scheme] * len(Ts))
    #z = zip(Ts, np.ones(len(Vs)) * std, Vs, np.ones(len(Ts)) * NPARTICLES,
    #         np.ones(len(Ts)) * PERIOD_RANGE, ['noise'] * len(Ts), [scheme] * len(Ts))
    results = pool.map(func=unwrap_f, iterable=z, chunksize=None)
    pool.close()
    pool.join()
    print "Time lapsed:", (time()-t1)
    np.save      ("results/"+runType+"_"+scheme+"_"+str(NPARTICLES)+"_"+str(PERIOD_RANGE), results)
    xxx = np.load("results/"+runType+"_"+scheme+"_"+str(NPARTICLES)+"_"+str(PERIOD_RANGE)+".npy")
    print xxx

def stdExperiment(runType, scheme,NPARTICLES,PERIOD_RANGE):
    #runPhaseDrift(10, generateSeries(T=10, std=2, rate=0.1))
    #return
    xr = np.arange(0, 13, 2)
    saveConfiguration(runType, std, xr, NPARTICLES, PERIOD_RANGE, scheme)
    t1 = time()
    pool = Pool(processes=min(getNumberOfCores(), NITERATION*len(xr)))
    Ts = np.ones(NITERATION*len(xr))*T
    Stds = np.repeat(xr, NITERATION)
    NPs = np.ones(len(Ts))*NPARTICLES
    PRs = np.ones(len(Ts))*PERIOD_RANGE
    Noises = np.ones(len(Ts)) * noise
    #z = zip(Ts, Vs, np.ones(len(Vs))*noise, np.ones(len(Ts))*NPARTICLES,
    #         np.ones(len(Ts))*PERIOD_RANGE,['std']*len(Ts), [scheme]*len(Ts))
    z = zip(Ts, Stds, Noises, NPs, PRs, [runType] * len(Ts), [scheme] * len(Ts))
    #print z
    results = pool.map(func=unwrap_f, iterable=z, chunksize=None)
    pool.close()
    pool.join()
    print "Time lapsed:", (time()-t1)
    np.save      ("results/"+runType+"_"+scheme+"_"+str(NPARTICLES)+"_"+str(PERIOD_RANGE), results)
    xxx = np.load("results/"+runType+"_"+scheme+"_"+str(NPARTICLES)+"_"+str(PERIOD_RANGE)+".npy")
    print xxx
