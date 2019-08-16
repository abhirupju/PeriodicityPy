import scipy.stats as stats
import math
import numpy as np
import platform

Scheme_Plain="Plain"
Scheme_Resample_ESS="Resample_ESS"
Scheme_Dynamics_Mixed_PHalf="Dynamics_Mixed_PHalf"
Scheme_Dynamics_Mixed_PTwoSided="Dynamics_Mixed_PTwoSided"
Scheme_Dynamics_ESS="Dynamics_ESS"
Scheme_Dynamics_Inverse_ESS="Dynamics_Inverse_ESS"
Scheme_Dynamics_MW_W="Dynamics_MW_W"
Scheme_Reweight_Noise="Reweight_Noise"
Scheme_Reweight_Std="Reweight_Std"
Scheme_Systematic_Resample="Systematic_Resample"
Scheme_Overlapping_Expon_Ranges="Overlapping_Expon_Ranges"

def boundedSample(rv, mmin, mmax, tryCount=10):
    xx = rv.rvs()
    #tryCount = 100
    if(xx > mmin and xx < mmax):
        return xx
    if (tryCount < 0):
        return 0
    tryCount -= 1
    return boundedSample(rv, mmin, mmax, tryCount)

def convertToTimeSeries(events, resolution=1):
    n = int(math.ceil(max(events)*resolution)+1)
    s = np.zeros(n)
    for i in range(0, len(events)):
        e = events[i]
        x = int(math.floor(e*resolution))
        s[x] = 1
    return s

def isDICE():
    if (platform.linux_distribution()[0] == "Ubuntu"):
        return False
    return True