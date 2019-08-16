#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import get_window
import os, sys, math, random
from optparse import OptionParser
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from time import time

def annotate(ax, mx, my, annotationText, dx=10, dy=-0.5):
    ax.plot([mx], [my], 'o', color='red')
    ax.annotate(annotationText, xy=(mx, my), xytext=(mx+dx, my+dy),arrowprops=dict(arrowstyle="->"), fontsize=20)

def getkneepointarr(r):
    kr, beforeknee = np.array([0]), True
    for i in range(1, len(r)):
        if (beforeknee == True):
            d = r[i-1] - r[i]
            if ( d > 0):
                #print "before knee:",i
                kr = np.append(kr, 0)
            else:
                #print "knee end:", i
                kr = np.append(kr, r[i])
                beforeknee = False
        else:
            #print "after knee:", i
            kr = np.append(kr, r[i])
    return kr

def calculateAutoCorrelation(ax, y):
    r = autocorr(y)
    dr = getkneepointarr(r)
    mx = np.argmax(dr)
    if (ax != None):
        ax.set_ylabel("Autocorrelation value")
        ax.set_xlabel("Potential Period")
        ax.plot(range(1,len(r)+1), r, color='black')
        ax.plot(range(1,len(dr)+1), dr, color='green')
        ax.set_title("T0:"+str(mx))
        my = max(dr)
        annotate(ax, mx, my, "Period:"+str(mx), 5, 0.05)
    #plt.show()
    return mx

def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size/2:]

def getPeriods(ax, x):
    return calculateAutoCorrelation(ax, x)


