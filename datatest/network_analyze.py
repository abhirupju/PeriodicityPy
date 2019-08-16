#!/usr/bin/python
from pcapfile import savefile
import numpy as np
import matplotlib.pyplot as plt
import sys,math,os
from itertools import groupby

basetimestamp = 0
lastfiletimestamp = 0
lastcount = 0

def convertToTimeSeq(events):
    n = int(max(events)+1)
    s = np.zeros(n)
    print "len(events):", len(events)
    for e in events:
        s[e] = 1
    return s

#consider all files joined
def parsepcapfilejoin(s, filename):
    global basetimestamp
    global lastfiletimestamp
    global lastcount
    testcap = open(filename, 'rb')
    capfile = savefile.load_savefile(testcap, verbose=True)
    count = 1
    prevtimestamp = capfile.packets[0].timestamp
    if (basetimestamp == 0):
        basetimestamp = prevtimestamp
    else:
        if (lastfiletimestamp == prevtimestamp):
            count = lastcount
    print "process:", filename, "len(capfile.packets):", len(capfile.packets)
    print basetimestamp
    ts = np.array([])
    for pkt in capfile.packets[1:]:
        timestamp = pkt.timestamp
        d = timestamp - prevtimestamp
        if (d == 0):
            count += 1
        else:
            #print prevtimestamp, "count:", count
            if (count > THRESH_PKT):
                s = np.append(s, prevtimestamp-basetimestamp)
            prevtimestamp = timestamp
            count = 1
    #print prevtimestamp, "count:", count
    if (count > THRESH_PKT):
        s = np.append(s, prevtimestamp-basetimestamp)
    lastfiletimestamp = timestamp
    lastcount = count
    s = np.unique(s)
    return s

def plotsignal(s):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(range(len(s)), s, color='black')
    ax.set_ylim([-0.2,1.2])
    ax.set_xlim([0,len(s)])
    ax.set_xlabel("Time ("+str(RESOLUTION)+" ms)", fontsize=25)
    ax.set_ylabel("Event >"+str(THRESH_PKT)+"in "+str(RESOLUTION)+" ms", fontsize=25)
    #ax.set_title("Time window:"+str(EPOCH_DIV)+" ms")
    plt.show()

RESOLUTION = 15 #ms
THRESH_PKT = 500 #univ:30, caida:3
RES_MS = 1
EPOCH_DIV = RESOLUTION*30 #3 sec #ms univ:5e3 caida:60e3

def getTimeStamp(pkt):
    return int((pkt.timestamp*1000 + int(pkt.timestamp_ms/(1000*RES_MS)))/RESOLUTION)

def parsepcapfilesingle(s, cs, filename):
    testcap = open(filename, 'rb')
    capfile = savefile.load_savefile(testcap, verbose=True)
    count = 1
    prevtimestamp = getTimeStamp(capfile.packets[0])
    basetimestamp = prevtimestamp
    print "filename:", filename, "len(capfile.packets):", len(capfile.packets)
    print basetimestamp
    ts = np.array([])
    for pkt in capfile.packets[1:]:
        timestamp = getTimeStamp(pkt)
        d = timestamp - prevtimestamp
        if (d == 0):
            count += 1
        else:
            #print "count:", count
            if (count > THRESH_PKT):
                #print prevtimestamp-basetimestamp, "count:", count
                s = np.append(s, prevtimestamp-basetimestamp)
                cs = np.append(cs, count)
            prevtimestamp = timestamp
            count = 1
    #print prevtimestamp, "count:", count
    if (count > THRESH_PKT):
        s = np.append(s, prevtimestamp-basetimestamp)
    lastfiletimestamp = timestamp
    lastcount = count
    s = np.unique(s)
    return s, cs

def parsepcap(filename):
    testcap = open(filename, 'rb')
    capfile = savefile.load_savefile(testcap, verbose=True)
    print "filename:", filename, "len(capfile.packets):", len(capfile.packets)
    ts = np.array([])
    for pkt in capfile.packets[1:]:
        timestamp = pkt.timestamp*1000 + pkt.timestamp_ms
        #print timestamp
        ts = np.append(ts, timestamp)
    return ts

def getax():
    fig = plt.figure()
    return fig.add_subplot(111)

def countFilesInDir(d):
    path, dirs, files = os.walk(d).next()
    return len(files)

def getCountStaistics(counts):
    print "mean count:", np.mean(counts), "std count:", np.std(counts)

def iterateDir(rootdir):
    filesubstr = "new_files"
    #"""
    events, counts = np.array([]), np.array([])
    d = countFilesInDir(rootdir)
    for i in range(0, d):
        filename = rootdir+filesubstr+str(i)
        print "analyzing file:", filename
        events = np.append(events, parsepcap(filename))
        #events, counts = parsepcapfilesingle(events, counts, filename)
        #events, counts = events.astype(int), counts.astype(int)
    #events = events - min(events)
    #np.save("/home/abhirup/Documents/Research/dataset/network_traffic/univ/univ2_trace/pt0/events", events)
    #plotsignal(convertToTimeSeq(events))
    return events

def getEvents():
    #for i in range(0,9):
    #    iterateDir("univ/univ2_trace/pt"+str(i)+"/")
    return iterateDir("/home/abhirup/Documents/Research/dataset/network_traffic/univ/univ2_trace/pt1/")
    #iterateDir("caida/dns20131002/")
    #plt.show()
    
if __name__ == "__main__":
    events = getEvents()
    np.save("network_data/events_pt1", events)