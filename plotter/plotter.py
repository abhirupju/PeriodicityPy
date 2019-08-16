import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math, os
from matplotlib.patches import Polygon
from scipy.ndimage.filters import gaussian_filter as gaussian_filter

matplotlib.rcParams['legend.numpoints'] = 1
matplotlib.rcParams['lines.linewidth'] = 3
matplotlib.rcParams['font.size'] = 20

def convertToTimeSeries(events, resolution=1):
    n = int(math.ceil(max(events)*resolution)+1)
    s = np.zeros(n)
    for i in range(0, len(events)):
        e = events[i]
        x = int(math.floor(e*resolution))
        s[x] = 1
    return s

def getax():
    fig = plt.figure()
    return fig.add_subplot(111)

def getNormalizedError(pred_periods, true_period):
    x = np.power(((pred_periods - true_period)/true_period),2)
    return x

def getSubError(pred_periods, true_period):
    xx = np.power(np.log((np.fabs(pred_periods - true_period)/true_period)),2)
    return xx

def getRatioError(pred_periods, true_period):
    xx = np.power(np.log((pred_periods/true_period)),2)
    print "pred_periods:", pred_periods
    print "errors:", xx
    return xx

def getRatioErrorBox(pred_periods, true_period):
    xx = np.fabs(np.log((pred_periods/true_period)))
    return xx

def getRMS(rs):
    rms = np.power(np.sum(rs, axis=0)/len(rs), 0.5)
    stderrs = np.std(rs, axis=0)/np.power(len(rs),0.5)
    return rms, stderrs

def convertToTimeSeq(events):
    n = int(max(events)+1)
    s = np.zeros(n)
    print "len(events):", len(events)
    for e in events:
        s[e] = 1
    return s

def plotData(data, errs, filename=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    width = max(data[:,0])/(len(data)*6)
    ax.bar(data[:,0],         data[:,1], width=width, yerr=errs[:,1], label="Particle Filter", color='b', edgecolor='black', ecolor='r')
    #TODO: Why this 1.5 times? :(
    xx = data[:,0]
    ax.bar(xx+1.5*width,   data[:,2], width=width,yerr=errs[:,2], label="Fourier Transform", color='g', edgecolor='black', ecolor='r', align='center')
    ax.bar(xx+2*width, data[:,3], width=width,yerr=errs[:,3], label="Histogram", color='m', edgecolor='black', ecolor='r')
    ax.bar(xx+3*width, data[:,4], width=width,yerr=errs[:,4], label="Autocorrelation", color='c', edgecolor='black', ecolor='r')
    #ax.bar(data[:,0]+4*width, data[:,5], width=width,yerr=errs[:,5], label="FFT AutoCorr combined", color='y', edgecolor='black', ecolor='r')
    ax.set_xticks(xx+1.5*width)
    #ax.set_ylim([0,30])
    ymin, ymax = ax.get_ylim()
    ax.set_ylim([0,ymax])
    #ax.set_yscale("log")
    print "filename:", filename
    if (filename.split("_")[0] == 'noise'):
        ax.set_xticklabels(data[:,0]*10)
        ax.set_xlabel("#Noise events : #Periodic events")
    elif (filename.split("_")[0] == 'std'):
        ax.set_xlabel("Standard deviation")
        ax.set_xticklabels(data[:,0])
    ax.set_ylabel("Error")
    #ax.set_title(filename[:filename.rfind('.')])
    xmin, xmax = ax.get_xlim()
    #ax.plot([xmin,xmax],[10,10])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    # Put a legend below current axis
    lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2)
    #if (filename == None):
    savefile(filename[:filename.rfind('.')], lgd)
    #plt.show()

def getFormattedData(d):
    spread = d
    center = np.ones(len(d)) * np.median(d)
    high = np.ones(len(d)) * np.percentile(d, 75)
    low = np.ones(len(d)) * np.percentile(d, 25)
    return np.concatenate((spread, center, high, low), 0)

def addArrays(x1,x2):
    if (len(x1) < len(x2)):
        return x1 + x2[:len(x1)]
    elif (len(x2) < len(x1)):
        return x2 + x1[:len(x2)]
    else:
        return x1 + x2

def cpuTimePlot(rootdir):
    pfdata, fftdata, segdata, autocorrdata = np.array([]), np.array([]), np.array([]), np.array([])
    iteration = 0
    for subdir, dirs, files in os.walk(rootdir):
        if (subdir==rootdir):
            for filename in files:
                xdata  = np.load(subdir+filename)
                xx = xdata
                if ("pf" in filename):
                    iteration+=1
                    if (len(pfdata) == 0):
                        pfdata = xx
                    else:
                        pfdata = addArrays(pfdata, xx)
                if ("fft" in filename):
                    if (len(fftdata) == 0):
                        fftdata = xx
                    else:
                        fftdata = addArrays(fftdata, xx)
                if ("seq" in filename):
                    if (len(segdata) == 0):
                        segdata = xx
                    else:
                        segdata = addArrays(segdata, xx)
                if ("auto" in filename):
                    if (len(autocorrdata) == 0):
                        autocorrdata = xx
                    else:
                        autocorrdata = addArrays(autocorrdata, xx)

    pfdata = pfdata/iteration
    fftdata = fftdata/iteration
    segdata = segdata/iteration
    autocorrdata = autocorrdata/iteration
    ax=getax()
    ax.plot(pfdata, label="Particle Filter")
    ax.plot(fftdata, label="Fourier Transform")
    ax.plot(segdata, label="Segmentation Algorithm")
    ax.plot(autocorrdata, label="Autocorrelation")
    ax.set_xlabel("Number of observed events")
    ax.set_ylabel("RMS of Relative Period \n Prediction Error($e$)", multialignment='center')
    #ax.set_xlim([0,150])
    #ax.set_xticks(range(0,151, 50))
    #ax.set_yticks([0,1,2,3])
    ax.set_yscale("log")
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    # Put a legend below current axis
    lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2)
    savefile('cputime', lgd)
    plt.show()

def boxPlotData4(data, filename=None):
    periods = np.unique(data[:,0])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hold(True)
    i = 0
    pdata = []
    for i in range(len(periods)):
        p = periods[i]
        d = np.array([])
        print "period:", p
        for row in range(len(data)):
            if (data[row][0] == p):
                if (len(d) == 0):
                    d = getRatioErrorBox(data[row], p)
                else:
                    d = np.vstack((d, getRatioErrorBox(data[row], p)))
        pf_data = getFormattedData(d[:,1])
        fft_data = getFormattedData(d[:,2])
        hist_data = getFormattedData(d[:,3])
        autocorr_data = getFormattedData(d[:,4])
        pf_data.shape = (-1, 1)
        fft_data.shape = (-1, 1)
        hist_data.shape = (-1, 1)
        autocorr_data.shape = (-1, 1)
        pdata.append(pf_data)
        pdata.append(fft_data)
        pdata.append(hist_data)
        pdata.append(autocorr_data)
    bp = ax.boxplot(pdata)
    #ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    mmxx = np.arange(len(periods))*4+4.5
    bottom, top = ax.get_ylim()
    bottom=-1
    for mmmxx in mmxx:
        ax.plot([mmmxx,mmmxx],[bottom,top], color='lightgray', alpha=0.5)
    #plt.show()
    #return
    boxColors = ['blue', 'green', 'c', 'm', 'y']
    numBoxes = 4*len(periods)
    medians = list(range(numBoxes))
    for i in range(numBoxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = list(zip(boxX, boxY))
        # Alternate between Dark Khaki and Royal Blue
        k = i % 4
        boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
        ax.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            plt.plot(medianX, medianY, 'k')
            medians[i] = medianY[0]
    
    # Set the axes ranges and axes labels
    ax.set_xlim(0.5, numBoxes + 0.5)
    ax.set_ylim(bottom, top)
    ax.set_xticks(np.arange(len(periods))*4 + 2.5)
    #ax.set_xticklabels(periods)
    xtickNames = plt.setp(ax, xticklabels=periods.astype(int))
    #ax.set_xlabel("$\\frac{Number\;of\;background\;events}{Number\;of\;periodic\;events}$")
    #ax.set_xlabel("Relative deviation of Periodic events $(\\frac{\sigma}{T})$", multialignment='center')
    ax.set_xlabel("True Period $(T)$")
    ax.set_ylabel("Relative Period \n Prediction Error($e$)", multialignment='center')
    #plt.setp(xtickNames, fontsize=40)
    import matplotlib.patches as mpatches
    p1 = mpatches.Rectangle((0, 0), 1, 1, color='blue')
    p2 = mpatches.Rectangle((0, 0), 1, 1, color='g')
    p3 = mpatches.Rectangle((0, 0), 1, 1, color='c')
    p4 = mpatches.Rectangle((0, 0), 1, 1, color='m')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    # Put a legend below current axis
    #lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2)

    lgd = ax.legend([p1, p2, p3, p4], ['Particle Filter', 'Fourier Transform', 'Histogram', 'Autocorrelation'],
               loc='upper center', bbox_to_anchor=(0.5, -0.3), fancybox=True, ncol=2, fontsize=20)
    savefile('period_boxplot', lgd)
    plt.show()

def boxPlotData5(data, filename=None):
    periods = np.unique(data[:,0])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hold(True)
    i = 0
    pdata = []
    for i in range(len(periods)):
        p = periods[i]
        d = np.array([])
        print "period:", p
        for row in range(len(data)):
            if (data[row][0] == p):
                if (len(d) == 0):
                    d = getRatioErrorBox(data[row], 10)
                else:
                    d = np.vstack((d, getRatioErrorBox(data[row],10)))
        pf_data = getFormattedData(d[:,1])
        fft_data = getFormattedData(d[:,2])
        hist_data = getFormattedData(d[:,3])
        autocorr_data = getFormattedData(d[:,4])
        fftautocorr_data = getFormattedData(d[:,5])
        pf_data.shape = (-1, 1)
        fft_data.shape = (-1, 1)
        hist_data.shape = (-1, 1)
        autocorr_data.shape = (-1, 1)
        fftautocorr_data.shape = (-1, 1)
        pdata.append(pf_data)
        pdata.append(fft_data)
        pdata.append(hist_data)
        pdata.append(autocorr_data)
        pdata.append(fftautocorr_data)
    bp = ax.boxplot(pdata)
    #ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    mmxx = np.arange(len(periods))*5+5.5
    bottom, top = ax.get_ylim()
    bottom=-1
    for mmmxx in mmxx:
        ax.plot([mmmxx,mmmxx],[bottom,top], color='lightgray', alpha=0.5)
    #plt.show()
    #return
    boxColors = ['blue', 'green', 'c', 'm', 'y']
    numBoxes = 5*len(periods)
    medians = list(range(numBoxes))
    for i in range(numBoxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = list(zip(boxX, boxY))
        # Alternate between Dark Khaki and Royal Blue
        k = i % 5
        boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
        ax.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            plt.plot(medianX, medianY, 'k')
            medians[i] = medianY[0]
    
    # Set the axes ranges and axes labels
    ax.set_xlim(0.5, numBoxes + 0.5)
    ax.set_ylim(bottom, top)
    ax.set_xticks(np.arange(len(periods))*5 + 2.5)
    #ax.set_xticklabels(periods)
    xtickNames = plt.setp(ax, xticklabels=periods)
    ax.set_xlabel("$\\frac{Number\;of\;background\;events}{Number\;of\;periodic\;events}$")
    #ax.set_xlabel("$\\frac{\sigma}{T}$", fontsize=40)
    #ax.set_xlabel("$T$")
    ax.set_ylabel("Relative Period \n Prediction Error($e$)", multialignment='center')
    #plt.setp(xtickNames, fontsize=40)
    import matplotlib.patches as mpatches
    p1 = mpatches.Rectangle((0, 0), 1, 1, color='blue')
    p2 = mpatches.Rectangle((0, 0), 1, 1, color='g')
    p3 = mpatches.Rectangle((0, 0), 1, 1, color='c')
    p4 = mpatches.Rectangle((0, 0), 1, 1, color='m')
    p5 = mpatches.Rectangle((0, 0), 1, 1, color='y')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    # Put a legend below current axis
    #lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2)

    lgd = ax.legend([p1, p2, p3, p4, p5], ['Particle Filter', 'Fourier Transform', 'Histogram', 'Autocorrelation', 'Fourier Transform \n Autocorrelation combination'],
               loc='upper center', bbox_to_anchor=(0.5, -0.3), fancybox=True, ncol=2, fontsize=20)
    savefile('noise_boxplot', lgd)
    plt.show()

def boxplotlegend():
    import matplotlib.patches as mpatches
    fig=plt.figure()
    ax = fig.add_subplot(111)
    p1 = mpatches.Rectangle((0, 0), 1, 1, color='blue')
    p2 = mpatches.Rectangle((0, 0), 1, 1, color='g')
    p3 = mpatches.Rectangle((0, 0), 1, 1, color='c')
    p4 = mpatches.Rectangle((0, 0), 1, 1, color='m')
    p5 = mpatches.Rectangle((0, 0), 1, 1, color='y')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    # Put a legend below current axis
    #lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2)

    lgd = ax.legend([p1, p2, p3, p4, p5], ['Particle Filter', 'Fourier Transform', 'Histogram', 'Autocorrelation', 'Fourier Transform \n Autocorrelation combination'],
               loc='upper center', bbox_to_anchor=(0.5, -0.3), fancybox=True, ncol=1, fontsize=20)
    savefile('legend_boxplot', lgd)
    plt.show()
    
def boxPlotnumParticles(data):
    numParticles = np.unique(data[:,0])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hold(True)
    i = 0
    pdata = []
    for i in range(len(numParticles)):
        p = numParticles[i]
        d = np.array([])
        print "numParticles:", p
        for row in range(len(data)):
            if (data[row][0] == p):
                if (len(d) == 0):
                    d = getRatioErrorBox(data[row], 10)
                else:
                    d = np.vstack((d, getRatioErrorBox(data[row],10)))
        pf_data = getFormattedData(d[:,1])
        pf_data.shape = (-1, 1)
        pdata.append(pf_data)
    bp = ax.boxplot(pdata)
    for box in bp['boxes']:
        box.set(linewidth=3)

    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(linewidth=3)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(linewidth=3)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(linewidth=3)
    ax.set_yticks([0.0,0.25,0.5,1.0,1.5])
    tlabels = []
    for i in numParticles:
        tlabels.append("$2^{"+str(int(np.log2(i)))+"}$")
    #print "tlabels:", tlabels
    ax.set_xticks(np.arange(1, len(numParticles)+1))
    xtickNames = plt.setp(ax, xticklabels=tlabels)
    ax.tick_params(axis='both', which='major', pad=10)
    #plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')
    ax.set_xlabel("Number of Particles")
    ax.set_ylabel("Relative error($\\xi$)", multialignment='center')
    savefile('numParticles_boxplot', None)
    #plt.show()

def periodChangePlot(legend_font_size):
    pfdata = np.load("../results/periodVary_pfperiods28.npy")
    fftdata = np.load("../results/periodVary_fft.npy")
    segdata = np.load("../results/periodVary_segment.npy")
    autocorrdata = np.load("../results/periodVary_autocorr.npy")

    info = np.load("../results/periodVary_info.npy")
    ax = getax()
    x = np.arange(0, len(pfdata), 5).astype(int)

    print "lens:", len(pfdata), len(fftdata), len(autocorrdata)

    #ax.plot(x, pf_data[x], 'o-', color='g', label='Estimated period')

    ax.set_xlim([0, 260])
    tr = np.array([])
    for r in info:
        print r[0], r[1]
        tr = np.append(tr, np.ones(int(r[1]))*r[0])
    ax.plot(tr,linestyle='--', color='black', label='True period')
    ax.plot(x, pfdata[x], '-', color='blue', label="Particle Filter")
    ax.plot(x, fftdata[x], '--',color='green', label="Fourier Transform")
    ax.plot(x, segdata[x], ':',color='red', label="Segmentation algorithm")
    ax.plot(x, autocorrdata[x], '-.',color='gray', label="Autocorrelation")

    ax.set_xlabel("Number of observations")
    ax.set_ylabel("Period")

    ax.yaxis.set_ticks(np.array([5, 10, 20, 30, 40]))

    ax.set_ylim([0, 45])

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    # Put a legend below current axis
    lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2, fontsize=legend_font_size)
    #if (filename == None):
    savefile("changePeriod")
    #plt.show()

def saveplotfile(figure, filename, transparent=False, lgd=None, extension="pdf"):
    print ("saving:", filename)
    if (extension == "pdf"):
        if (lgd!=None):
            figure.savefig(filename+'.pdf', format='pdf', dpi=1000, transparent=True, bbox_extra_artists=(lgd,), bbox_inches='tight')
        else:
            figure.savefig(filename+'.pdf', format='pdf', dpi=1000, transparent=True, bbox_inches='tight')
    elif (extension == "png"):
        if (lgd!=None):
            figure.savefig(filename+'.png', format='png', dpi=600, transparent=transparent, bbox_extra_artists=(lgd,), bbox_inches='tight')
        else:
            figure.savefig(filename+'.png', format='png', dpi=600, transparent=transparent, bbox_inches='tight')

def savefile(filename, lgd=None):
    if (lgd!=None):
        plt.savefig(filename+'.eps', format='eps', dpi=1600, bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.savefig(filename+'.png', format='png', dpi=1000, bbox_inches='tight')

def getErrors(T, data, errformula="ratio"):
    rates = np.unique(data[:,0])
    errors, stderrs = np.array([]), np.array([]) 
    for r in rates:
        print "rate:", r
        x = data[np.where(data[:,0]==r)]
        slicedx = x[0:len(x),1:len(x)]
        #print "rate:", r
        if (errformula == "ratio"):
            ers = getRatioError(slicedx, T)
        elif(errformula == "normalized"):
            ers = getNormalizedError(slicedx, T)
        else:
            ers = slicedx
        err, errbar = getRMS(ers)
        err = ers
        if (len(errors) == 0):
            errors = np.insert(err,0,r)
            stderrs = np.insert(errbar,0,r)
        else:
            errors = np.vstack((errors, np.insert(err,0,r)))
            stderrs = np.vstack((stderrs, np.insert(errbar,0,r)))
    return errors, stderrs

def onlinePlot(rootdir):
    pfdata, fftdata, segdata, autocorrdata = np.array([]), np.array([]), np.array([]), np.array([])
    iteration = 0
    for subdir, dirs, files in os.walk(rootdir):
        if (subdir==rootdir):
            for filename in files:
                xdata  = np.load(subdir+filename)
                #print "filename:", filename, xdata
                xx = getRatioError(xdata,10)
                #xx = xdata
                if ("pf" in filename):
                    iteration+=1
                    if (len(pfdata) == 0):
                        pfdata = xx
                    else:
                        pfdata = addArrays(pfdata, xx)
                if ("fft" in filename):
                    if (len(fftdata) == 0):
                        fftdata = xx
                    else:
                        fftdata = addArrays(fftdata, xx)
                if ("seg" in filename):
                    if (len(segdata) == 0):
                        segdata = xx
                    else:
                        segdata = addArrays(segdata, xx)
                if ("auto" in filename):
                    if (len(autocorrdata) == 0):
                        autocorrdata = xx
                    else:
                        autocorrdata = addArrays(autocorrdata, xx)
    #pfdata = pfdata/iteration
    #fftdata = fftdata/iteration
    #histdata = histdata/iteration
    #autocorrdata = autocorrdata/iteration
    x = np.arange(0, 162, 6)
    pfdata = np.power(pfdata/iteration,0.5)[x]
    fftdata = np.power(fftdata/iteration,0.5)[x]
    segdata = np.power(segdata/iteration,0.5)[x]
    autocorrdata = np.power(autocorrdata/iteration,0.5)[x]
    ax=getax()
    ax.plot(x, pfdata, 'o-', color='blue', label="Particle Filter")
    ax.plot(x, fftdata, '*--',color='green', label="Fourier Transform")
    ax.plot(x, segdata, 's:',color='red', label="Segmentation algorithm")
    ax.plot(x, autocorrdata, 'x-.',color='gray', label="Autocorrelation")
    ax.set_xlabel("Number of observations")
    ax.set_ylabel("RMS of Relative error($\\xi$)")
    ax.set_xlim([0,160])
    ax.set_xticks(range(0,161, 40))
    #ax.set_yticks([0,1,2,3])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    # Put a legend below current axis
    lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2)
    savefile('online', lgd)
    plt.show()


def stepratePFCDFPlot():
    baseFolder = "/home/abhirup/Documents/Research/dataset/steprate/population/signals_PFCDF/"
    baseFolder2 = "/home/abhirup/Documents/Research/dataset/steprate/population/signals/"
    
    pf15Error, pfError, pf2Error, pf4Error = np.array([]), np.array([]), np.array([]), np.array([])
    for folder in os.listdir(baseFolder):
        for fname in os.listdir(baseFolder+folder):
            if (fname.endswith("trimmed")==True):
                fname = fname+"_bins1_std"
                basePath = baseFolder+folder+"/"
                pf15periods = np.load(basePath+fname+"15.npy")
                pfperiods = np.load(basePath+fname+".npy")
                pf2periods = np.load(basePath+fname+"2.npy")
                pf4periods = np.load(basePath+fname+"4.npy")
                
                
                events, gtruthX, = np.array([]), np.array([])
                with open(baseFolder2+folder+"/"+fname) as infile:
                    for line in infile:
                            events = np.append(events, float(line.replace('\n','').replace('\r','')))
                xt = np.loadtxt('/home/abhirup/Documents/Research/dataset/steprate/population/'+fname[:fname.find("trimmed")]+"ground_truth", delimiter='\t')
                prevlx = 0
                for x in xt:
                    lx = np.where(events > x[0]/1000.0)[0]
                    if (len(lx) > 0):
                        xx = min(lx)
                    else:
                        xx = len(events)
                    gtruthX = np.append(gtruthX, np.ones(xx-prevlx)*(x[1]/1000.0))
                    prevlx = xx
                print len(pf15periods), len(pfperiods), len(gtruthX)
                """
                gtruthX = gtruthX[2:]
                pf15periods = pf15periods[2:]
                pfperiods = pfperiods[2:]
                pf2periods = pf2periods[2:]
                pf4periods = pf4periods[2:]
                """
                for i in range(len(pf15periods)):
                    pf15Error = np.append(pf15Error, getRatioErrorBox(pf15periods[i], gtruthX[i]))
                    pfError = np.append(pfError, getRatioErrorBox(pfperiods[i], gtruthX[i]))
                    pf2Error = np.append(pf2Error, getRatioErrorBox(pf2periods[i], gtruthX[i]))
                    pf4Error = np.append(pf4Error, getRatioErrorBox(pf4periods[i], gtruthX[i]))
                print len(pf15Error), len(pfError), len(pf2Error)
    print "--------------------------"
    print len(pf15Error), len(pfError), len(pf2Error)
    pf15Error_sorted = np.sort(pf15Error)
    p_pf15Error = 1. * np.arange(len(pf15Error)) / (len(pf15Error) - 1)
                        
    pfError_sorted = np.sort(pfError)
    p_pfError = 1. * np.arange(len(pfError)) / (len(pfError) - 1)
    
    pf2Error_sorted = np.sort(pf2Error)
    p_pf2Error = 1. * np.arange(len(pf2Error)) / (len(pf2Error) - 1)
    
    pf4Error_sorted = np.sort(pf4Error)
    p_pf4Error = 1. * np.arange(len(pf4Error)) / (len(pf4Error) - 1)
    
    ax = getax()
    ax.plot(pf15Error_sorted,p_pf15Error*100.0, color='black', label='$1.5 \\times \sigma(signal)$')
    ax.plot(pfError_sorted,p_pfError*100.0, color='blue', label='$\sigma(signal)$')
    ax.plot(pf2Error_sorted,p_pf2Error*100.0, color='g', label='$\sigma(signal)/2$')
    ax.plot(pf4Error_sorted,p_pf4Error*100.0, color='r', label='$\sigma(signal)/4$')
    
    ax.set_xlabel('Relative Error($\\xi$)')
    ax.set_ylabel('% of observations')
    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    # Put a legend below current axis
    lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2, fontsize=20)
    #if (filename == None):
    savefile("stepDataPFCDF")

    #plt.show()
    
def steprateCDFPlot(legend_font_size):
    baseFolder = "/home/abhirup/Documents/Research/dataset/steprate/population/signals/"
    pfError, fftError, autoCorrError, segError = np.array([]), np.array([]), np.array([]), np.array([])
    for folder in os.listdir(baseFolder):
        for fname in os.listdir(baseFolder+folder):
            if ("std2" in fname and "npy" not in fname and "png" not in fname and "_period" not in fname):

                basePath = baseFolder+folder+"/"
                estimatedperiods = np.load(basePath+fname+".npy")
                fftperiods = np.load(basePath+fname+"_FFT.npy")
                autocorrperiods = np.load(basePath+fname+"_AutoCorr.npy")
                segperiods = np.load(basePath+fname+"_SEG.npy")
                
                events, gtruthX, = np.array([]), np.array([])
                with open(basePath+fname) as infile:
                    for line in infile:
                            events = np.append(events, float(line.replace('\n','').replace('\r','')))
                xt = np.loadtxt('/home/abhirup/Documents/Research/dataset/steprate/population/'+fname[:fname.find("trimmed")]+"ground_truth", delimiter='\t')
                prevlx = 0
                for x in xt:
                    lx = np.where(events > x[0]/1000.0)[0]
                    if (len(lx) > 0):
                        xx = min(lx)
                    else:
                        xx = len(events)
                    gtruthX = np.append(gtruthX, np.ones(xx-prevlx)*(x[1]/1000.0))
                    prevlx = xx
                
                gtruthX = gtruthX[2:]
                estimatedperiods = estimatedperiods[2:]
                for i in range(len(gtruthX)):
                    pfError = np.append(pfError, getRatioErrorBox(estimatedperiods[i], gtruthX[i]))
                    fftError = np.append(fftError, getRatioErrorBox(fftperiods[i], gtruthX[i]))
                    autoCorrError = np.append(autoCorrError, getRatioErrorBox(autocorrperiods[i], gtruthX[i]))
                    segError = np.append(segError, getRatioErrorBox(segperiods[i], gtruthX[i]))
    pfError_sorted = np.sort(pfError)
    p_pfError = 1. * np.arange(len(pfError)) / (len(pfError) - 1)
    
    fftError_sorted = np.sort(fftError)
    p_fftError = 1. * np.arange(len(fftError)) / (len(fftError) - 1)
    
    segError_sorted = np.sort(segError)
    p_segError = 1. * np.arange(len(segError)) / (len(segError) - 1)
    
    autoCorrError_sorted = np.sort(autoCorrError)
    p_autoCorrError = 1. * np.arange(len(autoCorrError)) / (len(autoCorrError) - 1)
    
    ax = getax()
    
    ax.plot(pfError_sorted, (p_pfError)*100.0, color='blue', label='Particle Filter')
    ax.plot(fftError_sorted, p_fftError*100.0, '--', color='g', label='Fourier Transform')
    ax.plot(autoCorrError_sorted, p_autoCorrError*100.0, '-.', color='grey', label='Auto-correlation')
    ax.plot(segError_sorted, p_segError*100.0, ':', color='r', label='Segmentation Algorithm')
    
    ax.set_xlabel('Relative Error($\\xi$)')
    ax.set_ylabel('% of observations')
    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height])
    # Put a legend below current axis
    lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2, fontsize=legend_font_size)
    #if (filename == None):
    savefile("stepDataCDF")

    #plt.show()

def steprateTrace11Plot(legend_font_size):
    baseFolder = "/home/abhirup/Documents/Research/dataset/steprate/population/signals/"
    for folder in os.listdir(baseFolder):
        #print "folder:"+folder
        for fname in os.listdir(baseFolder+folder):
            if ("std" in fname and "npy" not in fname and "png" not in fname and "_period" not in fname):
                
                fname = "sensor_event_file_1479914987957_trimmed_bins1_std"
                folder="11"
                
                #Plot Signal
                ax = getax()
                rawsignal = np.loadtxt(baseFolder+folder+"/sensor_event_file_1479914987957_trimmed")
                magsignal = []
                for s in rawsignal:
                    magsignal.append(math.sqrt(s[1]**2+s[2]**2+s[3]**2))
                normalizedsignal = magsignal-np.mean(magsignal)
                sigma=4
                filteredsignal = gaussian_filter(normalizedsignal, sigma)
                std = np.std(filteredsignal)
                ts = (rawsignal[:,0]-rawsignal[0,0])/1000.0
                ax.plot(ts, filteredsignal, linewidth=1, color='black')
                ax.plot(ts, np.ones(len(ts))*std*1.5, '--', color='black', label='$1.5 \\times \sigma(signal)$')
                ax.plot(ts, np.ones(len(ts))*std/1.0, '--', color='blue', label='$\sigma(signal)$')
                ax.plot(ts, np.ones(len(ts))*std/2.0, '--', color='g', label='$\sigma(signal)/2$')
                ax.plot(ts, np.ones(len(ts))*std/4.0, '--', color='red', label='$\sigma(signal)/4$')
                ax.set_yticklabels([])
                ax.set_xlabel('Time (Seconds)')
                ax.set_xlim([10, 35])
                ax.set_ylim([min(filteredsignal), max(filteredsignal)])
                box = ax.get_position()
                ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
                lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2)
                #savefile("stepDataTrace11")
                #return

                basePath = baseFolder+folder+"/"
                print "Processing:"+basePath+fname
                estimatedperiods = np.load(basePath+fname+".npy")
                fftperiods = np.load(basePath+fname+"_FFT.npy")
                autocorrperiods = np.load(basePath+fname+"_AutoCorr.npy")
                segperiods = np.load(basePath+fname+"_SEG.npy")
                
                events, gtruthX, = np.array([]), np.array([])
                print "eventFile:", (basePath+fname) 
                with open(basePath+fname) as infile:
                    for line in infile:
                            events = np.append(events, float(line.replace('\n','').replace('\r','')))
                xt = np.loadtxt('/home/abhirup/Documents/Research/dataset/steprate/population/'+fname[:fname.find("trimmed")]+"ground_truth", delimiter='\t')
                prevlx = 0
                for x in xt:
                    lx = np.where(events > x[0]/1000.0)[0]
                    if (len(lx) > 0):
                        xx = min(lx)
                    else:
                        xx = len(events)
                    print "Add ", x[1], "in", prevlx, xx
                    gtruthX = np.append(gtruthX, np.ones(xx-prevlx)*(x[1]/1000.0))
                    prevlx = xx
                pfError, fftError, autoCorrError, segError = np.array([]), np.array([]), np.array([]), np.array([])
                gtruthX = gtruthX[2:]
                estimatedperiods = estimatedperiods[2:]
                print "Lengths:", len(events), len(gtruthX), len(estimatedperiods), len(fftperiods), len(autocorrperiods), len(segperiods)
                for i in range(len(gtruthX)):
                    pfError = np.append(pfError, getRatioErrorBox(estimatedperiods[i], gtruthX[i]))
                    fftError = np.append(fftError, getRatioErrorBox(fftperiods[i], gtruthX[i]))
                    autoCorrError = np.append(autoCorrError, getRatioErrorBox(autocorrperiods[i], gtruthX[i]))
                    segError = np.append(segError, getRatioErrorBox(segperiods[i], gtruthX[i]))
            
                #ax.plot(gtruthX, '--', color='g', label='Ground Truth')
                
                fig = plt.figure(figsize=(6,5))
                ax = fig.add_subplot(111)
                ax.plot(pfError, '-*', color='blue', label='Particle Filter')
                ax.plot(fftError, '-', marker='*', color='g', label='Fourier Transform')
                ax.plot(autoCorrError, '-.', marker='x', color='grey', label='Autocorrelation')
                ax.plot(segError, '--',  marker='s', color='r', label='Segmentation algorithm')
                ax.set_ylabel("Relative Error($\\xi$)")
                ax.set_xlabel("Number of observations")
                ax.set_xlim([0, len (gtruthX)])
                ax.xaxis.set_ticks([0, 20, 40, 60, 80])
                box = ax.get_position()
                ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width*1.2, box.height * 1.1])
                lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2, fontsize=legend_font_size)
                #savefile("stepError")

                #fig.savefig("Trace11_error.png")   # save the figure to file
                #plt.close(fig)
                
                matplotlib.rcParams['font.size'] = 19
                
                fig = plt.figure(figsize=(6,5))
                ax = fig.add_subplot(111)
                stride = max( int(len(gtruthX) / 10), 1)
                ax.plot(gtruthX, '--', marker='^', markersize=15, markevery=stride, color='black', label='Ground Truth')
                #ax.plot(estimatedperiods, color='blue', label='Particle Filter')
                #ax.plot(fftperiods, '--', marker='*', color='g', label='Fourier Transform')
                #ax.plot(autocorrperiods, '-.', marker='x', color='grey', label='Autocorrelation')
                #ax.plot(segperiods, ':', marker='s', color='r', label='Segmentation algorithm')
                
                ax.plot(estimatedperiods, color='blue', label='Particle Filter')
                ax.plot(fftperiods, color='g', label='Fourier Transform')
                ax.plot(autocorrperiods, color='brown', label='Autocorrelation')
                ax.plot(segperiods, color='r', label='Segmentation algorithm')
                
                ax.set_ylabel("Estimated Period (seconds)")
                ax.set_xlabel("Number of observations")
                ax.set_xlim([0, len (gtruthX)])
                ax.yaxis.set_ticks([0.5, 1.0, 2.0, 3.0, 4.0])
                ax.xaxis.set_ticks([0, 20, 40, 60, 80])
                ax.set_ylim([0,4.5])
                #box = ax.get_position()
                #ax.set_position([box.x0-0.5, box.y0 + box.height * 0.1, box.width, box.height*0.9])
                #lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2, fontsize=15)
                savefile("stepEstimatedPeriod-full")
                #plt.close(fig)
                #plt.show()
                return

def stepratePlot():
    baseFolder = "/home/abhirup/Documents/Research/dataset/steprate/population/signals/"
    for folder in os.listdir(baseFolder):
        #print "folder:"+folder
        for fname in os.listdir(baseFolder+folder):
            if ("std" in fname and "npy" not in fname and "png" not in fname and "_period" not in fname):
                
                #fname = "sensor_event_file_1479914926610_trimmed_bins1_std"
                #folder="10"

                basePath = baseFolder+folder+"/"
                print "Processing:"+basePath+fname
                estimatedperiods = np.load(basePath+fname+".npy")
                fftperiods = np.load(basePath+fname+"_FFT.npy")
                autocorrperiods = np.load(basePath+fname+"_AutoCorr.npy")
                segperiods = np.load(basePath+fname+"_SEG.npy")
                
                events, gtruthX, = np.array([]), np.array([])
                print "eventFile:", (basePath+fname) 
                with open(basePath+fname) as infile:
                    for line in infile:
                            events = np.append(events, float(line.replace('\n','').replace('\r','')))
                xt = np.loadtxt('/home/abhirup/Documents/Research/dataset/steprate/population/'+fname[:fname.find("trimmed")]+"ground_truth", delimiter='\t')
                prevlx = 0
                for x in xt:
                    lx = np.where(events > x[0]/1000.0)[0]
                    if (len(lx) > 0):
                        xx = min(lx)
                    else:
                        xx = len(events)
                    print "Add ", x[1], "in", prevlx, xx
                    gtruthX = np.append(gtruthX, np.ones(xx-prevlx)*(x[1]/1000.0))
                    prevlx = xx
                pfError, fftError, autoCorrError, segError = np.array([]), np.array([]), np.array([]), np.array([])
                gtruthX = gtruthX[2:]
                estimatedperiods = estimatedperiods[2:]
                print "Lengths:", len(events), len(gtruthX), len(estimatedperiods), len(fftperiods), len(autocorrperiods), len(segperiods)
                for i in range(len(gtruthX)):
                    pfError = np.append(pfError, getRatioErrorBox(estimatedperiods[i], gtruthX[i]))
                    fftError = np.append(fftError, getRatioErrorBox(fftperiods[i], gtruthX[i]))
                    autoCorrError = np.append(autoCorrError, getRatioErrorBox(autocorrperiods[i], gtruthX[i]))
                    segError = np.append(segError, getRatioErrorBox(segperiods[i], gtruthX[i]))
            
                #ax.plot(gtruthX, '--', color='g', label='Ground Truth')
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(pfError, '-*', color='blue')
                ax.plot(fftError, '-', color='green')
                ax.plot(autoCorrError, '-^', color='c')
                ax.plot(segError, '--', color='m')
                ax.set_ylabel("Relative Error($\\xi$)")
                fig.savefig(basePath+fname.split('_')[-1]+"_error.png")   # save the figure to file
                plt.close(fig)
                
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(gtruthX, '--', color='black')
                ax.plot(estimatedperiods, '-*', color='blue')
                ax.plot(fftperiods, '-', color='green')
                ax.plot(autocorrperiods, '-^', color='c')
                ax.plot(segperiods, '--', color='m')
                ax.set_ylabel("Estimated period")
                fig.savefig(basePath+fname.split('_')[-1]+"_period.png")   # save the figure to file
                plt.close(fig)
                
                
                #savefile(basePath+fname+"_npy")
                #plt.legend(loc='best')
                #plt.show()
                #return

def pulsePredictErrorPlot():
    baseFolder="/home/abhirup/Documents/Research/dataset/pulserate/signal/"
    for folder in os.listdir(baseFolder):
        folder="4"
        events = np.loadtxt(baseFolder+folder+"/s"+folder+".csv")[2:]
        pf = np.load(baseFolder+folder+"/s"+folder+".csv.npy")[2:]
        autocorr = np.load(baseFolder+folder+"/s"+folder+".csv_AutoCorr.npy")
        fft = np.load(baseFolder+folder+"/s"+folder+".csv_FFT.npy")
        seg = np.load(baseFolder+folder+"/s"+folder+".csv_SEG.npy")
        
        print len(events), len(pf), len(autocorr), len(fft), len(seg)
        
        f, ax = plt.subplots()
        predErrors, autocorrErrors, fftErrors, segErrors = [],[],[],[]
        for i in range(len(events)-1):
            #ax.axvline(x=events[i], ymin=0, ymax=0.8, color='black')
            #ax.axvline(x=events[i] + pf[i], ymin=0, ymax=0.8, color='red')
            predErrors.append(abs(events[i] + pf[i] - events[i+1]))
            
            autocorrErrors.append(abs(events[i] + autocorr[i-2] - events[i+1]))
            fftErrors.append(abs(events[i] + fft[i-2] - events[i+1]))
            segErrors.append(abs(events[i] + seg[i-2] - events[i+1]))
        print len(predErrors), len(autocorrErrors), len(fftErrors), len(segErrors)
        #ax.set_yticklabels([])
        ax.plot(events[1:], predErrors, color='blue', label='Particle Filter')
        ax.plot(events[1:], fftErrors,color='green', label='Fourier Transform')
        ax.plot(events[1:], autocorrErrors,color='gray', label='Autocorrelation')
        ax.plot(events[1:], segErrors, color='red', label='Segmentation Algorithm')
        plt.xlabel('Time (seconds)')
        plt.ylabel('|Predicted-actual| sec')
        plt.xlim([min(events)+1, max(events)])
        plt.yticks([1, 2, 4, 6])
        #plt.ylabel('|predicted - actual| (secs)')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        # Put a legend below current axis
        lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2, fontsize=18)

        savefile('pulse-prediction-error', lgd)
        
        #imgName=baseFolder+"s"+folder+"_predictions.png"
        #plt.tight_layout()
        #plt.savefig(imgName)
        break
    #plt.show()

def plot_different_distributions():
    normal_period = np.load('../results/diffdistributions/normal.npy')
    laplace_period = np.load('../results/diffdistributions/laplace.npy')
    gamma_period = np.load('../results/diffdistributions/gamma.npy')
    
    normal_errors = np.fabs(np.log((normal_period/10)))
    laplace_errors = np.fabs(np.log((laplace_period/10)))
    gamma_errors = np.fabs(np.log((gamma_period/10)))
    f, ax = plt.subplots()
    ax.plot(range(len(normal_errors)), normal_errors, color='blue', label='Normal distribution')
    ax.plot(range(len(laplace_errors)), laplace_errors, '--', color='green', label='Laplace distribution')
    ax.plot(range(len(gamma_errors)), gamma_errors, '-.', color='red', label='Gamma distribution')
    
    ax.set_ylabel("Relative error($\\xi$)")
    ax.set_xlabel("Number of observations")
    
    ax.set_xlim([0,45])
    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, ncol=2, fontsize=18)
    saveplotfile(f, 'different-distributions', lgd)
    plt.show()
    
if __name__ == "__main__":
    matplotlib.rcParams['font.size'] = 22
    #stepratePlot()
    #steprateCDFPlot(legend_font_size=17)
    #stepratePFCDFPlot()
    
    #pulsePredictErrorPlot()
    #steprateTrace11Plot(legend_font_size=17)
    
    plot_different_distributions()
    #periodChangePlot(legend_font_size=17)
    #onlinePlot("../results/online/")
    #plotNetworkResult()
    #plotNeuronResult()
    #boxPlotnumParticles(np.load("numParticles.npy"))
    #data = np.load("../results/final/period_256.npy")
    #boxPlotData4(data, "boxperiod")
    #errors, stderrs = getErrors(10, data, "ratio")
    #plotData(errors, stderrs, "std")

