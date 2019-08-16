from exp import bgNoiseLessTest
from exp import bgNoiseTest

#matplotlib.rcParams['legend.numpoints'] = 1
#matplotlib.rcParams['lines.linewidth'] = 3
#matplotlib.rcParams['font.size'] = 30

def main(runType,  nParticles, prange, scheme):
    #readconfig.init("config/noiseparams")
    t = bgNoiseTest.experiment(runType, nParticles, prange, scheme)
    #pfbr.NoiseLessExperiment()
    #pf.NoiseLessExperiment()