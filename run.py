import main
from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options] runType",
                          version="%prog 1.0")

parser.add_option("-r", "--run", dest='runType', default="diffdistribution", help="specify noise/std/period")
parser.add_option("-s", "--scheme", dest='scheme', default="Reweight_Noise|Systematic_Resample", help="specify scheme")
parser.add_option("-n", "--num", dest='number', default=256, help="number of particles", type=int)
parser.add_option("-p", "--periodrange", dest='prange', default=100000, help="period range", type=float)
(options, args) = parser.parse_args()

print options.runType, int(options.number), int(options.prange), options.scheme

main.main(options.runType, int(options.number), int(options.prange), options.scheme)
