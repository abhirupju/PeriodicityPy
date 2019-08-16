from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import sys
import subprocess

args = sys.argv[1:]

# Make a `cleanall` rule to get rid of intermediate and library files
if "clean" in args:
    print "Deleting cython files..."
    # Just in case the build directory was created by accident,
    # note that shell=True should be OK here because the command is constant.
    subprocess.Popen("rm -rf build", shell=True, executable="/bin/bash")
    subprocess.Popen("find . -name '*.pyc' -exec rm -rf {} \;", shell=True, executable="/bin/bash")
    subprocess.Popen("find . -name '*.c' -exec rm -rf {} \;", shell=True, executable="/bin/bash")
    subprocess.Popen("find . -name '*.so' -exec rm -rf {} \;", shell=True, executable="/bin/bash")
    # Now do a normal clean
    sys.argv[1] = "clean"
    exit(0)

# We want to always use build_ext --inplace
if args.count("build_ext") > 0 and args.count("--inplace") == 0:
    sys.argv.insert(sys.argv.index("build_ext")+1, "--inplace")

extensions=[
    Extension("main",                                      ["main.pyx"]),
    Extension("algo.autocorrelation",                      ["algo/autocorrelation.pyx"]),
    Extension("algo.fftperiod",                            ["algo/fftperiod.pyx"]),
    Extension("algo.histogram",                            ["algo/histogram.pyx"]),
    Extension("algo.ParticleFilterBuildRelease",           ["algo/ParticleFilterBuildRelease.pyx"]),
    Extension("algo.pdist",           ["algo/pdist.pyx"]),
    Extension("algo.fftautocorrmixed",                     ["algo/fftautocorrmixed.pyx"]),
    Extension("algo.segment_periodicity_gcd",              ["algo/segment_periodicity_gcd.pyx"]),
    Extension("algo.segment_periodicity_incomplete_paper", ["algo/segment_periodicity_incomplete_paper.pyx"]),
    Extension("exp.bgNoiseLessTest",                       ["exp/bgNoiseLessTest.pyx"]),
    Extension("exp.bgNoiseTest",                           ["exp/bgNoiseTest.pyx"]),
    Extension("gen.buildrelease",                          ["gen/buildrelease.pyx"]),
    Extension("gen.lockstep",                              ["gen/lockstep.pyx"]),
    Extension("utils.util",                                ["utils/util.pyx"]),
    Extension("utils.util",                                ["utils/readconfig.pyx"]),
]

setup(
  name = 'PeriodicityPy',
  cmdclass = {'build_ext': build_ext},
  ext_modules = cythonize(extensions),
)
