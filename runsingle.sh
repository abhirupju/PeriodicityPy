#!/bin/bash

#nice -5 python run.py -r std -n 256 -p 1000000 -s Reweight_Noise_Systematic_Resample
#nice -5 python run.py -r noise -n 256 -p 1000000 -s Reweight_Noise_Systematic_Resample
#nice -5 python run.py -r period -n 256 -p 1000000 -s Reweight_Noise_Systematic_Resample


#nice -5 python run.py -r numParticles -n 500 -p 1000000 -s Reweight_Noise_Systematic_Resample


#nice -5 python run.py -r online -n 256 -p 1000000 -s Reweight_Noise_Systematic_Resample

python2 run.py

#nice -5 python run.py -r cputime -n 256 -p 1000000 -s Reweight_Noise_Systematic_Resample
