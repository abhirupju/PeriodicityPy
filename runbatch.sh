#!/bin/bash

# Kill all jombies if things go wrong:
# kill -9 `ps aux | grep s1563028 | awk '{if ($11~"python")print $2}'`

echo "
#-> Plain
#    -Period range: 0-30
#    -Prior: period - uniform
#    -Resample: every step
#    -Dynamics: static. t-dist scale=0.1
"
python run.py -r noise -n 1000 -p 30 &
#python run.py -r noise -n 5000 -p 100000 &
python run.py -r std -n 1000 -p 30 &
#python run.py -r std -n 5000 -p 100000 &
wait
echo "
#-> Resample_ESS
#    -Period range: 0-30
#    -Prior: period - uniform
#    -Resample: ess based N/2
#    -Dynamics: static. t-dist scale=0.1
"
python run.py -r noise -n 1000 -p 30 -s Resample_ESS &
#python run.py -r noise -n 5000 -p 100000 -s Resample_ESS &
python run.py -r std -n 1000 -p 30 -s Resample_ESS &
#python run.py -r std -n 5000 -p 100000 -s Resample_ESS &
wait
echo "
#-> Dynamics_Mixed_PHalf
#    -Period range: 0-30
#    -Prior: period - uniform
#    -Resample: every step
#    -Dynamics: static. Mixed dist. pi*(oldP/2) + (1-pi)*(t-dist scale=0.1)
"
python run.py -r noise -n 1000 -p 30 -s Dynamics_Mixed_PHalf &
#python run.py -r noise -n 5000 -p 100000 -s Dynamics_Mixed_PHalf &
python run.py -r std -n 1000 -p 30 -s Dynamics_Mixed_PHalf &
#python run.py -r std -n 5000 -p 100000 -s Dynamics_Mixed_PHalf &
wait
echo "
#-> Dynamics_Mixed_PHalf_PDouble
#    -Period range: 0-30
#    -Prior: period - uniform
#    -Resample: every step
#    -Dynamics: static. Mixed dist. pi1*(oldP/2) + pi2*(oldP*2) * (1-pi1-pi2)*(t-dist scale=0.1)
"
python run.py -r noise -n 1000 -p 30 -s Dynamics_Mixed_PTwoSided &
#python run.py -r noise -n 5000 -p 100000 -s Dynamics_Mixed_PHalf_PDouble &
python run.py -r std -n 1000 -p 30 -s Dynamics_Mixed_PTwoSided &
#python run.py -r std -n 5000 -p 100000 -s Dynamics_Mixed_PHalf_PDouble &
wait
echo "
#-> Dynamics_ESS
#    -Period range: 0-30
#    -Prior: period - uniform
#    -Resample: every step
#    -Dynamics: static. t-dist scale=1/ESS
"
python run.py -r noise -n 1000 -p 30 -s Dynamics_ESS &
#python run.py -r noise -n 5000 -p 100000 -s Dynamics_ESS &
python run.py -r std -n 1000 -p 30 -s Dynamics_ESS &
#python run.py -r std -n 5000 -p 100000 -s Dynamics_ESS &
wait
echo "
#-> Dynamics_MW_W
#    -Period range: 0-30
#    -Prior: period - uniform
#    -Resample: every step
#    -Dynamics: static. dynamics only particles with (w < max_weight) with t-dist scale=0.1
"
python run.py -r noise -n 1000 -p 30 -s Dynamics_MW_W &
#python run.py -r noise -n 5000 -p 100000 -s Dynamics_MW_W &
python run.py -r std -n 1000 -p 30 -s Dynamics_MW_W &
#python run.py -r std -n 5000 -p 100000 -s Dynamics_MW_W &
wait

#nParticle: 1000
#nSamples: 100
#nIterations: 10
#Period: 10

#Test Each with:
#1. Prior_Prange_uniform: Uniform with range 0-30
#2. Prior_Punbound_expon: expon scale=1, range 0-1e10
#3. NParticle_5000: 5000

#** Save configuration and numpy result with scheme name

