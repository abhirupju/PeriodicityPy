#!/bin/bash

python run.py -r noise -n 1000 -p 30 -s Dynamics_Mixed_PTwoSided_Dynamics_MW_W
python run.py -r noise -n 1000 -p 30 -s Dynamics_Mixed_PTwoSided_Resample_ESS
python run.py -r noise -n 1000 -p 30 -s Dynamics_Mixed_PTwoSided_Resample_ESS_Dynamics_MW_W

python run.py -r noise -n 1000 -p 30 -s Dynamics_ESS_Dynamics_MW_W
python run.py -r noise -n 1000 -p 30 -s Dynamics_ESS_Resample_ESS
python run.py -r noise -n 1000 -p 30 -s Dynamics_ESS_Resample_ESS_Dynamics_MW_W

