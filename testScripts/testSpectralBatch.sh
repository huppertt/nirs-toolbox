#!/bin/sh

module load cuda42

/home/pkg/`hostname`/MATLAB/R2012b/bin/matlab -nodesktop < /home/barker/PhD_Data/+nirs/testScripts/testSpectral.m &

wait

exit