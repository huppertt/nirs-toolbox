#!/bin/sh

module load cuda42

/home/pkg/`hostname`/MATLAB/R2014a/bin/matlab -nodesktop < /home/barker/MATLAB/+nirs/demo/sphere_demo.m &

wait

exit