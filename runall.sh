#!/bin/bash
set -e
NUM_THREADS=8
MATLAB_CMD="\
    maxNumCompThreads($NUM_THREADS); \
    app_dc.checkerboard_example.full_benchmark; \
    app_dc.checkerboard_example.postprocess; \
    app_dc.checkerboard_example.postprocess2; \
    app_dc.checkerboard_example.postprocess3; \
    app_dc.checkerboard_example.postprocess4; \
    app_dc.checkerboard_example.plot_configs; \
"
set -x
nohup matlab -nodisplay -logfile checkerboard-matlab.log -batch "$MATLAB_CMD" > /dev/null &
MATLAB_PID=$!
nohup ./memusg $MATLAB_PID > checkerboard-memusg.log &
