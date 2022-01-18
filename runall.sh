#!/bin/bash
set -e
MATLAB_CMD="\
    app_dc.checkerboard_example.full_benchmark; \
    app_dc.checkerboard_example.postprocess; \
    app_dc.checkerboard_example.postprocess2; \
    app_dc.checkerboard_example.plot_configs; \
"
set -x
matlab -nodisplay -logfile checkerboard-matlab.log -batch "$MATLAB_CMD"
