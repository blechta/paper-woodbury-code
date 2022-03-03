#!/bin/bash

MEMUSG_LOG=checkerboard-memusg.log
GNUPLOT_CMD="\
    set xdata time; \
    set timefmt '%s'; \
    set format x '%H:%M'; \
    p '$MEMUSG_LOG' u 4:6 w l; \
"
gnuplot -p -e "$GNUPLOT_CMD"
