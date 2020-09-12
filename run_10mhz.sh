#!/bin/bash

if [ $# -lt 1 ] ; then
    SIM_TIME="10"
else
    SIM_TIME=$1
fi

rm gpssim.bin
./gnss_sim -e brdc3540.14n \
    -s 10000000 \
    -l 35.681298,139.766247,10.0 \
    -d $SIM_TIME
