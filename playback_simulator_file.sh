#!/bin/bash

/usr/lib/uhd/examples/tx_samples_from_file \
    --freq 1575420000 \
    --rate 10000000 \
    --gain 5 \
    --file "gpssim.bin" \
    --ref external
