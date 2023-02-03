#!/bin/bash

set -e

sinter collect \
    --circuits out/circuits/*.stim \
    --save_resume_filepath out/stats.csv \
    --metadata_func auto \
    --decoders pymatching \
    --max_shots 1_000_000 \
    --max_errors 1000 \
    --processes 4
