#!/bin/bash

set -e

PYTHONPATH=src tools/fuse_xz_data \
    --stats out/stats.csv \
    > out/fused_stats.csv
