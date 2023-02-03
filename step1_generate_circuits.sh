#!/bin/bash

set -e

PYTHONPATH=src parallel --ungroup tools/gen_circuits \
    --out_dir out/circuits \
    --distance 3 5 7 9 11 13 15 \
    --noise_model SI1000 \
    --noise_strength {} \
    --style 4-CZ 3-CZ 4-ISWAP 3-ISWAP 3-CZ-wiggle 3-ISWAP-wiggle WIGGLING-CZ GLIDING-CZ SLIDING-CZ 3-CZ_MZZ \
    --basis X Z \
    ::: 0.0001 0.0002 0.0003 0.0005 0.0008 0.001 0.002 0.003 0.004 0.005 0.008 0.01

PYTHONPATH=src parallel --ungroup tools/gen_circuits \
    --out_dir out/circuits \
    --distance 3 5 7 9 11 13 15 \
    --noise_model UniformDepolarizing \
    --noise_strength {} \
    --style 4-CX 3-CX 3-CXSWAP 4-CXSWAP 3-CX-wiggle 3-CXSWAP-wiggle WIGGLING-CX GLIDING-CX SLIDING-CX 3-CX_MXX_MZZ \
    --basis X Z \
    ::: 0.0001 0.0002 0.0003 0.0005 0.0008 0.001 0.002 0.003 0.004 0.005 0.008 0.01

PYTHONPATH=src parallel --ungroup tools/gen_circuits \
    --out_dir out/circuits \
    --distance 4 6 8 10 12 14\
    --noise_model UniformDepolarizing \
    --noise_strength {} \
    --style TORIC-4-CX TORIC-3-CX_MXX_MZZ TORIC-3_HEAVY-CX TORIC-3_SEMI_HEAVY-CX \
    --basis X Z \
    ::: 0.0001 0.0002 0.0003 0.0005 0.0008 0.001 0.002 0.003 0.004 0.005 0.008 0.01
