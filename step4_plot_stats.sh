#!/bin/bash

set -e

mkdir -p out/assets/regen/footprint
mkdir -p out/assets/regen/plot

PYTHONPATH=src ./tools/gen_figures_latex > out/assets/regen/circuit_figures.tex

# gen circuit schedule diagrams & convert to png
PYTHONPATH=src ./tools/gen_schedule_figures
ls out/assets/regen/schedules | grep -P "\.svg$" | parallel inkscape -o out/assets/regen/schedules/{}.png -w 1024 -h 1024 out/assets/regen/schedules/{}


PYTHONPATH=src parallel tools/plot_fan \
    out/fused_stats.csv \
    --unit quop \
    --basis XZ \
    --style {} \
    --out out/assets/regen/plot/fan_{}.png \
    ::: TORIC-3-CX_MXX_MZZ


#make extended benchmarking plots for each circuit
PYTHONPATH=src parallel tools/plot_fan \
    out/fused_stats.csv \
    --unit quop \
    --basis XZ \
    --style {} \
    --out out/assets/regen/plot/fan_{}.png \
    ::: 4-CX 4-CZ 3-CX 3-CZ 3-CX-wiggle 3-CZ-wiggle 4-CXSWAP 4-ISWAP 3-CXSWAP 3-CXSWAP-wiggle 3-ISWAP-wiggle WIGGLING-CX WIGGLING-CZ GLIDING-CX GLIDING-CZ SLIDING-CX SLIDING-CZ TORIC-4-CX TORIC-3_HEAVY-CX TORIC-3_SEMI_HEAVY-CX 3-CX_MXX_MZZ 3-CZ_MZZ TORIC-3-CX_MXX_MZZ

PYTHONPATH=src parallel tools/plot_classic \
    out/fused_stats.csv \
    --unit round \
    --basis XZ \
    --style {} \
    --out out/assets/regen/plot/classic_{}.png \
    ::: 4-CX 4-CZ 3-CX 3-CZ 3-CX-wiggle 3-CZ-wiggle 4-CXSWAP 4-ISWAP 3-CXSWAP 3-CXSWAP-wiggle 3-ISWAP-wiggle WIGGLING-CX WIGGLING-CZ GLIDING-CX GLIDING-CZ SLIDING-CX SLIDING-CZ TORIC-4-CX TORIC-3_HEAVY-CX TORIC-3_SEMI_HEAVY-CX 3-CX_MXX_MZZ 3-CZ_MZZ TORIC-3-CX_MXX_MZZ


# make footprint plots
PYTHONPATH=src tools/plot_footprint \
    out/fused_stats.csv \
    --unit teraquop \
    --basis XZ \
    --filter_func "metadata['style'] in ['4-CZ', '3-CZ']" \
    --order_func "['4-CZ', '3-CZ'].index(metadata['style'])" \
    --label_func "metadata['style']" \
    --title "Teraquop footprints for standard and hex-grid circuits" \
    --out out/assets/regen/footprint/footprint_hex.png

PYTHONPATH=src tools/plot_footprint \
    out/fused_stats.csv \
    --unit teraquop \
    --basis XZ \
    --filter_func "metadata['style'] in ['4-CZ', '4-ISWAP']" \
    --order_func "['4-CZ', '4-ISWAP'].index(metadata['style'])" \
    --label_func "metadata['style']" \
    --title "Teraquop footprints for standard and ISWAP circuits" \
    --out out/assets/regen/footprint/footprint_iswap.png

PYTHONPATH=src tools/plot_footprint \
    out/fused_stats.csv \
    --unit teraquop \
    --basis XZ \
    --filter_func "metadata['style'] in ['4-CZ', 'WIGGLING-CZ', 'GLIDING-CZ',  'SLIDING-CZ']" \
    --order_func "['4-CZ', 'WIGGLING-CZ', 'GLIDING-CZ',  'SLIDING-CZ'].index(metadata['style'])" \
    --label_func "metadata['style']" \
    --title "Teraquop footprints for standard and walking circuits" \
    --out out/assets/regen/footprint/footprint_walking.png

PYTHONPATH=src tools/plot_footprint \
    out/fused_stats.csv \
    --unit teraquop \
    --basis XZ \
    --filter_func "metadata['style'] in ['4-CZ', '3-CZ', '3-CZ-wiggle', '4-ISWAP', '3-ISWAP', '3-ISWAP-wiggle', 'WIGGLING-CZ', 'GLIDING-CZ', 'SLIDING-CZ']" \
    --order_func "['4-CZ', '3-CZ', '3-CZ-wiggle', '4-ISWAP', '3-ISWAP', '3-ISWAP-wiggle', 'WIGGLING-CZ', 'GLIDING-CZ', 'SLIDING-CZ', '3-CZ_MZZ'].index(metadata['style'])" \
    --label_func "metadata['style']" \
    --title "Teraquop footprints for all planar CZ and ISWAP circuits" \
    --out out/assets/regen/footprint/footprint_si1000.png

PYTHONPATH=src tools/plot_footprint \
    out/fused_stats.csv \
    --unit teraquop \
    --basis XZ \
    --filter_func "metadata['style'] in ['4-CX', '3-CX', '4-CXSWAP', '3-CXSWAP', 'WIGGLING-CX', '3-CX-wiggle', '3-CXSWAP-wiggle']" \
    --order_func "['4-CX', '3-CX', '4-CXSWAP', '3-CXSWAP', 'WIGGLING-CX', '3-CX-wiggle', '3-CXSWAP-wiggle', '3-CX_MXX_MZZ'].index(metadata['style'])" \
    --label_func "metadata['style']" \
    --title "Teraquop footprints for all planar CX and CXSWAP circuits" \
    --out out/assets/regen/footprint/footprint_uniform_depolarizing.png

PYTHONPATH=src tools/plot_footprint \
    out/fused_stats.csv \
    --unit teraquop \
    --basis XZ \
    --filter_func "metadata['style'] in ['TORIC-4-CX', 'TORIC-3_HEAVY-CX', 'TORIC-3_SEMI_HEAVY-CX']" \
    --order_func "['TORIC-4-CX', 'TORIC-3_HEAVY-CX', 'TORIC-3_SEMI_HEAVY-CX'].index(metadata['style'])" \
    --label_func "metadata['style']" \
    --title "Teraquop footprints for all toric circuits" \
    --out out/assets/regen/footprint/footprint_toric.png

PYTHONPATH=src tools/plot_footprint \
    out/fused_stats.csv \
    --unit teraquop \
    --basis XZ \
    --filter_func "metadata['style'] in ['3-CZ_MZZ', '3-CX_MXX_MZZ', 'TORIC-3-CX_MXX_MZZ']" \
    --order_func "['3-CZ_MZZ', '3-CX_MXX_MZZ', 'TORIC-3-CX_MXX_MZZ'].index(metadata['style'])" \
    --label_func "metadata['style']" \
    --title "Teraquop footprints for all hybrid entanglement circuits" \
    --out out/assets/regen/footprint/footprint_hybrid.png
