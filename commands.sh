#!/usr/bin/env bash
set -eou pipefail

python3 scripts/tral_detect.py \
    -f data/primates_SLC6A3_reformat.fa \
    -s DNA \
    -o results/ \
    -p -1 \
    -d XSTREAM TRF T-REKS

python3 scripts/tral_score.py \
    -d results/raw \
    -m phylo_gap01

python3 scripts/tral_filter.py \
    -i results/raw/ \
    -o results/filtered/ \
    -m phylo_gap01 \
    -p 0.05 \
    -d 0.05 \
    -l 7

python3 scripts/tral_refine.py \
    -r results/filtered/ \
    -o results/refined/ \
    -f data/primates_SLC6A3_reformat.fa \
    -s DNA \
    -c phylo_gap01 \
    -p 0.05 \
    -d 0.05


