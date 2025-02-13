#!/bin/bash

# This script resizes the condensed phase configurations to create a slab geometry
# for direct coexistence simulations. The initial condensed phase is placed in contact
# with a vacuum region. Example slab size: Lz = 180\sigma (condensed phase occupies 36\sigma), 
# while Lx = Ly = 20\sigma for both condensed phase and slab.

set -e  # exit on first error

seqs_file="seqs-ps.lst"  # file containing sequence IDs

while read -r seqid _; do
    echo "seqID: $seqid"
    cd "$seqid" || { echo "dir $seqid not found"; continue; }

    # Generate resized slab configurations
    python3 ../../get_starting_phase_keyfiles/gen_resize_config.py \
        --config-list final-phase-0.atom --output starting-slab-before-cycling-lz-180 phase.json
    
    # uncomment the line below for generating an additional slab size (e.g., Lz = 252\sigma)
    # python3 ../../get_starting_phase_keyfiles/gen_resize_config.py \
    #    --config-list final-phase-0.atom --output starting-slab-before-cycling-lz-252 phase.json

    cd ..
done < "$seqs_file"
