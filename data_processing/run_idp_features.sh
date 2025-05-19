#!/bin/bash

# Define contact cutoff distances for each model
HP_CUTOFFS=("2.0" "2.5" "3.0" "3.5" "4.0")
IDP_CUTOFFS=("16" "20" "24" "28" "32")
IDP_CUTOFFS=("24")

print_section_header() {
    echo ""
    echo "================================================================================"
    echo "$1"
    echo "================================================================================"
    echo ""
}

# Run for IDP model
print_section_header "Processing IDP model..."
for seqid in {2..75}; do
    for cutoff in "${IDP_CUTOFFS[@]}"; do
        echo "Running: python3 generate_contact_map_descriptors.py --seqid $seqid --model idp --cutoff $cutoff"
        python3 generate_contact_map_descriptors.py --seqid "$seqid" --model idp --cutoff "$cutoff"
    done
done
