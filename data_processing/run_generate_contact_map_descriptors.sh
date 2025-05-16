#!/bin/bash

# Define contact cutoff distances for each model
HP_CUTOFFS=("2.0" "2.5" "3.0" "3.5" "4.0")
IDP_CUTOFFS=("08" "12" "16" "20" "24" "28" "32")

print_section_header() {
    echo ""
    echo "================================================================================"
    echo "$1"
    echo "================================================================================"
    echo ""
}

# Run for HP model with b2_value=400
print_section_header "Processing HP model (b2_value=400)..."
for seqid in {1..75}; do
    for cutoff in "${HP_CUTOFFS[@]}"; do
        echo "Running: python3 generate_contact_map_descriptors.py --seqid $seqid --model hp --b2_value 400 --cutoff $cutoff"
        python3 generate_contact_map_descriptors.py --seqid "$seqid" --model hp --b2_value 400 --cutoff "$cutoff"
    done
done

# Run for HP model with b2_value=1000
print_section_header "Processing HP model (b2_value=1000)..."
for seqid in {1..72}; do
    for cutoff in "${HP_CUTOFFS[@]}"; do
        echo "Running: python3 generate_contact_map_descriptors.py --seqid $seqid --model hp --b2_value 1000 --cutoff $cutoff"
        python3 generate_contact_map_descriptors.py --seqid "$seqid" --model hp --b2_value 1000 --cutoff "$cutoff"
    done
done

# Run for IDP model
print_section_header "Processing IDP model..."
for seqid in {1..75}; do
    for cutoff in "${IDP_CUTOFFS[@]}"; do
        echo "Running: python3 generate_contact_map_descriptors.py --seqid $seqid --model idp --cutoff $cutoff"
        python3 generate_contact_map_descriptors.py --seqid "$seqid" --model idp --cutoff "$cutoff"
    done
done
