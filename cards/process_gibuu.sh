#!/bin/bash

# Script to process GiBUU files through nuisflat

# Create output directories if they don't exist
mkdir -p /exp/icarus/data/users/rvizarr/nuisance/GiBUU/fhc_Nu14
mkdir -p /exp/icarus/data/users/rvizarr/nuisance/GiBUU/fhc_Nu-14

# Process fhc_Nu14 directory
echo "Processing fhc_Nu14 files..."
INPUT_DIR1="/exp/icarus/data/users/jskim/Generators/GiBUU/r2025_02/ICARUS_NuMI/250728_GiBUU_MoreStatsFromEAF/Target_Ar/fhc_Nu14"
OUTPUT_DIR1="/exp/icarus/data/users/rvizarr/nuisance/GiBUU/fhc_Nu14"

for infile in ${INPUT_DIR1}/output_GiBUU_*.root; do
    # Skip .nuisflat.root files
    if [[ "$infile" == *".nuisflat.root"* ]]; then
        continue
    fi

    # Extract the number from the filename (e.g., output_GiBUU_1000.root -> 1000)
    filename=$(basename "$infile")
    number=$(echo "$filename" | sed 's/output_GiBUU_\(.*\)\.root/\1/')

    outfile="${OUTPUT_DIR1}/output_GiBUU_${number}.nuisflat.root"

    echo "Processing: $infile -> $outfile"
    nuisflat -c GENIE.card \
        -f GenericFlux \
        -i "GiBUU:${infile}" \
        -o "${outfile}"
done

# Process fhc_Nu-14 directory
echo "Processing fhc_Nu-14 files..."
INPUT_DIR2="/exp/icarus/data/users/jskim/Generators/GiBUU/r2025_02/ICARUS_NuMI/250728_GiBUU_MoreStatsFromEAF/Target_Ar/fhc_Nu-14"
OUTPUT_DIR2="/exp/icarus/data/users/rvizarr/nuisance/GiBUU/fhc_Nu-14"

for infile in ${INPUT_DIR2}/output_GiBUU_*.root; do
    # Skip .nuisflat.root files
    if [[ "$infile" == *".nuisflat.root"* ]]; then
        continue
    fi

    # Extract the number from the filename
    filename=$(basename "$infile")
    number=$(echo "$filename" | sed 's/output_GiBUU_\(.*\)\.root/\1/')

    outfile="${OUTPUT_DIR2}/output_GiBUU_${number}.nuisflat.root"

    echo "Processing: $infile -> $outfile"
    nuisflat -c GENIE.card \
        -f GenericFlux \
        -i "GiBUU:${infile}" \
        -o "${outfile}"
done

echo "Done!"