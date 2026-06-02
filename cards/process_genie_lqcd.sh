#!/bin/bash

# Script to process GENIE AR23_20i_01_000 (FA from LQCD, arXiv:2601.02676) through nuisflat
# Raw .root files are on /pnfs (tape-backed)
# Path layout: .../AR23_20i_01_000/ICARUS/260123_3_04_00_AR23_FASeries/Target_Ar/fhc_Nu±14/

# Create output directories if they don't exist
mkdir -p /exp/icarus/data/users/rvizarr/nuisance/GENIE_LQCD/fhc_Nu14
mkdir -p /exp/icarus/data/users/rvizarr/nuisance/GENIE_LQCD/fhc_Nu-14

# Process fhc_Nu14 directory
echo "Processing fhc_Nu14 files..."
INPUT_DIR1="/pnfs/icarus/persistent/users/jskim/Generators/GENIE/AR23_20i_01_000/ICARUS/260123_3_04_00_AR23_FASeries/Target_Ar/fhc_Nu14"
OUTPUT_DIR1="/exp/icarus/data/users/rvizarr/nuisance/GENIE_LQCD/fhc_Nu14"

for infile in ${INPUT_DIR1}/output_GENIE_*.root; do
    # Skip .nuisflat.root files
    if [[ "$infile" == *".nuisflat.root"* ]]; then
        continue
    fi

    # Extract the number from the filename (e.g., output_GENIE_0.root -> 0)
    filename=$(basename "$infile")
    number=$(echo "$filename" | sed 's/output_GENIE_\(.*\)\.root/\1/')

    outfile="${OUTPUT_DIR1}/output_GENIE_${number}.nuisflat.root"

    echo "Processing: $infile -> $outfile"
    nuisflat -c GENIE.card \
        -f GenericFlux \
        -i "GENIE:${infile}" \
        -o "${outfile}"
done

# Process fhc_Nu-14 directory
echo "Processing fhc_Nu-14 files..."
INPUT_DIR2="/pnfs/icarus/persistent/users/jskim/Generators/GENIE/AR23_20i_01_000/ICARUS/260123_3_04_00_AR23_FASeries/Target_Ar/fhc_Nu-14"
OUTPUT_DIR2="/exp/icarus/data/users/rvizarr/nuisance/GENIE_LQCD/fhc_Nu-14"

for infile in ${INPUT_DIR2}/output_GENIE_*.root; do
    # Skip .nuisflat.root files
    if [[ "$infile" == *".nuisflat.root"* ]]; then
        continue
    fi

    # Extract the number from the filename
    filename=$(basename "$infile")
    number=$(echo "$filename" | sed 's/output_GENIE_\(.*\)\.root/\1/')

    outfile="${OUTPUT_DIR2}/output_GENIE_${number}.nuisflat.root"

    echo "Processing: $infile -> $outfile"
    nuisflat -c GENIE.card \
        -f GenericFlux \
        -i "GENIE:${infile}" \
        -o "${outfile}"
done

echo "Done!"