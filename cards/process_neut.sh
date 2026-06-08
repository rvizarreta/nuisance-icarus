#!/bin/bash

# Create output directories
mkdir -p /exp/icarus/data/users/rvizarr/nuisance/NEUT/fhc_Nu14
mkdir -p /exp/icarus/data/users/rvizarr/nuisance/NEUT/fhc_Nu-14

# Process fhc_Nu14
echo "Processing fhc_Nu14 files..."
INPUT_DIR1="/pnfs/icarus/persistent/users/jskim/Generators/NEUT/6.1.3/ICARUS/260510_NEUT_Regen/Target_Ar/fhc_Nu14"
OUTPUT_DIR1="/exp/icarus/data/users/rvizarr/nuisance/NEUT/fhc_Nu14"

for infile in ${INPUT_DIR1}/output_NEUT_*.root; do
    if [[ "$infile" == *".nuisflat.root"* ]]; then
        continue
    fi
    filename=$(basename "$infile")
    number=$(echo "$filename" | sed 's/output_NEUT_\(.*\)\.root/\1/')
    outfile="${OUTPUT_DIR1}/output_NEUT_${number}.nuisflat.root"
    echo "Processing: $infile -> $outfile"
    nuisflat -c GENIE.card -f GenericFlux -i "NEUT:${infile}" -o "${outfile}"
done

# Process fhc_Nu-14
echo "Processing fhc_Nu-14 files..."
INPUT_DIR2="/pnfs/icarus/persistent/users/jskim/Generators/NEUT/6.1.3/ICARUS/260510_NEUT_Regen/Target_Ar/fhc_Nu-14"
OUTPUT_DIR2="/exp/icarus/data/users/rvizarr/nuisance/NEUT/fhc_Nu-14"

for infile in ${INPUT_DIR2}/output_NEUT_*.root; do
    if [[ "$infile" == *".nuisflat.root"* ]]; then
        continue
    fi
    filename=$(basename "$infile")
    number=$(echo "$filename" | sed 's/output_NEUT_\(.*\)\.root/\1/')
    outfile="${OUTPUT_DIR2}/output_NEUT_${number}.nuisflat.root"
    echo "Processing: $infile -> $outfile"
    nuisflat -c GENIE.card -f GenericFlux -i "NEUT:${infile}" -o "${outfile}"
done

echo "Done!"