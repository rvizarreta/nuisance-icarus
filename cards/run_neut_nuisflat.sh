#!/bin/bash
# HTCondor wrapper for NEUT nuisflat
# Usage: run_neut_nuisflat.sh <flux_dir> <file_number>
# Example: run_neut_nuisflat.sh fhc_Nu14 0

FLUX_DIR=$1    # fhc_Nu14 or fhc_Nu-14
FILE_NUM=$2    # 0, 1, 2, ...

INPUT_BASE="/pnfs/icarus/persistent/users/jskim/Generators/NEUT/6.1.3/ICARUS/260510_NEUT_Regen/Target_Ar"
OUTPUT_BASE="/exp/icarus/data/users/rvizarr/nuisance/NEUT"
CARDS_DIR="/exp/icarus/app/users/rvizarr/nuisance-icarus/nuisance-icarus/cards"

INFILE="${INPUT_BASE}/${FLUX_DIR}/output_NEUT_${FILE_NUM}.root"
OUTFILE="${OUTPUT_BASE}/${FLUX_DIR}/output_NEUT_${FILE_NUM}.nuisflat.root"

# Set up environment
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup genie v3_04_00 -q e20:prof

# Source NEUT
source /exp/icarus/app/users/rvizarr/nuisance-icarus/neut/build/Linux/bin/setup.NEUT.sh

# Source nuisance-icarus
source /exp/icarus/app/users/rvizarr/nuisance-icarus/build/Linux/setup.sh

# Run
cd ${CARDS_DIR}
echo "Processing: ${INFILE} -> ${OUTFILE}"
nuisflat -c GENIE.card \
    -f GenericFlux \
    -i "NEUT:${INFILE}" \
    -o "${OUTFILE}"

echo "Exit code: $?"