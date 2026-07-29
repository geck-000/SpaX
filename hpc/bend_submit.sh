#!/bin/bash -l
# Submit the extended bending-only study from the SHARED test_rve dir.
# It scopes all file ops to this study's BND_L* IDs (job list GlobalJobList_bend,
# post_parts_bend, results_bending.csv), so it coexists with other runs there.
# Run on the login node:  bash bend_submit.sh
set -e
WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
mkdir -p "$WORKDIR/logs"
cd "$WORKDIR"
sbatch bend_01_generate.sh
