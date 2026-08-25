#!/bin/bash
#SBATCH --job-name=slabconv
#SBATCH --account=project_2019020
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=%x-%A_%a.out
#SBATCH --array=0-5
#
# Abaqus C3D4H on the mesh-convergence cell, CSC Roihu.
#
# The companion to elements_ccx/tests/slabconv.sh, which runs the SAME decks
# through CalculiX (plain C3D4 and F-barES-FEM-T4 c=1).  make_slabconv.py emits
# both variants from one geometry -- m_abq.inp with C3D4H in the inclusion,
# m_ccx.inp with C3D4 -- so the only difference between what the two codes see
# is the element keyword.  That is the whole point: it makes C3D4H a fourth arm
# on an identical mesh rather than a reference from another campaign with
# another packing.
#
# Stage the decks first, from the laptop:
#   rsync -a out_slabconv/ roihu:/scratch/project_2019020/slabconv/
#
# then
#   sbatch --array=0-5 hpc/submit_slabconv_abaqus.sh
#
# WORKDIR   where the n*/ directories live (default /scratch/<acct>/slabconv)
# CASE      which sweep directory under WORKDIR (default kg500_one_xsym)
set -euo pipefail

source /usr/share/lmod/lmod/init/bash
export MODULEPATH=${MODULEPATH:-/appl/modulefiles}
module load abaqus/2026

WORKDIR=${WORKDIR:-/scratch/project_2019020/slabconv}
CASE=${CASE:-kg500_one_xsym}
NS=(10 20 30 40 50 60)
n=${NS[${SLURM_ARRAY_TASK_ID:-0}]}

# Abaqus writes a lot of scratch; keep it off the shared filesystem.
export TMPDIR=${LOCAL_SCRATCH:-$TMPDIR}

for st in und drn; do
    d="$WORKDIR/$CASE/n$n/$st"
    [ -f "$d/m_abq.inp" ] || { echo "missing $d/m_abq.inp"; exit 1; }
    cd "$d"
    # scratch= keeps the .odb and temporaries local; interactive=off so the
    # job does not return before the solve finishes and the array element
    # exits while Abaqus is still writing.
    abaqus job=m_abq input=m_abq.inp \
           cpus="${SLURM_CPUS_PER_TASK:-16}" \
           scratch="$TMPDIR" \
           interactive
    echo "== $CASE n=$n $st done: $(ls -la m_abq.dat 2>/dev/null | awk '{print $5}') bytes"
done
