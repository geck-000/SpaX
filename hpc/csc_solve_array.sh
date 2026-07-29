#!/bin/bash -l
#SBATCH --job-name=csc_solve
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --account=project_XXXXXX
#SBATCH --ntasks=1
#SBATCH --nodes=1
# Generic single-deck Abaqus solver for a Slurm array. The wrapper that submits
# this (submit_si2nd_l400.sh / submit_failure.sh) sets --partition, --cpus-per-task,
# --mem, --time, --array and exports WORKDIR + JOBLIST. One array task = one deck.
unset SLURM_GTIDS
# CSC's /etc/profile.d/zz-csc-env.sh bails out for non-interactive shells (leaving
# `module` as the Tcl env-modules that CANNOT parse Abaqus's .lua Lmod file ->
# "Magic cookie missing"). Requesting non-interactive init sources Lmod + StdEnv
# so `module load abaqus/2026` works regardless of how the job was submitted.
export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load abaqus/2026
WORKDIR=${WORKDIR:?set WORKDIR}
JOBLIST=${JOBLIST:?set JOBLIST}
cd "$WORKDIR" || exit 1

JOBNAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
[ -z "$JOBNAME" ] && { echo "ERROR: no job for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }
[ -f "${JOBNAME}.odb" ] && { echo "SKIP: ${JOBNAME}.odb exists"; exit 0; }
[ ! -f "${JOBNAME}.inp" ] && { echo "ERROR: ${JOBNAME}.inp missing"; exit 1; }

echo "===== ${JOBNAME}  cpus=${SLURM_CPUS_PER_TASK}  start $(date) ====="
SCR="${LOCAL_SCRATCH:-$WORKDIR/scratch_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}}"
mkdir -p "$SCR"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

abaqus job="${JOBNAME}" cpus=$SLURM_CPUS_PER_TASK scratch="$SCR" \
       ask_delete=OFF mp_mode=threads memory="90%" interactive
RC=$?
echo "Abaqus exit: ${RC}"

# Keep only this job's deck + ODB; drop .dat/.sta/.msg/.com/.prt/.sim/.lck.
for f in "${JOBNAME}".*; do
    [ -e "$f" ] || continue
    case "${f##*.}" in inp|odb) ;; *) rm -f "$f" 2>/dev/null ;; esac
done
[ "$SCR" != "$LOCAL_SCRATCH" ] && rm -rf "$SCR"
echo "===== ${JOBNAME}  done $(date)  (odb: $([ -f ${JOBNAME}.odb ] && echo yes || echo NO)) ====="
exit $RC
