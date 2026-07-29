#!/bin/bash -l
#SBATCH --job-name=bnd_solve
#SBATCH --output=logs/bnd_solve_%A_%a.out
#SBATCH --error=logs/bnd_solve_%A_%a.err
#SBATCH --account=project_XXXXXX
#SBATCH --cpus-per-task=16
#SBATCH --mem=1000G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=36:00:00
#SBATCH --partition=hugemem
#SBATCH --gres=nvme:1500
# ============ RESOURCES — READ THIS (Puhti) ============
# Bending is RAM-bound (Lesicar constraints + direct solver). L/d=12 OOM'd at
# 40 GB. Rough need: L/d=12 ~150 GB, L/d=14 ~250 GB, L/d=16 ~400 GB+.
# Puhti has NO "bigmem": standard 384 GB nodes are --partition=large --mem=360G
# (fits L/d=12, maybe 14); the BigMem nodes are --partition=hugemem (needed for
# L/d=16, safe for all). To save the scarce hugemem allocation, run sizes
# SEPARATELY: filter GlobalJobList_bend by size (grep L096 / L112 / L128) and put
# L/d=12(/14) on --partition=large --mem=360G, only L/d=16 on hugemem.
# --gres=nvme:1500 -> Abaqus scratch on node-local NVMe ($LOCAL_SCRATCH); the
# out-of-core solver files are huge and must NOT sit on Lustre /scratch.
# =======================================================
# NOTE: --array is set on the command line by bend_01_generate.sh.
unset SLURM_GTIDS
module load abaqus/2025
WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
cd "$WORKDIR" || exit 1

JOBNAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" GlobalJobList_bend)
[ -z "$JOBNAME" ] && { echo "ERROR: no job for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }
[ -f "${JOBNAME}.odb" ] && { echo "SKIP: ${JOBNAME}.odb exists"; exit 0; }
[ ! -f "${JOBNAME}.inp" ] && { echo "ERROR: ${JOBNAME}.inp missing"; exit 1; }

echo "===== ${JOBNAME}  cpus=${SLURM_CPUS_PER_TASK}  start $(date) ====="
# Abaqus scratch on fast node-local NVMe (fall back to a WORKDIR dir if NVMe
# wasn't granted, so the script still runs without --gres).
SCR="${LOCAL_SCRATCH:-$WORKDIR/scratch_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}}"
mkdir -p "$SCR"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

abaqus job="${JOBNAME}" cpus=$SLURM_CPUS_PER_TASK scratch="$SCR" \
       ask_delete=OFF mp_mode=threads memory="90%" interactive
RC=$?
echo "Abaqus exit: ${RC}"

# Keep only THIS job's deck + ODB; drop its .dat/.sta/.msg/.com/.prt/.sim/.lck.
for f in "${JOBNAME}".*; do
    [ -e "$f" ] || continue
    case "${f##*.}" in inp|odb) ;; *) rm -f "$f" 2>/dev/null ;; esac
done
[ "$SCR" != "$LOCAL_SCRATCH" ] && rm -rf "$SCR"
echo "===== ${JOBNAME}  done $(date)  (odb: $([ -f ${JOBNAME}.odb ] && echo yes || echo NO)) ====="
exit $RC
