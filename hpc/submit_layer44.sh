#!/bin/bash -l
# Section 4.4 in layered form: the size sweep of 4.4.1 and the nonlinear
# reload of 4.4.2, both currently resting on pocket-and-channel cells.
#
# Three decks, generated together and solved together, but postprocessed apart
# because the size sweep and the nonlinear reload extract different things.
# The homogeneous control has to be solved at the same sizes and the same
# element size as the layered cells or it measures nothing, which is why it
# rides in the same generation array rather than being run separately.
set -e
WORKDIR=/scratch/project_2019020/test_rve
cd "$WORKDIR"; mkdir -p logs

ACCT=project_2019020
MAXARR=${SPAX_MAX_ARRAY:-25}
# generate_array.sh defaults this to a placeholder path; numpy and gmsh live here.
export PYTHONUSERBASE=/projappl/project_2019020/spax_py

gen () {   # $1 deck  $2 outdir  $3 mesh order
  local N=$(($(wc -l < "$1") - 1))
  sbatch --parsable --account=$ACCT \
    --array=1-${N}%${MAXARR} \
    --export=ALL,WORKDIR=$WORKDIR,CSV=$1,OUTDIR=$2,SPAX_MESH_ORDER=$3,PYTHONUSERBASE=$PYTHONUSERBASE \
    generate_array.sh
}

# Quadratic elements for the bending sweep, matching the pocket sweep it has to
# be compared against; linear for the nonlinear reload, matching that one.
G1=$(gen rve_eringen_layer.csv       out_erglayer  2)
G2=$(gen rve_eringen_layer_homog.csv out_erglayerh 2)
G3=$(gen rve_nlgeom_layer.csv        out_nlglayer  1)
echo "gen: $G1 (15 layered) $G2 (5 control) $G3 (9 nlgeom)"

cat > collect_layer44.sh <<'EOS'
#!/bin/bash -l
#SBATCH --job-name=l44_collect
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=small
#SBATCH --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:30:00
set -e
WORKDIR=/scratch/project_2019020/test_rve
cd "$WORKDIR"
for d in out_erglayer out_erglayerh out_nlglayer; do
  find $d -name 'Job-*.inp' -exec mv -t "$WORKDIR" {} + 2>/dev/null || true
  find $d -name '*_periodic_pairs.csv' -exec mv -t "$WORKDIR" {} + 2>/dev/null || true
done

ls Job-ERGL_*.inp Job-ERGLH_*.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_erglayer
ls Job-NLGL_*.inp                 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_nlglayer
NE=$(wc -l < GlobalJobList_erglayer); NN=$(wc -l < GlobalJobList_nlglayer)
echo "size sweep: $NE decks (expect 20)   nlgeom: $NN decks (expect 9)"
[ "$NE" -ge 1 ] || { echo "ERROR: no ERGL decks"; exit 1; }
[ "$NN" -ge 1 ] || { echo "ERROR: no NLGL decks"; exit 1; }

# The L=0.72 cells are 72 elements to an edge at quadratic order, so the sweep
# gets the memory and the wall time its largest member needs, not its median.
SE=$(sbatch --parsable --account=project_2019020 \
  --partition=small --cpus-per-task=4 --mem=32G --time=01:30:00 \
  --array=1-${NE}%20 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_erglayer \
  csc_solve_array.sh)
echo "solve size sweep: $SE"

SN=$(sbatch --parsable --account=project_2019020 \
  --partition=small --cpus-per-task=4 --mem=16G --time=01:00:00 \
  --array=1-${NN}%9 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_nlglayer \
  csc_solve_array.sh)
echo "solve nlgeom: $SN"

PE=$(sbatch --parsable --account=project_2019020 --dependency=afterany:${SE} \
  --export=ALL,WORKDIR=$WORKDIR,CSV=rve_eringen_layer.csv,RESULTS=results_eringen_layer.csv \
  postprocess_firstorder.sh)
PH=$(sbatch --parsable --account=project_2019020 --dependency=afterany:${SE} \
  --export=ALL,WORKDIR=$WORKDIR,CSV=rve_eringen_layer_homog.csv,RESULTS=results_eringen_layer_homog.csv \
  postprocess_firstorder.sh)
echo "post sweep: $PE / $PH"

PN=$(sbatch --parsable --account=project_2019020 --dependency=afterany:${SN} \
  --export=ALL,WORKDIR=$WORKDIR,CSV=rve_nlgeom_layer.csv,SUMM=results_nlgeom_layer.csv,CURVES=curves_nlgeom_layer.csv \
  postprocess_nlgeom.sh)
echo "post nlgeom: $PN"
EOS
chmod +x collect_layer44.sh

C=$(sbatch --parsable --account=$ACCT \
  --dependency=afterany:${G1}:${G2}:${G3} collect_layer44.sh)
echo "collect+solve+post chained as $C"
