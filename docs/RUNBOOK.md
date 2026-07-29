# Cluster runbook

How to take a campaign from parameter deck to result table on an HPC cluster.
The examples use Slurm and the scripts in `../hpc/`, which were written for CSC
Roihu; adapt the account, partition and module lines for your own site.

For the toolkit itself — parameters, environment variables, output columns —
see the repository `README.md`. For operational pitfalls, see `USER_DOCS.md`.

---

## 0. What runs where

| Stage | Needs Abaqus? | Where |
|-------|---------------|-------|
| Generate decks from a parameter CSV | no | laptop or cluster |
| Solve the decks | **yes** | cluster |
| Extract ODBs → `results_*.csv` | **yes** (`abaqus python`) | cluster |
| Analyse CSVs → figures and numbers | no | laptop |

Only the middle two stages need a licence, so generation and the whole analysis
stage run locally. Generation is the memory-hungry stage; solving is the
wall-clock-hungry one.

---

## 1. Stage the campaign

```bash
# 1. build the parameter deck
cd studies && python3 make_<campaign>.py          # -> ../params/rve_<campaign>.csv

# 2. generate the Abaqus input decks locally
cd .. && SPAX_SEED=<seed> python3 SpaX_Standalone.py \
         params/rve_<campaign>.csv out_<campaign>/

# 3. copy the decks and the pipeline to the cluster
export WORKDIR=/scratch/<project>/<workdir>
rsync -az out_<campaign>/ <cluster>:$WORKDIR/ --include='Job-*' --exclude='*'
rsync -az SpaX_Standalone.py SpaX_PostProcess.py \
          params/rve_<campaign>.csv <cluster>:$WORKDIR/
```

Both `SpaX_Standalone.py` and `SpaX_PostProcess.py` must be present in
`WORKDIR`: the post-processor imports the deck reader from the generator.

## 2. Solve

```bash
ssh <cluster>
cd $WORKDIR
WORKDIR=$WORKDIR bash submit_<campaign>.sh
```

Every script in `../hpc/` honours `WORKDIR` from the environment, so nothing
needs editing. The Slurm account cannot be taken from a variable in a `#SBATCH`
directive — either edit the `--account=` line once, or override per submission:

```bash
sbatch --account=<your_account> submit_<campaign>.sh
```

**Sizing.** These are small jobs. A linear `C3D4` column slice (~50k elements)
solves in minutes on 4 cores; a quadratic `C3D10H` bending deck (~78k tets) wants
8 cores, ~32 GB and a couple of hours. Ask for 8 cores rather than more — on most
billing models that keeps the charge equal to the cores used — and set a
walltime that caps a stuck job rather than one that never trips.

## 3. Post-process

```bash
WORKDIR=$WORKDIR bash postprocess_<campaign>.sh
```

which runs, per deck,

```bash
abaqus python SpaX_PostProcess.py <params.csv> <odb_dir> results_<campaign>.csv
```

Array tasks write `results_<campaign>_<i>.csv` partials; union them with

```bash
python3 SpaX_PostProcess.py --merge parts_dir/ results_<campaign>.csv
```

## 4. Pull back and analyse

```bash
rsync -az <cluster>:$WORKDIR/results_<campaign>.csv results/
cd results && python3 ../analysis/<analyzer>.py
```

Analysers read their inputs by bare filename, so run them from `results/`.

---

## 5. Housekeeping

- **Delete ODBs once the results CSV is pulled.** A campaign is several GB of
  ODBs and the CSV holds everything the analysis needs. Cluster scratch is
  usually shared and often purged on a timer.
- **Decks are reproducible**, so there is no need to archive `out_*/`. The
  parameter CSV plus `SPAX_SEED` regenerates them.
- **Check the generation summary before submitting**: `GENERATION COMPLETE: N
  RVEs` and the per-deck node counts. Periodic meshing retries at high inclusion
  counts are normal and appear in the log as tracebacks — see `USER_DOCS.md`.

## 6. When a solve fails

Identify genuine failures by the **absence of an `Abaqus exit:` line** in the
job log, not by whether an `.odb` exists. A job killed at walltime leaves a
truncated ODB behind, and the solve array skips any deck whose ODB is present —
so a naive resubmit silently skips exactly the jobs that failed. Delete those
ODBs first, then resubmit.
