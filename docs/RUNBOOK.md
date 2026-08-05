# Cluster runbook

How to take a campaign from parameter deck to result table on an HPC cluster.
The examples use Slurm and the scripts in `../hpc/`, which are the real
submissions from one site kept as worked examples — adapt the account,
partition and module lines to your own.

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

If generation itself is run as a Slurm array rather than locally, read
`USER_DOCS.md` §5 first — the per-task seed must be derived from the **global**
row index, or every task packs the same microstructure and the campaign's
replicates are not independent.

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

**Sizing.** Most of these are small jobs. A linear `C3D4` slice (~50k elements)
solves in minutes on 4 cores; a quadratic `C3D10H` bending deck (~78k tets)
wants 8 cores, ~32 GB and a couple of hours. Ask for 8 cores rather than more —
on most billing models that keeps the charge equal to the cores used — and set a
walltime that caps a stuck job rather than one that never trips.

**Bending is the exception, and it is memory-bound rather than core-bound.** The
mesh grows as `(L/d)³` and the direct solver dominates, so requirements climb
steeply with cell size. Keep the solver's out-of-core scratch on **node-local**
storage and never on a shared parallel filesystem — point Abaqus `scratch=` at
the node-local path the scheduler grants, and request enough of it. Before
scaling a bending sweep up, read `USER_DOCS.md` §6: the size effect these solves
were originally chasing turned out not to be microstructural.

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

Regenerate the manuscript figures with

```bash
cd results && PYTHONPATH=../analysis python3 ../analysis/make_rev_figs.py
```

and check that nothing has drifted with the two self-verifying harnesses, which
re-derive every published number from the committed tables and fail loudly on a
mismatch:

```bash
cd results
python3 ../analysis/verify_column.py
python3 ../analysis/verify_sizeeffect.py
```

---

## 5. Housekeeping

**Verify first, then delete. The order is the whole point.**

1. Confirm every result file on the cluster has a local counterpart, and that
   the local copy is **byte-identical** — not merely present. Compare checksums,
   not filenames.
2. Only then remove the bulk:

```bash
cd $WORKDIR
find . -maxdepth 1 -name '*.odb'     -delete
find . -maxdepth 1 -name 'Job-*.inp' -delete
```

Keep the `results_*.csv` on the cluster. They are small and they are the record
of what actually ran.

Three things make the ordering matter rather than merely tidy:

- **A campaign is orders of magnitude larger than its results.** Several hundred
  solves can leave tens of gigabytes of ODBs behind, against a few hundred
  kilobytes of CSV. Cluster scratch is shared and usually purged on a timer, so
  leaving it is antisocial as well as risky.
- **A re-run post-process can leave two files covering the same runs with
  different values.** The later one is normally authoritative, but "normally" is
  not "always" — check timestamps and checksums against what is committed before
  deleting the ODBs that would let you re-extract.
- **Decks are reproducible; ODBs are not cheap to recreate.** The parameter CSV
  plus `SPAX_SEED` regenerates `out_*/`, so there is no need to archive it. But
  if a campaign might need re-extraction with different output fields, pull what
  you need before the ODBs go.

For generation output directories the equivalent is to keep only the decks:

```bash
find out_<campaign> -type f ! -name '*.inp' -delete
```

Check the generation summary before submitting anything: `GENERATION COMPLETE:
N RVEs` and the per-deck node counts. Periodic meshing retries at high inclusion
counts are normal and appear in the log as tracebacks — see `USER_DOCS.md`.

## 6. When a solve fails

Identify genuine failures by the **absence of an `Abaqus exit:` line** in the
job log, not by whether an `.odb` exists. A job killed at walltime leaves a
truncated ODB behind, and the solve array skips any deck whose ODB is present —
so a naive resubmit silently skips exactly the jobs that failed. Delete those
ODBs first, then resubmit.
