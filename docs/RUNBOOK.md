# SpaX — runbook

Everything needed to run a campaign end to end. For what the toolkit is and
why it is built this way, see [`USER_DOCS.md`](USER_DOCS.md).

Cluster instructions here are deliberately generic: they assume Slurm and
nothing else. The scripts in `hpc/` are the real ones from the machine this
work was run on, kept as worked examples of how a campaign was actually
submitted — read them for the shape of a submission, not as something to run
unmodified. Set `WORKDIR`, the Slurm account and the partition names to match
your own site.

---

## 1. Environment

```bash
pip install numpy gmsh          # generation and meshing
pip install matplotlib pandas   # analysis and figures
pip install pyvista             # optional, 3-D renders
```

Abaqus is needed **only** to solve decks and to read ODBs. Generation and the
entire analysis stage run without it and without a licence.

---

## 2. The cycle

Generation → solve → post-process. Each stage writes files the next one reads,
so stages can be run on different machines.

### 2.1 Generate

```bash
cd params && python3 ../studies/make_<campaign>.py     # -> rve_<campaign>.csv
cd ..
SPAX_SEED=<seed> python3 SpaX_Standalone.py params/rve_<campaign>.csv out_<campaign>/
```

This writes one Abaqus input deck per RVE per load case into `out_<campaign>/`.

**Seeding.** `SPAX_SEED` fixes the packing, and the seed is resolved *per row
index*. The same seed with a different number of rows therefore gives different
packings. Two campaigns are comparable packing-for-packing only if their
parameter CSVs have the same shape. When generation is split across parallel
array tasks, the per-task seed must be derived from the task's **global row
index**, not from its index within its own slice — otherwise every task
generates the same packing and "replicates" come out near-identical. See
`hpc/generate_array.sh` for the derivation.

### 2.2 Solve

Locally, for a small campaign:

```bash
abaqus job=Job-<id>-<mode> input=Job-<id>-<mode>.inp cpus=4 interactive
```

On a cluster, submit the decks as a Slurm array. The pattern used throughout
`hpc/` is: stage the decks and the three core modules into `$WORKDIR`, build a
job list, size the array to it, and chain post-processing behind the solve with
a dependency.

```bash
rsync -az out_<campaign>/ <host>:$WORKDIR/ --include='Job-*' --exclude='*'
ssh <host> "cd $WORKDIR && bash submit_<campaign>.sh"
```

### 2.3 Post-process

```bash
abaqus python SpaX_PostProcess.py params/rve_<campaign>.csv out_<campaign>/ results_<campaign>.csv
```

Then pull `results_*.csv` (and `post_*/` if the full tensor was requested) back
to the repository's `results/`.

---

## 3. Reproducing the published figures

All analyzers are plain Python and read the CSVs already committed under
`results/`. Run them **from** `results/`, which is where they look for inputs:

```bash
cd results && PYTHONPATH=../analysis python3 ../analysis/make_rev_figs.py
```

Two self-checking harnesses re-derive every quoted number and fail loudly if
one drifts:

```bash
cd results
python3 ../analysis/verify_column.py       # depth-graded column
python3 ../analysis/verify_sizeeffect.py   # bending size effect
```

---

## 4. Housekeeping — delete decks and ODBs after pulling results

**This is policy, not a suggestion.** Cluster scratch is a shared, finite
resource, and a campaign's ODBs dwarf everything else: a single set of solves
can leave tens of gigabytes behind, while the extracted results are a few
hundred kilobytes of CSV.

Once results are extracted **and pulled back and verified**:

```bash
cd $WORKDIR
find . -maxdepth 1 -name '*.odb'     -delete
find . -maxdepth 1 -name 'Job-*.inp' -delete
```

Keep the `results_*.csv` on the cluster — they are small and are the record of
what ran.

For generation output directories, the equivalent is to keep only the decks:

```bash
find out_<campaign> -type f ! -name '*.inp' -delete
```

Order matters. Verify *before* deleting: confirm every result file on the
cluster has a local counterpart, and that the local copy is byte-identical, not
merely present. A re-run post-process can leave two files covering the same
runs with different values; the later one is normally authoritative, but check
rather than assume.

Decks are cheap to regenerate from `params/`; ODBs are not, so if a campaign
might need re-extraction with different fields, pull what you need first.

---

## 5. Traps

Read these before running anything.

- **A solve array that skips a deck when its `.odb` merely exists** is testing
  for the file, not for a completed solve. A job killed at walltime leaves a
  truncated ODB, so a naive resubmit silently skips exactly the jobs that
  failed. Identify failures by the absence of an `Abaqus exit:` line in the
  log, and delete those ODBs before resubmitting.

- **Strain must be read from the deck, never assumed.** The tensor extractor
  once defaulted to a hard-coded applied strain while the decks prescribe a
  fixed *displacement*. The two coincided exactly at one box size, so it was
  silently correct for every campaign until the cell size changed, then
  silently rescaled the whole tensor. It now reads the deck and raises rather
  than guessing. Ratios were always immune, since a common factor cancels.

- **Generation forks one mesher per core** unless `SPAX_GEN_WORKERS` says
  otherwise. At around 7×10⁵ elements each worker needs ~3 GB, so on a 16 GB
  machine set `SPAX_GEN_WORKERS=2`, and `=1` above 10⁶ elements. Killing the
  parent leaves **orphaned mesher children** still holding the memory — and
  their command line reads `GmshPeriodic`, not `Standalone`, so the obvious
  `pkill` pattern misses them.

- **Periodic meshing retries are normal** at high inclusion counts, where more
  inclusions cross the periodic faces. A traceback in the log is usually the
  retry mechanism working. Check the reported count of completed RVEs and the
  per-deck node counts before concluding the decks are bad.

- **Module systems change under you.** If `module load abaqus` fails inside a
  batch job while working interactively, source a snapshot of the environment
  instead. Copy the working form from a script known to be current rather than
  from the oldest one in the directory.

- **Bending is memory-bound**, far more than the uniaxial cases: the mesh grows
  as `(L/d)³` and the direct solver is the constraint. Keep the solver's
  out-of-core scratch on node-local storage, never on a shared parallel
  filesystem. Before scaling up, read §4.3 of `USER_DOCS.md` — the size effect
  this was chasing turned out to be a kinematic artefact, so larger bending
  cells are unlikely to be worth their cost.

- **Scatter is the population standard deviation** (`ddof=0`) everywhere. The
  `pandas` and `statistics` defaults are `ddof=1`.
