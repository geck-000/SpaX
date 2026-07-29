# RVE-size convergence study (statistical representativeness)

Determines how large the periodic RVE must be **relative to the inclusions** for
the homogenised properties to be statistically representative.

## Design

Everything is held fixed except the box size `L`:

- material: ice matrix (E=9.43 GPa, ν=0.33) + **10% void + 10% soft brine** (K/G),
  mild Z-elongated ellipsoids (sphericity 0.75) — a representative sea-ice config
- inclusion size fixed: `r_avg=0.04` (diameter `d=0.08`), `r_std=0.01`
- mesh resolution fixed: `L_mesh=0.033` (so elements-per-inclusion is constant)
- **5 box sizes × 5 random realizations = 25 RVEs**

| L | L/d | L/r | ~N inclusions |
|---|---|---|---|
| 0.25 | 3.1 | 6.3 | 12 |
| 0.35 | 4.4 | 8.8 | 32 |
| 0.45 | 5.6 | 11.3 | 68  ← current production size |
| 0.55 | 6.9 | 13.8 | 124 |
| 0.65 | 8.1 | 16.3 | 205 |

Each RVE runs 3 load cases → **E_eff** (uniaxial), **G_eff** (shear), and the
**MCST length scale `l`** (bending). 25 × 3 = **75 Abaqus solves**.

## Run on CSC

Use the existing pipeline (`01_generate.sh` → solve → post → merge) in a **clean,
dedicated work dir** so it doesn't mix with the main parametric run:

```bash
# 1. dedicated work dir
STUDY=/scratch/project_XXXXXX/rve_size_study
mkdir -p "$STUDY/logs"

# 2. copy the deck + pipeline scripts + analyzer, point the scripts at $STUDY
cp rve_size_study.csv analyze_rve_study.py Spatium_*.py "$STUDY"/
cp 0[1-4]_*.sh submit_all.sh "$STUDY"/
sed -i "s|/scratch/project_XXXXXX/test_rve|$STUDY|g" "$STUDY"/0[1-4]_*.sh "$STUDY"/submit_all.sh
cd "$STUDY"

# 3. pick the study CSV (the dir has only one CSV, but be explicit) + fixed seed
export SPAX_CSV=rve_size_study.csv
export SPAX_SEED=1          # reproducible realizations (one distinct seed per row)

# 4. submit: generate -> solve array (75) -> postprocess (25) -> merge -> results.csv
sbatch 01_generate.sh
```

`01_generate.sh` auto-sizes the solve array to the 75 jobs and the post array to
the 25 RVEs, and chains them with SLURM dependencies. The largest boxes (L=0.65,
~205 inclusions, ~0.75 M elements) are the heavy ones — the 6 h solve walltime is
plenty.

## Analyse

Once `results.csv` is produced:

```bash
python3 analyze_rve_study.py results.csv .
```

Outputs:
- **`rve_convergence.png`** — mean ± 1 std of E_eff / G_eff / l / D_ratio vs `L/d`,
  with the ±2% band around the largest-RVE mean.
- **`rve_cov.png`** — seed-to-seed coefficient of variation vs `L/d`.
- printed table + the recommended representative size per property
  (smallest `L/d` where the mean has plateaued to ±2% **and** CoV < 2%).

The RVE is "statistically fine" at the smallest `L/d` where both the mean stops
moving and the scatter falls under tolerance. Thresholds are tunable at the top
of `analyze_rve_study.py` (`COV_TOL`, `PLATEAU_TOL`).

---

## Extended deck for the length scale `l`

The first study found E_eff/G_eff converge by `L/d ≈ 6–7`, but the MCST length
scale `l` was **not converged** (mean still rising, CoV 15–35% at `L/d=8`).
`rve_size_study_lengthscale.csv` pushes further and doubles the realizations:

- **4 sizes × 10 seeds = 40 RVEs** (120 solves)
- `L/d = 6, 8, 10, 12` (`L = 0.48, 0.64, 0.80, 0.96`), `N ≈ 83, 196, 382, 660`

Run it **exactly like the first study** (`SPAX_CSV=rve_size_study_lengthscale.csv`,
`sbatch 01_generate.sh`), then `analyze_rve_study.py results.csv`.

**This deck is heavy — the big boxes need more resources than the defaults:**

- **Mesh size grows as `(L/d)³`** → the `L=0.96` RVE is ~1.9 M elements. For
  `L ≥ 0.80`, raise the solve allocation in `02_solve.sh`, e.g.
  `--mem-per-cpu=8G` (and consider `--cpus-per-task=20`); the bending case is the
  hungriest. The `L=1.12+` sizes (~3–6 M elements) were deliberately left out as
  they overrun a 40 GB node.
- **Generation is slower** (packing is O(N²); meshing 1–2 M-element periodic RVEs
  takes minutes each). Give `01_generate.sh` more `--time` / `--cpus-per-task` if
  it bumps the 4 h walltime.
- **More mesh retries** are likely at large `N` (more inclusions crossing the
  periodic faces). `SPAX_MAX_RETRIES` defaults to 6; the 10 seeds/size give
  redundancy if a couple of realizations skip.

**Stitch the two studies** for a full `L/d = 3 → 12` curve by concatenating the
result files before analysing (drop the duplicate header):

```bash
cp results.csv results_lengthscale.csv          # the extended run's output
tail -n +2 results_lengthscale.csv >> results_first.csv   # results_first.csv = the L/d 3–8 study
python3 analyze_rve_study.py results_first.csv .
```

If `l` is *still* climbing at `L/d=12`, the converged value is beyond practical
single-RVE FE here — report `l` as an ensemble mean ± standard error, or move to
a larger HPC allocation for `L/d ≥ 14`.

---

## Hybrid approach (recommended): first-order small, bending-only big

`E_eff`/`G_eff` are size-independent (their *mean* converges on small RVEs), and
`l` only needs them through `D_classical = E*/(1-ν*²)·L⁴/12`. So there is no point
solving uniaxial+shear on the expensive big boxes — each load case is a separate
factorisation of the same mesh, so **solving only bending on the big RVEs is ~3×
cheaper**. The length scale is then

    l² = (D_rve − D_classical) / (G*·L²)

with `D_rve` from each big-RVE bending solve and the **converged** `E*, ν*, G*`
taken from the small-RVE first-order study (more correct than a noisy per-RVE
`E_eff` anyway). Validated to reproduce the full per-RVE `l` (e.g. L/d=8.1:
0.108 vs 0.106).

**Stage 1 — first-order (reuse the first study).** The original `rve_size_study.csv`
run already solved uniaxial+shear on small RVEs → its `results.csv` *is* the
first-order source. (Or run a cheap dedicated deck with `Kappa=0`, which skips
bending entirely.) Keep it as e.g. `results_firstorder.csv`.

**Stage 2 — bending-only on the big boxes.** Run `rve_size_study_lengthscale.csv`
but **filter the solve to bending jobs only** — one extra `grep` on the job list
in `01_generate.sh` so the solve array is sized to the 40 `-ben` jobs, not 120:

```bash
# in 01_generate.sh, change the GlobalJobList line to:
ls *.inp 2>/dev/null | sed 's/\.inp$//' | sort | grep -- '-ben$' > GlobalJobList
```

The big RVEs are still *meshed* once (unavoidable), and their `-utx`/`-ss13`
decks are written but never solved. Post-processing then finds only the bending
ODB per RVE → `E_eff`/`G_eff` come out `MISSING` (expected) and **`D_rve` is
written to `results.csv`**. Save it as e.g. `results_bending.csv`.

**Stage 3 — combine.** Pure Python, no Abaqus:

```bash
python3 hybrid_length_scale.py results_firstorder.csv results_bending.csv .
```

Prints `E*, ν*, G*` (with scatter), the per-size `l` mean ± std + CoV, a plateau
verdict, and saves `hybrid_length_scale.png`. `FO_MIN_LD` at the top of the
script sets the minimum `L/d` of first-order RVEs used for the converged moduli.
