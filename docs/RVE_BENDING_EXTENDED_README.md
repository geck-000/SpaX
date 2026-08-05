# Extended bending-only study (push the length scale `l` to convergence)

The hybrid analysis showed `E_eff`/`G_eff` converge on small RVEs but the MCST
length scale `l` is **still climbing** at `L/d=10` (0.077 → 0.105 → 0.143, no
plateau). This deck extends the bending-only side to **`L/d = 12, 14, 16`** so
`l` can be pushed toward its limit, reusing the small-RVE first order via
`hybrid_length_scale.py`.

Bending-only is what makes this affordable: on the big boxes only the bending
load case is solved (~3× cheaper than all three), and the first-order moduli come
from the small-RVE study, not these expensive RVEs.

## Files in this package (copy ALL to the CSC work dir)

| file | role |
|---|---|
| `rve_bending_extended.csv` | the deck: 3 sizes × 8 seeds = **24 RVEs**, `Mode2` empty (bending-only), `Kappa=0.11` |
| `bend_submit.sh` | one-shot submit (login node) |
| `bend_01_generate.sh` | mesh + drop `-utx` decks + build bending job list + submit solve/post/merge |
| `bend_02_solve.sh` | **the heavy stage** — bending solve on big/hugemem nodes |
| `bend_03_postprocess.sh` | extract `D_rve` per RVE (E/G come out `MISSING`, expected) |
| `bend_04_merge.sh` | union partials → `results_bending.csv` |

Also needs (already in the repo on CSC): `Spatium_Standalone.py`,
`Spatium_GmshPeriodic.py`, `Spatium_PostProcess.py`, and — for the final
analysis — `hybrid_length_scale.py` + the small-RVE first-order `results.csv`.

| L | L/d | ~N inclusions | ~elements |
|---|---|---|---|
| 0.96 | 12 | 660 | ~1.9 M |
| 1.12 | 14 | 1048 | ~3.1 M |
| 1.28 | 16 | 1565 | ~4.6 M |

## ⚠️ Feasibility — read before submitting

These are **research-scale** solves. `L/d=12` already **OOM'd at 40 GB**; the
bending case is RAM-bound (Lesicar second-order constraints + direct solver).
Rough memory: `L/d=12` ~150 GB, `L/d=14` ~250 GB, `L/d=16` ~400 GB+.

- **Abaqus scratch on node-local NVMe.** `bend_02` requests `--gres=nvme:1500`
  and points Abaqus `scratch="$LOCAL_SCRATCH"`: the direct solver's out-of-core
  files are huge and must NOT sit on Lustre `/scratch`. hugemem nodes have
  up to ~6 TB NVMe — raise the `nvme` request if a solve runs out of scratch. The
  script falls back to a WORKDIR scratch dir if NVMe wasn't granted.
- `L/d=16` (~4.6 M elements) is the stretch goal — **start with 12 and 14**; if
  16 won't solve even on hugemem, the direct solver can't take it and you'd need
  the iterative solver or a coarser mesh (a separate change).
- Generation (`bend_01`) also meshes 2–5 M-element RVEs — slow and memory-heavy;
  `SPAX_GEN_WORKERS=4` limits concurrency. Bump `--time`/`--mem-per-cpu` if needed.
- 8 seeds/size is a balance for `l`'s ~17 % CoV; trim or raise in the deck to
  taste (cost scales linearly).

## Run

```bash
# 1. dedicated work dir + the code + this package
STUDY=/scratch/project_XXXXXX/rve_bending_extended
mkdir -p "$STUDY/logs"
cp rve_bending_extended.csv bend_*.sh "$STUDY"/
cp Spatium_*.py hybrid_length_scale.py "$STUDY"/

# 2. set WORKDIR/account/partitions inside bend_0*.sh to match your system
#    (they default to /scratch/project_XXXXXX/rve_bending_extended).

# 3. submit
cd "$STUDY"
./bend_submit.sh        # generate -> solve(array) -> post(array) -> merge -> results_bending.csv
```

## Analyse (off the cluster)

Combine with the **small-RVE first-order** `results.csv` (from the original
`rve_size_study.csv` run):

```bash
python3 hybrid_length_scale.py results.csv results_bending.csv .
```

Stitch with the earlier L/d=6–10 bending data for the full `l(L/d)` curve by
concatenating the bending result files first (drop the duplicate header):

```bash
cp results_lscale.csv all_bending.csv
tail -n +2 results_bending.csv >> all_bending.csv
python3 hybrid_length_scale.py results.csv all_bending.csv .
```

`l` is converged once its mean stops rising with `L/d` and the per-size CoV
settles. If it is still climbing at `L/d=16`, report `l` as an ensemble
mean ± standard error and note that the couple-stress limit is beyond reach of
single-RVE FE for this microstructure.
