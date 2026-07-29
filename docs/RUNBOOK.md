# CSC Puhti runbook — L=0.40 second-order + failure-onset study

Budget on hand: **1671 BU**. Both campaigns together cost **well under ~250 BU**
worst-case (realistically ~80–130 BU). No hugemem, no large partition — these are
small jobs. `small` partition throughout; the failure solves could even use `test`.

WORKDIR on Puhti = `/scratch/project_XXXXXX/test_rve` (the shared dir the other
studies used). Account `project_XXXXXX`, `module load abaqus/2025`.

## 0. One-time: copy files to CSC
From the local repository checkout:
```bash
WD=<cluster>:/scratch/project_XXXXXX/test_rve
# scripts + extractor (this bundle)
scp csc_solve_array.sh submit_si2nd_l400.sh postprocess_si2nd_l400.sh \
    submit_failure.sh postprocess_failure.sh failure_extract.py  $WD/
# the engines must already be there; refresh to this version to be safe
scp Spatium_Standalone.py Spatium_PostProcess.py  $WD/
# L=0.40 decks (already generated locally)
scp out_si2nd/Job-SI2_L400_*.inp  $WD/
# column-slice uniaxial decks (already generated locally)
scp out_column/Job-ICE_z*-utx.inp  $WD/
# CSVs needed by post-processing
scp rve_seaice_2nd.csv results_si2nd.csv  $WD/   # results_si2nd has the L240/L320 rows
```

## A. L=0.40 second-order bending (the actual A1 blocker)
12 decks = 4 seeds × {utx, ss13, ben}, quadratic C3D10H, ~78k tets each.
```bash
ssh <cluster>; cd /scratch/project_XXXXXX/test_rve
bash submit_si2nd_l400.sh        # solve array (small, 8 cpu, 32G, 2h, %6) -> post (rows 9-12)
```
Resources: `--partition=small --cpus-per-task=8 --mem=32G --time=02:00:00`. 8 cores
keeps billing == cores; 2 h caps a stuck job. If a `-ben` deck times out, resubmit
that one with `--time=04:00:00` (still fits `small`).

Then pull results back and finalize the 3-size MCST slope **locally** (pandas only):
```bash
scp <cluster>:/scratch/project_XXXXXX/test_rve/results_si2nd.csv ./   # now 12 rows
python3 Spatium_PostProcess.py analyze eq19 results_si2nd.csv
```
Expected: Eq.19 slope stays ≈0 / negative at the 3rd size → **no MCST length scale**
confirmed (closes the paper's second-order section).

## B. Strength / failure-onset mapping (the new study)
10 column slices ICE_z05..z95, uniaxial only, linear C3D4 (~50k elem) — tiny.
```bash
ssh <cluster>; cd /scratch/project_XXXXXX/test_rve
bash submit_failure.sh           # solve array (small, 4 cpu, 12G, 30 min) -> failure_extract
scp <cluster>:/scratch/project_XXXXXX/test_rve/results_failure.csv ./
# run with the SAME python env you use for analyze_studies.py (needs pandas+matplotlib):
python analyze_failure.py results_failure.csv       # -> study_failure.png + first-failure depth
```
`failure_extract.py` reports, per slice, the matrix max-principal SCF and the
Mohr-Coulomb demand (φ default 30°, set `SPAX_MC_PHI_DEG`) as percentiles normalised
by the macro stress. `analyze_failure.py` converts P99 → first-failure macro stress
σ_fail(z) = σ_t/SCF_p99 (tensile) and 2c·cosφ/MCnorm_p99 (MC); the depth with the
lowest σ_fail fails first. The depth **ranking is strength-independent** — σ_t (def
1 MPa) and c (def 0.6 MPa) only set the MPa scale, so the headline ("warm channelled
base cracks first") needs no strength assumption. Expected: minimum at the porous
base (z85/z95), consistent with SCF_p99≈5.6 and 27% of matrix above SCF 2 there.

## Notes
- All decks request `S,E,LE,EVOL`, so `failure_extract.py` gets correct element
  volumes (→ correct macro S11_bar and volume-weighted volfracs).
- Column slices use box size L=0.50 (hard-coded in postprocess_failure.sh).
- ODBs are pruned to deck+ODB by the solver; delete ODBs after post-processing.
- For more scatter on B, regenerate ICE_z* decks with extra seeds and resubmit;
  still only a few BU each.
