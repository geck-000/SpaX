# tensors/ — per-slice 6×6 elasticity tensors (earlier campaigns)

Full stiffness matrices written by

```bash
abaqus python ../SpaX_PostProcess.py elasticity <odb_dir> <out.csv> <L> <run_id>
```

one `elasticity_tensor_<run_id>.csv` per RVE. These are the earlier campaigns —
the `ICE_z*` graded column and the `BTEN_z*` base slices.

The later, manuscript-facing tensor sets live in a subdirectory each, because
each is a single coherent ensemble rather than an accumulation: `column/` (the
ten-slice column), `basetensor_seeds/` (five packings of the warm base at
`L=0.50`) and `bt80/` (the same base at `L=0.80`).

Keep them in subdirectories. The analysers glob the top level of `tensors/`
non-recursively — `analyze_studies.py` averages `elasticity_tensor_BTEN_*.csv`
over *depth slices*, and flattening `basetensor_seeds/` up would silently fold
five same-depth packing replicates into that average.

Aggregate a set into a depth profile with `../analysis/aggregate_coltensor.py`,
or into an ensemble statement with `../analysis/aggregate_basetensor_seeds.py`.
