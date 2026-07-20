# hpc/ — HPC batch scripts (CSC Roihu / Slurm)

Slurm scripts used to run the campaign on CSC Roihu. They are archival examples
of the generate → solve → post-process pipeline at cluster scale; paths and
account/project fields are specific to that environment and should be adapted
before reuse.

| File | Role |
|------|------|
| `submit_firstorder.sh` | Slurm array: solve the first-order (uniaxial/shear) decks. |
| `submit_coltensor.sh` | Solve the full-column 6×6 tensor decks. |
| `submit_brineK.sh` | Solve the brine `K(T)` study [#5] decks. |
| `submit_nlgeom.sh` | Solve the large-deformation study [#8] decks. |
| `submit_tilt.sh` | Solve the channel-tilt study [#6] decks. |
| `postprocess_firstorder.sh` | `abaqus python SpaX_PostProcess.py` over the solved first-order ODBs. |
| `postprocess_coltensor.sh` | Extract the per-slice 6×6 elasticity tensors. |
| `postprocess_nlgeom.sh` | Extract nlgeom reference-point curves. |

**Note.** After the CSC v2026_03 software-stack change, `module load abaqus`
fails in batch; these scripts source a snapshot env (`~/abaqus_env.sh`) instead.
They assume the core modules (`SpaX_Standalone.py`, `SpaX_PostProcess.py`) and
the decks are staged in the job working directory.
